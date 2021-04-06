import juregrid3d
import netCDF4 as nc
import argparse
import numpy as np


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpfile', metavar='str', nargs=1, type=str)
    parser.add_argument('posfile', metavar='str', nargs=1, type=str)
    parser.add_argument('outfile', metavar='str', nargs=1, type=str)
    parser.add_argument('--interp', default="NEAREST", type=str)
    parser.add_argument('--scale', default=None)
    parser.add_argument('--reduce', default=1, type=int)
    parser.add_argument('--vars', default="CO", type=str)
    return parser.parse_args()


def copy_attr(source, target):
    for attr in source.ncattrs():
        setattr(target, attr, getattr(source, attr))


def regrid(options):
    inf = nc.Dataset(options.inpfile[0])

    vertcoor = inf.exp_VERTCOOR_name
    icoords = {}
    for var in ["LON", "LAT", vertcoor.upper()]:
        icoords[var] = inf.variables[var][:][::options.reduce]

    delta = inf.variables[vertcoor.upper() + "_DELTA"][:]
    minimum = getattr(inf, "exp_POS_" + vertcoor + "_min")
    grid = np.asarray([minimum] + list(np.cumsum(delta)))
    scale = inf.exp_POS_r_high

    idata = []
    for var in options.vars.split(","):
        idata.append((var, np.ascontiguousarray(inf.variables[var][:][::options.reduce], dtype=np.double)))

    pnf = nc.Dataset(options.posfile[0])
    ocoords = {}
    for var in ["LON", "LAT", vertcoor.upper()]:
        ocoords[var] = pnf.variables[var][:]
    pnf.close()

    newdata = juregrid3d.regrid(
        ocoords[vertcoor.upper()], ocoords["LON"], ocoords["LAT"],
        icoords[vertcoor.upper()], icoords["LON"], icoords["LAT"],
        [x[1] for x in idata],
        interp=options.interp, levels=grid, scale=scale)

    onf = nc.Dataset(options.outfile[0], "w", format="NETCDF4_CLASSIC")
    copy_attr(inf, onf)

    onf.createDimension("NPARTS", len(ocoords["LON"]))

    for var in ocoords:
        ncvar = onf.createVariable(var, "f4", ["NPARTS"])
        ncvar[:] = ocoords[var]
        copy_attr(inf.variables[var], ncvar)
    for odata, (var, _) in zip(newdata, idata):
        ncvar = onf.createVariable(var, "f4", ["NPARTS"])
        ncvar[:] = odata
        copy_attr(inf.variables[var], ncvar)

    onf.close()
    inf.close()


def _main():
    options = parse_options()
    regrid(options)


if __name__ == "__main__":
    _main()
