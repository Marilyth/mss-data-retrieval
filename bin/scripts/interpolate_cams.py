from __future__ import print_function
import juregrid3d
import numpy as np
import netCDF4 as nc
import argparse


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpfile', metavar='<input CAMS file>', type=str)
    parser.add_argument('outfile', metavar='<output netCDF file>', type=str)
    parser.add_argument('levels', metavar='<levels>', type=str,
                        help="comma separated list, e.g. 100,200,300")
    parser.add_argument('--theta', action="store_true", default=False,
                        help="use theta levels instead of pressure")
    return parser.parse_args()


def _main():
    options = parse_options()
    inf = nc.Dataset(options.inpfile)
    press = (inf.variables["a"][:][np.newaxis, :, np.newaxis, np.newaxis] *
             inf.variables["p0"][0] +
             inf.variables["b"][:][np.newaxis, :, np.newaxis, np.newaxis] *
             inf.variables["ps"][:][:, np.newaxis, :, :]) * 1e-2
    levels = np.asarray([float(_x) for _x in options.levels.split(",")], dtype=np.double)
    print("Computing interpolation weights")
    if options.theta:
        vc_name = "theta"
        vc_standard_name = "atmosphere_potential_temperature_coordinate"
        vc_units = "K"
        if "THETA" in inf:
            theta = inf.variables["THETA"][:]
        else:
            theta = inf.variables["TEMP"][:] * (1000. / press) ** (287.04 / 1004.64)
        interpolate = juregrid3d.RectilinearInterpolate4D(theta, 1, levels)
        vc_positive = "up"
    else:
        vc_name = "press"
        vc_standard_name = "atmosphere_pressure_coordinate"
        vc_units = "hPa"
        interpolate = juregrid3d.RectilinearInterpolate4D(np.log(press), 1, np.log(levels))
        vc_positive = "down"

    onf = nc.Dataset(options.outfile, "w", format=inf.file_format)

    for attname in inf.ncattrs():
        setattr(onf, attname, getattr(inf, attname))

    for dim in ["lat", "lon", "time"]:
        onf.createDimension(dim, len(inf.dimensions[dim]))

    onf.createDimension(vc_name, len(levels))
    ovar = onf.createVariable(vc_name, "f8", (vc_name))
    ovar[:] = levels
    ovar.units = vc_units
    ovar.standard_name = vc_standard_name
    ovar.positive = vc_positive

    for name, ivar in inf.variables.items():
        if len(ivar.shape) < 4:
            if not all([_x in onf.dimensions for _x in ivar.dimensions]):
                print("Skipping     ", name)
                continue
            print("Copying      ", name)
            ovar = onf.createVariable(name, ivar.dtype, ivar.dimensions)
            ovar[:] = ivar[:]

        else:
            print("Interpolating", name)
            ovar = onf.createVariable(name, ivar.dtype, ("time", vc_name, "lat", "lon"))
            ovar[:] = interpolate(ivar[:])

        for attname in ivar.ncattrs():
            setattr(ovar, attname, getattr(ivar, attname))

    if options.theta:
        ovar = onf.createVariable("pressure", "f8", ("time", vc_name, "lat", "lon"))
        ovar.units = "hPa"
        ovar.standard_name = "air_pressure"
        ovar[:] = np.exp(interpolate(np.log(press)))

    if not any(onf.variables[_x].standard_name == "ertel_potential_vorticity"
               for _x in onf.variables if "standard_name" in onf.variables[_x].ncattrs()):
        ovar = onf.createVariable("PV", "f4", ("time", vc_name, "lat", "lon"))
        ovar.standard_name = "ertel_potential_vorticity"
    if not any(onf.variables[_x].standard_name == "air_potential_temperature"
               for _x in onf.variables if "standard_name" in onf.variables[_x].ncattrs()):
        ovar = onf.createVariable("THETA", "f4", ("time", vc_name, "lat", "lon"))
        ovar.standard_name = "air_potential_temperature"

    if onf.variables["time"][0] != 0:
        datetimes = nc.num2date(onf.variables["time"][:], onf.variables["time"].units)
        new_time_unit = "hours since " + datetimes[0].isoformat()
        onf.variables["time"][:] = nc.date2num(datetimes, new_time_unit)
        onf.variables["time"].units = new_time_unit
    inf.close()
    onf.close()


_main()
