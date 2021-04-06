import juregrid3d
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.basemap
import argparse
import numpy as np
from bisect import bisect


def logspace(start, stop, num=50):
    return np.exp(np.linspace(np.log(start), np.log(stop), num, dtype=np.double)),


ALTITUDE_UNITS = {
    "K": "K",
    "Z": "K",
    "M": "km",
    "P": "hPa"
}

ALTITUDE_VARS = {
    "K": "THETA",
    "Z": "ZETA",
    "M": "GPH",
    "P": "PRESS_INIT"
}

ALTITUDE_RANGES = {
    "K": logspace(100., 2500.),
    "Z": logspace(100., 2500.),
    "M": np.arange(1., 40., 1, dtype=np.double),
    "P": logspace(500., 0.1)
}


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile', metavar='str', nargs=1, type=str)
    parser.add_argument('--altitude', default="K", type=str)
    parser.add_argument('--level', default="380K", type=str)
    parser.add_argument('--entity', default="CO", type=str)
    parser.add_argument('--interp', default="NEAREST", type=str)
    parser.add_argument('--scale', default=None)
    parser.add_argument('--cmap', default="Spectral_r", type=str)
    parser.add_argument('--save', default="", type=str)
    parser.add_argument('--projection', default="cyl", type=str)
    parser.add_argument('--show', dest='show', action='store_true')
    parser.add_argument('--no-show', dest='show', action='store_false')
    parser.set_defaults(show=True)
    parser.add_argument('--lonrange', default="", type=str)
    parser.add_argument('--latrange', default="", type=str)
    parser.add_argument('--altrange', default="", type=str)
    parser.add_argument('--entrange', default="", type=str)
    parser.add_argument('--reduce', default=1, type=int)

    options = parser.parse_args()
    assert options.show or options.save

    return options


def read_data(options):
    nf = nc.Dataset(options.datafile[0])
    data = {"CLAMS": {}, "altitude": ALTITUDE_VARS[options.altitude]}
    for var in ["LON", "LAT", options.entity, data["altitude"]]:
        if var in nf.variables:
            data["CLAMS"][var] = nf.variables[var][:][::options.reduce]
    if "GPH" in data["CLAMS"]:
        data["CLAMS"]["GPH"] /= 9.81 * 1000.

    if options.altitude in ["Z", "K"]:
        delta = nf.variables[data["altitude"] + "_DELTA"][:]
        minimum = getattr(nf, "exp_POS_" + data["altitude"].lower() + "_min")
        data["grid"] = np.asarray([minimum] + list(np.cumsum(delta)))
        data["scale"] = nf.exp_POS_r_high
    nf.close()
    return data


def regrid(options, data):
    lonrange = np.arange(-180., 180.)
    latrange = np.arange(-90., 90.)
    altrange = ALTITUDE_RANGES.get(options.altitude, None)
    if options.lonrange:
        lonrange = np.linspace(*[float(_x) for _x in options.lonrange.split(",")], dtype=np.double)
    if options.latrange:
        latrange = np.linspace(*[float(_x) for _x in options.latrange.split(",")], dtype=np.double)
    if options.altrange:
        altrange = np.linspace(*[float(_x) for _x in options.altrange.split(",")], dtype=np.double)

    if options.level[-1] == options.altitude:
        level = float(options.level[:-1])
        data["lonr"], data["latr"] = np.meshgrid(lonrange, latrange)
        data["altr"] = np.ones_like(data["lonr"]) * level
        # Reduce number of points to speed up computation
        sel = None
        if "grid" in data:
            idx = bisect(data["grid"], level)
            sel = ((data["CLAMS"][data["altitude"]] > data["grid"][max(0, idx - 1)]) &
                   (data["CLAMS"][data["altitude"]] < data["grid"][min(len(data["grid"]) - 1, idx + 1)]))
        elif options.altitude == "M":
            sel = (data["CLAMS"]["GPH"] > level - 5) & (data["CLAMS"]["GPH"] < level + 5)
        elif options.altitude == "T":
            sel = (data["CLAMS"]["THETA"] > level / 3.) & (data["CLAMS"]["THETA"] < level * 3.)
        elif options.altitude == "Z":
            sel = (data["CLAMS"]["ZETA"] > level / 3.) & (data["CLAMS"]["ZETA"] < level * 3.)
        elif options.altitude == "P":
            sel = (data["CLAMS"]["PRESS_INIT"] > level / 3.) & (data["CLAMS"]["PRESS_INIT"] < level * 3.)
        if sel is not None:
            for var in data["CLAMS"]:
                data["CLAMS"][var] = data["CLAMS"][var][sel]
    elif options.level[-1] == "E":
        level = float(options.level[:-1])
        data["latr"], data["altr"] = np.meshgrid(latrange, altrange)
        data["lonr"] = np.ones_like(data["latr"]) * level
    elif options.level[-1] == "N":
        level = float(options.level[:-1])
        data["lonr"], data["altr"] = np.meshgrid(lonrange, altrange)
        data["latr"] = np.ones_like(data["lonr"]) * level
    else:
        assert False, options.level[-1]

    kwargs = {"interp": options.interp}
    if "grid" in data:
        kwargs["levels"] = data["grid"]
        kwargs["scale"] = data["scale"]
    if options.scale is not None:
        kwargs["scale"] = float(options.scale)

    data["regridded"] = juregrid3d.regrid(
        data["altr"], data["lonr"], data["latr"],
        data["CLAMS"][data["altitude"]],
        data["CLAMS"]["LON"],
        data["CLAMS"]["LAT"],
        data["CLAMS"][options.entity], **kwargs)
    return data


def plot_basemap(options, data):
    (lon0, lon1), (lat0, lat1) = [(_x.min(), _x.max()) for _x in (data["lonr"], data["latr"])]
    if options.projection == "cyl":
        return mpl_toolkits.basemap.Basemap(
            fix_aspect=False,
            llcrnrlon=lon0, llcrnrlat=lat0,
            urcrnrlon=lon1, urcrnrlat=lat1,
            resolution='c', projection='cyl')
    elif options.projection == "stereo":
        hor_dist = 800 * 110 * (lon1 - lon0) * np.cos(np.deg2rad(0.5 * (lat0 + lat1)))
        ver_dist = 800 * 110 * (lat1 - lat0)
        return mpl_toolkits.basemap.Basemap(
            width=hor_dist,
            height=ver_dist,
            lat_ts=np.mean([lat0, lat1]),
            lat_0=np.mean([lat0, lat1]),
            lon_0=np.mean([lon0, lon1]),
            resolution='c', projection='stere')
    else:
        assert False, options.projection


def plot(options, data):
    basemap = plt
    vmin, vmax, norm = None, None, None
    if options.level[-1] == options.altitude:
        basemap = plot_basemap(options, data)
        basemap.drawcoastlines(color=[0.7, 0.7, 0.7])
        lon_cover = data["lonr"].max() - data["lonr"].min()
        lat_cover = data["latr"].max() - data["latr"].min()
        griddings = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0]
        lon_grid, lat_grid = [
            griddings[np.argmin([round(abs(_cover / _x - 7)) for _x in griddings])]
            for _cover in (lon_cover, lat_cover)]
        basemap.drawparallels(np.arange(-90., 91., lat_grid), labels=[1, 0, 0, 0])
        basemap.drawmeridians(np.arange(-180., 361., lon_grid), labels=[0, 0, 0, 1])
        x, y = basemap(data["lonr"], data["latr"])
        xlabel, ylabel = "", ""  # longitude (E)", "latitude (N)"
    elif options.level[-1] == "E":
        x, y = data["latr"], data["altr"]
        xlabel, ylabel = "latitude (N)", "altitude ({})".format(ALTITUDE_UNITS[options.altitude])
    elif options.level[-1] == "N":
        x, y = data["lonr"], data["altr"]
        xlabel, ylabel = "longitude (N)", "altitude ({})".format(ALTITUDE_UNITS[options.altitude])
    else:
        assert False, "unknown level {}".format(options.level[-1])
    if options.entrange:
        entrange = np.linspace(*[float(_x) for _x in options.entrange.split(",")], dtype=np.double)
        vmin, vmax = entrange.min(), entrange.max()
        norm = matplotlib.colors.BoundaryNorm(entrange, 255)

    basemap.pcolor(x, y, data["regridded"],
                   cmap=options.cmap, rasterized=True,
                   vmin=vmin, vmax=vmax, norm=norm)
    basemap.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if options.save:
        plt.savefig(options.save)
    if options.show:
        plt.show()


def _main():
    options = parse_options()
    data = read_data(options)
    data = regrid(options, data)
    plot(options, data)


if __name__ == "__main__":
    _main()
