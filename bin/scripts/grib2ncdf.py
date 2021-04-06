"""
This converts a gribfile from grib format to NetCDF format and also creates a
second file interpolated to a regular GPH grid.
The grib file in question is required to contain GPH for the surface and the
logarithm of surface pressure. Full pressure and GPH fields will be constructed.
All contained variables will be used to fill the 4-D field. All fields will be
reinterpolated to specified grid.

Copyright (C) 2019 by Forschungszentrum Juelich GmbH
Author(s): Joern Ungermann
"""

import argparse
import datetime
import juregrid3d
import logging
import netCDF4 as nc
import numpy as np
import os
import pygrib
import sys
import tqdm


LOG = logging.getLogger('juregrid3d')
ECMWF_META = [
    "class", "type", "stream", "typeOfProcessedData",
    "experimentVersionNumber", "gridType", "g2grid"]
VARNAME_TRANSLATE = {"T": "TEMP",
                     "D": "DIVERGENCE"}
STANDARD_NAME_TRANSLATE = {
    "10 metre U wind component": "surface_eastward_wind",
    "10 metre V wind component": "surface_northward_wind",
    "Boundary layer height": "atmosphere_boundary_layer_thickness",
    "Low cloud cover": "low_cloud_area_fraction",
    "Medium cloud cover": "medium_cloud_area_fraction",
    "High cloud cover": "high_cloud_area_fraction",
    "Sea-ice cover": "sea_ice_area_fraction",
    "Specific cloud ice water content": "specific_cloud_ice_water_content",
    "Specific cloud liquid water content": "specific_cloud_liquid_water_content",
    "Fraction of cloud cover": "cloud_area_fraction_in_atmosphere_layer",
}
GPH_FAC = {"m": 1000., "km": 1., "m^2s^-2": 9806.65}
PRESSURE_FAC = {"hPa": 1., "Pa": 100.}
NETCDF_PARAMS = {
    "NETCDF3_64BIT_OFFSET": {"fill_value": np.nan},
    "NETCDF4_CLASSIC": {"fill_value": np.nan},
    "NETCDF4_COMPRESSED": {"zlib": 1, "shuffle": 1, "fletcher32": 1, "fill_value": np.nan}
}


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'gribfile', metavar='filename', type=str)
    parser.add_argument(
        '--format', '-f', metavar='netcdf-format', type=str,
        choices=["NETCDF4_CLASSIC", "NETCDF4_COMPRESSED", "NETCDF3_64BIT_OFFSET"], default="NETCDF4_CLASSIC",
        help="Select NetCDF format (NETCDF4_COMPRESSED is NETCDF4_CLASSIC, "
             "but with compression, shuffle, checksum, etc.)")
    parser.add_argument(
        '--output', '-o', metavar='filename', type=str,
        help="Select a custom path and filename for netCDF output (exclusive with --output-path)")
    parser.add_argument(
        '--output-path', '-p', metavar='path', type=str,
        help="Select a custom directory for netCDF output (exclusive with --output)")
    parser.add_argument(
        '--logfile', metavar='logfilename', type=str,
        help="Use a logfile with increased verbosity compared to commandline output")
    parser.add_argument(
        '--levels', metavar='level,level,level', default="10,20,30", type=str,
        help="List of levels. Example: 10,20,30")
    parser.add_argument(
        '--addgrib', metavar='filename,filename,filename', type=str,
        help="Add further grib files with same gridding to the output. "
             "Specify files separated with ',' without spaces.")
    parser.add_argument(
        '--gphgrib', metavar='filename', type=str,
        help="Add a grib file, from which the GPH of the surface layer shall be retrieved."
             "To be used for forecast files without that field, where the associated"
             "analysis file should be used.")
    parser.add_argument(
        '--pressure-units', metavar='units', type=str, default="hPa",
        choices=['hPa', "Pa"],
        help="Select units to store pressure in. Options: ['hPa', 'Pa']. Default='hPa'")
    parser.add_argument(
        '--gph-units', metavar='units', type=str, default="km",
        choices=['km', 'm', 'm^2s^-2'],
        help="Select units to store gph in. Options: ['km', 'm', 'm^2s^-2']. Default='km'")
    parser.add_argument(
        '--time-units', metavar='units', type=str,
        help="Select units to store time in")
    parser.add_argument(
        '--dimension-order', metavar='order', type=str, default="ecmwf",
        choices=["ecmwf", "clams"],
        help="Choose order of dimensions (ecmwf: time, lev, lat, lon; clams: time, lon, lat, lev)")
    parser.add_argument(
        '--level-type', metavar='type', type=str, default="gph",
        choices=["gph", "pressure", "level", "theta"],
        help="Choose level type to interpolate to (gph, pressure, theta, level)")
    parser.add_argument(
        '--version', action='version',
        version=f'%(prog)s {juregrid3d.__version__} {juregrid3d.__revision__}')
    return parser.parse_args()


def get_press(a, b, ps):
    """
    Compute pressure in hPa at current ECMWF level from a, b, and surface pressure
    """
    return (a + b * ps) / 100.


def get_values(latlon, message):
    olatlon = message.latlons()
    if np.allclose(latlon, olatlon):
        return message.values
    else:
        assert np.allclose(latlon[0], olatlon[0])
        assert np.allclose(latlon[1], olatlon[1] + 180)
        assert latlon[1][0, 0] == 0
        idx = abs(olatlon[1] - 0).argmin()
        olon2 = np.column_stack([olatlon[1][:, idx:], olatlon[1][:, :idx]])
        olon2[olon2 < 0] = olon2[olon2 < 0] + 360
        assert np.allclose(latlon[1], olon2)
        return np.column_stack([message.values[:, idx:], message.values[:, :idx]])


def determineLocalGravity(z, lat):
    """
    Compute gravity according to 1967 Geodetic Reference System formula
    plus free air correction...
    """
    return 9.780318 * (1 + 0.0053024 * (np.sin(np.deg2rad(lat)) ** 2) -
                       0.0000058 * (np.sin(2 * np.deg2rad(lat)) ** 2)) - 3.086e-3 * z


def get_gph(p, t, z0, ps=None, q=None):
    """
    Construct geopotential altitude from pressure and temperature field and
    altitude of first level.
    """
    R = 287.8950831024931  # value of clams program
    g = 9.80665
    zs = np.empty_like(t)
    if q is not None:
        t = t * (1. + 0.607717 * q)
    # distance between ground and first full level.
    Rg = 1e-3 * (R / g)
    zs[0, ...] = z0
    if ps is not None:
        zs[0, ...] -= Rg * np.log(p[0] / ps) * t[0]
    for i in tqdm.tqdm(range(1, len(t)), ascii=True, desc="Computing GPH"):
        # Clams uses an isothermal atmosphere of average temperature between levels
        zs[i, ...] = zs[i - 1] - Rg * np.log(p[i] / p[i - 1]) * (t[i - 1] + t[i]) * 0.5
    return zs


# def get_gph(p, t, ps, z0, q=None):
#     """
#     Construct geopotential altitude from pressure and temperature field and
#     altitude of first level.
#     """
#     R = 287.8950831024931  # value of clams program
#     g = 9.80665
#     zs = np.empty_like(t)
#     zs[0, ...] = z0
#     if q is not None:
#         t *= (1. + 0.607717 * q)
#     phl = (p[1:] + p[:-1]) * 0.5
#     log_ph_p0 = np.empty(list(zs.shape[1:]))
#     log_p1_ph = np.empty(list(zs.shape[1:]))
#     Rg = (R / g) * 1e-3
#     zs[0, ...] = z0 - Rg * np.log(p[0] / ps) * t[0]
#     for i in tqdm.tqdm(range(1, len(t)), ascii=True, desc="Computing GPH"):
#         t0, t1 = t[i - 1], t[i]
#         p0, p1, ph = p[i - 1], p[i], phl[i - 1]
#         log_ph_p0[:] = np.log(ph / p0)
#         log_p1_ph[:] = np.log(p1 / ph)
#         zs[i, ...] = zs[i - 1] - Rg * (log_ph_p0 * t0 + log_p1_ph * t1)
#     return zs


def write(options, meta, latlon, typ, levels, data):
    """
    Write data to a netcdf file
    """
    lat, lon = latlon
    ncname = options.output
    if options.output_path is not None:
        ncname = os.path.join(
            options.output_path,
            f'{os.path.basename(options.gribfile).replace(".grib", "")}_{typ}.nc')
    if ncname is None:
        ncname = f'{options.gribfile.replace(".grib", "")}_{typ}.nc'

    LOG.info("Writing '%s'", ncname)
    with nc.Dataset(ncname, "w", format=options.format.replace("COMPRESSED", "CLASSIC")) as nf:
        for entry in [x for x in meta if x not in ["time"]]:
            setattr(nf, entry, meta[entry])

        history = datetime.datetime.now().isoformat() + ": " + " ".join(sys.argv)
        nf.history = history
        nf.date_created = datetime.datetime.now().isoformat()
        nf.date_modified = datetime.datetime.now().isoformat()
        nf.juregrid3d_version = juregrid3d.__version__
        nf.juregrid3d_revision = juregrid3d.__revision__

        nf.createDimension("time", len(meta["time"][1]))
        ncvar = nf.createVariable("time", "d", ("time"))

        if options.time_units is None:
            ncvar.units = "hours since " + meta["time"][0].isoformat()
        else:
            ncvar.units = options.time_units
        ncvar.standard_name = "time"
        ncvar.long_name = "Time"
        ncvar[:] = nc.date2num(meta["time"][1], ncvar.units)
        if typ == "gph":
            typ = "height"

        nf.createDimension(typ, len(levels))
        if typ == "height":
            ncvar = nf.createVariable(typ, "f4", (typ))
            ncvar.standard_name = "atmosphere_altitude_coordinate"
            ncvar.units = options.gph_units
            ncvar.long_name = "Geopotential Height"
            ncvar[:] = levels
        elif typ == "pressure":
            ncvar = nf.createVariable(typ, "f4", (typ))
            ncvar.units = options.pressure_units
            ncvar.standard_name = "atmosphere_pressure_coordinate"
            ncvar.long_name = "Pressure"
            ncvar[:] = levels
        elif typ == "theta":
            ncvar = nf.createVariable(typ, "f4", (typ))
            ncvar.units = "K"
            ncvar.standard_name = "atmosphere_potential_temperature_coordinate"
            ncvar.long_name = "Potential temperature"
            ncvar[:] = levels
        elif typ == "level":
            ncvar = nf.createVariable(typ, "i2", (typ))
            ncvar.standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
            ncvar[:] = levels
        else:
            raise RuntimeError
        nf.createDimension("lat", lat.shape[0])
        ncvar = nf.createVariable("lat", "f4", ("lat"))
        ncvar[:] = lat[:, 0]
        ncvar.units = "degree_north"
        ncvar.standard_name = "latitude"
        ncvar.long_name = "Latitude"
        nf.createDimension("lon", lat.shape[1])
        ncvar = nf.createVariable("lon", "f4", ("lon"))
        ncvar[:] = lon[0, :]
        ncvar.units = "degree_east"
        ncvar.standard_name = "longitude"
        ncvar.long_name = "Longitude"
        for name in tqdm.tqdm(data, ascii=True, desc="Writing"):
            if len(data[name]["values"].shape) == 4:
                if options.dimension_order == "ecmwf":
                    ncvar = nf.createVariable(
                        data[name]["varname"], "f4", ("time", typ, "lat", "lon"),
                        **NETCDF_PARAMS[options.format])
                    ncvar[:] = data[name]["values"]
                    assert np.all(ncvar[:].shape == data[name]["values"].shape)
                else:
                    ncvar = nf.createVariable(
                        data[name]["varname"], "f4", ("time", "lon", "lat", typ),
                        **NETCDF_PARAMS[options.format])
                    ncvar[:] = data[name]["values"].swapaxes(1, 3)
                    assert np.all(ncvar[:].shape == data[name]["values"].swapaxes(1, 3).shape)
            else:
                if options.dimension_order == "ecmwf":
                    ncvar = nf.createVariable(
                        data[name]["varname"], "f4", ("time", "lat", "lon"),
                        **NETCDF_PARAMS[options.format])
                    ncvar[:] = data[name]["values"]
                    assert np.all(ncvar[:].shape == data[name]["values"].shape)
                else:
                    ncvar = nf.createVariable(
                        data[name]["varname"], "f4", ("time", "lon", "lat"),
                        **NETCDF_PARAMS[options.format])
                    ncvar[:] = data[name]["values"].swapaxes(1, 2)
                    assert np.all(ncvar[:].shape == data[name]["values"].swapaxes(1, 2).shape)
            ncvar.units = data[name]["units"]
            if data[name]["cfname"] != "unknown":
                ncvar.standard_name = data[name]["cfname"]
            elif name in STANDARD_NAME_TRANSLATE:
                ncvar.standard_name = STANDARD_NAME_TRANSLATE[name]
            ncvar.long_name = name


def collect_variables(latlon, data, minlev, unames, names, usteps, steps, gr):
    """
    Retrieve all information of variables 'unames' from names/gr and store it
    into 'data'.
    """
    for name in unames:
        LOG.info("Collecting '%s'", name)
        var = np.empty_like(data["Pressure"]["values"], dtype=np.float32)
        var[:] = np.nan
        units, cfname, varname = "unknown", "unknown", name
        for i in range(1, gr.messages + 1):
            if names[i - 1] != name:
                continue
            m = gr.message(i)
            step = m.step
            if m.level == 0:
                if len(var.shape) == 4:
                    shape = data["Pressure"]["values"].shape
                    var = np.empty((shape[0], shape[2], shape[3]), dtype=np.float32)
                    var[:] = np.nan
                var[usteps.index(step), :, :] = get_values(latlon, m)
            else:
                var[usteps.index(step), m.level - minlev, :, :] = get_values(latlon, m)
            units = m.units
            cfname = m.cfName
            varname = m.cfVarName
        varname = varname.upper()
        varname = VARNAME_TRANSLATE.get(varname, varname)
        if len(var.shape) == 4:
            var = var[:, ::-1, :, :]
        data[name] = {"values": var, "units": units, "cfname": cfname,
                      "varname": varname.upper()}


def read_grib(options):
    """
    Read a grib file, fill in missing fields and return structure containing
    all the data. Eats your memory for breakfast.
    """
    def open_grib(filename):
        handle = pygrib.open(filename)

        names, levels, steps = [], [], []
        for i in range(1, handle.messages + 1):
            m = handle.message(i)
            names.append(m.name)
            levels.append(m.level)
            steps.append(m.step)
        return handle, names, levels, steps

    LOG.info("Opening '%s'", options.gribfile)
    gr, names, levels, steps = open_grib(options.gribfile)

    usteps = sorted(list(set(steps)))

    # unames sind alle variablen
    unames = set(names)
    print(unames, levels)
    unames.remove("Logarithm of surface pressure")
    if not options.gphgrib:
        if "Geopotential" not in unames:
            raise RuntimeError("GRIB does not contain GPH message. Forecast File? Use gphgrib option!")
        unames.remove("Geopotential")

    # surface stuff has level0 0. skip that here!
    minlev, maxlev = max(1, min(levels)), max(levels) + 1
    m = gr.message(names.index("Logarithm of surface pressure") + 1)
    latlon = m.latlons()

    # pv attribute contains hybrid factors on level boundaries. compute values within levels.
    # Halbiere m.pv? pv = hyai + hybi; a = hyam, b = hybm
    print(len(m.values), len(get_values(latlon, m)))
    print(m.topLevel)
    a, b = np.split(m.pv, 2)
    print(a, b)
    # Berechne alle Zwischenpunkte von a und b
    a = 0.5 * (a[1:] + a[:-1])
    b = 0.5 * (b[1:] + b[:-1])
    print(a, b)

    # Speicher alle variablen die beim downloaden angegeben wurden in m? Wieso?
    meta = {}
    for att in ECMWF_META:
        try:
            print(att)
            meta[att] = m[att]
        except RuntimeError:
            LOG.error("Attribute '%s' was not in GRIB file.", att)
    print(meta)
    meta["time"] = datetime.datetime(m.year, m.month, m.day, m.hour, m.minute, m.second)
    meta["time"] = meta["time"], [meta["time"] + datetime.timedelta(hours=_x) for _x in usteps]

    LOG.info("Constructing Pressure")
    press = np.empty([len(usteps), maxlev - minlev] + list(m.values.shape), dtype=np.double)
    print(f"Press: {press}")
    # Berechne pressure für logarithm of surface pressure variable
    # ps = e^m[lat,lon]
    # ToDo: pressure: [timestep, levelidx, lat, lon] = get_press(a[l - 1], b[l - 1], ps)?
    # cdo ml2pl
    for idx, name in enumerate(names):
        if name == "Logarithm of surface pressure":
            m = gr.message(idx + 1)
            timeidx = usteps.index(m.step)
            ps = np.exp(get_values(latlon, m))
            print(ps[:10])
            for l in tqdm.tqdm(list(range(minlev, maxlev)),
                               ascii=True, desc="Computing pressure"):
                press[timeidx, l - minlev, :, :] = get_press(a[l - 1], b[l - 1], ps)
    data = {"Pressure": {"values": press[:, ::-1, :, :] * PRESSURE_FAC[options.pressure_units],
                         "units": options.pressure_units,
                         "cfname": "air_pressure", "varname": "PRESS"}}
    collect_variables(latlon, data, minlev, unames, names, usteps, steps, gr)

    # Geopotential Height same dimension as pressure, [timestep, levelidx, lat, lon]
    gphs = np.empty_like(press, dtype=np.double)
    if options.gphgrib:
        LOG.info("Retrieving GPH from '%s'", options.gphgrib)
        gphgr, gphnames = open_grib(options.gphgrib)[:2]
    else:
        gphgr, gphnames = gr, names
    print(options.gphgrib, gphgr, gphnames)
    
    # tl;dr gphs = get_gph(based_pressure, temperature, geopotential/g, ps/100, specific humidity)
    # cdo gheight
    geopotential_found = False
    for idx, name in enumerate(gphnames):
        if name == "Geopotential":
            LOG.info("Found Geopotential %s %s %s", idx, name, m.typeOfLevel)
            m = gphgr.message(idx + 1)
            try:
                timeidx = usteps.index(m.step)
            except ValueError:
                if len(usteps) == 1:
                    timeidx = 0
                    LOG.warning("Using a non-fitting gphgrib file. Proceed with care!")
                else:
                    raise RuntimeError("gphgrib grib file does not fit to grib data")
            gph_ps = None
            if m.typeOfLevel == 'surface':
                if geopotential_found:
                    continue
                # Happens only once
                gph_ps = ps / 100
            if "Specific humidity" in data:
                LOG.info("Using Q for GPH")
                q = data["Specific humidity"]["values"][timeidx, ...]
            else:
                LOG.warning("No Q for GPH!!")
                q = None
            #get_gph(normed_pressure, normed_temperature, geopotential/g, ps/100, specific humidity)
            gphs[timeidx, ...] = get_gph(
                data["Pressure"]["values"][timeidx, ...] / PRESSURE_FAC[options.pressure_units],
                data["Temperature"]["values"][timeidx, ...],
                get_values(latlon, m) / 9806.65, ps=gph_ps, q=q)
            geopotential_found = True

    if not geopotential_found:
        LOG.critical("Found no surface or hybrid geopotentials!")
        exit(1)

    if options.gphgrib:
        gphgr.close()

    LOG.info("Constructing GPH")
    data["GPH"] = {"values": gphs * GPH_FAC[options.gph_units], "units": options.gph_units,
                   "cfname": "geopotential_height", "varname": "GPH"}

    if options.level_type == "theta":
        temp = data["Temperature"]["values"]
        pres = data["Pressure"]["values"] / PRESSURE_FAC[options.pressure_units]  # -> hPa
        data["THETA"] = {"values": temp * (1000. / pres) ** (287.04 / 1004.64),
                         "units": "K", "cfname": "air_potential_temperature", "varname": "THETA"}

    gr.close()

    if options.addgrib:
        gribs = options.addgrib.split(",")
        for grib in gribs:
            LOG.info("Opening '%s' for additional data", grib)
            gr, names, levels, steps = open_grib(grib)
            unames = [
                _x for _x in set(names)
                if _x not in data and _x not in ["Logarithm of surface pressure", "Geopotential"]]
            collect_variables(data, minlev, unames, names, usteps, steps, gr)
            gr.close()

    return meta, latlon, np.arange(minlev, maxlev)[::-1], data

# Ersetzen
def interpolate(options, latlon, data):
    """
    Reinterpolate vertically(!) in GPH using juregrid3d
    """
    LOG.info("Preparing interpolation")
    levels = np.asarray([float(_x) for _x in options.levels.split(",")])
    key = {"theta": "THETA", "gph": "GPH", "pressure": "Pressure"}[options.level_type]
    vert = data[key]["values"]
    vert = vert.astype(np.double)
    if options.level_type == "theta":  # fix non-monotonous profiles
        err = vert[:, 1:, :, :] < vert[:, :-1, :, :]
        while err.sum() > 0:
            vert[:, 1:, :, :][err] = vert[:, :-1, :, :][err] + 1e-8
            err = vert[:, 1:, :, :] < vert[:, :-1, :, :]

    del data[key]  # remove key so no variable in output file is created
    # Extrapolated values are set to NaN
    if options.level_type == "pressure":
        interpolate = juregrid3d.RectilinearInterpolate4D(np.log(vert), 1, np.log(levels))
    else:
        interpolate = juregrid3d.RectilinearInterpolate4D(vert, 1, levels)
    for name in data:
        LOG.info("Interpolating %s (%s)", name, data[name]["varname"])
        values = data[name]["values"]
        if len(values.shape) == 4:
            data[name]["values"] = interpolate(values)
    return latlon, levels, data


def _main():
    options = parse_options()
    juregrid3d.misc.setup_logging(options.logfile)

    if options.output is not None and options.output_path is not None:
        LOG.fatal("Choose either output OR output-path, not both!")
        sys.exit(1)

    meta, latlon, levels, data = read_grib(options)
    print(meta)
    if options.level_type != "level":
        latlon, levels, data = interpolate(options, latlon, data)
    write(options, meta, latlon, options.level_type, levels, data)
    LOG.info("Done")


if __name__ == "__main__":
    _main()
