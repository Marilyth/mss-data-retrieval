"""
Copyright (C) 2011 by Forschungszentrum Juelich GmbH
Author(s): Joern Ungermann

Please see docstring of main().
"""
from __future__ import print_function

import argparse
import sys
import os
import numpy as np
import netCDF4 as nc
import datetime
import juregrid3d
import pint


UR = pint.UnitRegistry()


def convert_to(value, from_unit, to_unit, default=1.):
    try:
        value_unit = UR.Quantity(value, UR(from_unit))
        result = value_unit.to(to_unit).magnitude
    except pint.UndefinedUnitError:
        result = value * default
    except pint.DimensionalityError:
        if UR(to_unit).to_base_units().units == UR.m:
            try:
                result = (value_unit / UR.Quantity(9.81, "m s^-2")).to(to_unit).magnitude
            except pint.DimensionalityError:
                result = value * default
        else:
            result = value * default
    return result


def swap_axes_read(variable):
    """
    Swaps the axis of the variable according to model for reading.
    """
    if np.ma.is_masked(variable):
        return variable.filled(np.nan)
    else:
        return variable


def get_altitudes(options):
    """
    Provides an default altitude grid. The supplied one seems to be appropriate
    given the hybrid pressure levels of ECMWF and the typical demands of
    airborne instrument retrievals.
    """
    if options.levels is not None:
        return np.asarray(options.levels.split(","))
    else:
        if options.level_type == "theta":
            return np.hstack([np.arange(300, 400, 10),
                              np.arange(400, 800, 40),
                              np.arange(800, 1500, 100),
                              np.arange(1500, 2100, 200),
                              [2100, 2350, 2600, 3000]])
        else:
            return np.hstack([np.arange(0, 20, 0.25),
                              np.arange(20, 40, 2),
                              np.arange(40, 61, 4)])


def reinterpolate(options, source, target, needs_conversion, add_var):
    """
    Interpolate data to new altitude corrdinate

    Parameters
    ----------
    options : [type]
        options as provided by argparse
    source : [type]
        input data
    target : [type]
        reinterpolated data
    needs_conversion : [list]
        which variables shall be interpolated to new height axis
    add_var : [type]
        any variables to be added to the output data

    Raises
    ------
    NotImplementedError
        [description]

    Descriptions
    ------------
    """
    if options.level_type == "theta":
        vcoord = swap_axes_read(source.variables[options.vert_var][:])
        vcoord = convert_to(vcoord, source.variables[options.vert_var].units, options.vert_unit)
        height = target.variables["theta"][:]
    elif options.level_type == "pressure":
        vcoord = swap_axes_read(source.variables[options.vert_var][:])
        vcoord = convert_to(vcoord, source.variables[options.vert_var].units, options.vert_unit)
        vcoord = np.log(vcoord)
        height = np.log(target.variables["pressure"][:])
    else:
        vcoord = swap_axes_read(source.variables[options.vert_var][:])
        vcoord = convert_to(vcoord, source.variables[options.vert_var].units, options.vert_unit)
        height = target.variables["height"][:]

    vcoord = np.ascontiguousarray(vcoord, dtype=np.double)
    height = np.ascontiguousarray(height, dtype=np.double)
    interpolate = juregrid3d.RectilinearInterpolate4D(vcoord, options.vertical_axis, height)
    if "gph" in add_var:
        ovar = target.variables["GPH"]
        new_shape = list(vcoord.shape)
        new_shape[options.vertical_axis] = len(source.variables["height"][:])
        gph = np.zeros(new_shape, dtype=np.double)
        gph[:] = source.variables["height"][:]
        ovar[:] = np.transpose(interpolate(gph), options.transpose)

    for var in needs_conversion:
        ivar = swap_axes_read(source.variables[var][:])
        ovar = target.variables[var]
        ovar[:] = np.transpose(interpolate(ivar), options.transpose)

    if "press" in add_var:
        ovar = target.variables["PRESS"]
        if all([_x in source.variables for _x in ("hyam", "hybm", "PS")]):
            press = (source.variables["hyam"][:][np.newaxis, np.newaxis, np.newaxis, :] * 100000. +
                     source.variables["hybm"][:][np.newaxis, np.newaxis, np.newaxis, :] *
                     source.variables["PS"][:][:, :, :, np.newaxis]) * 1e-2
            transpose = [0, 1, 2]
            transpose.insert(options.vertical_axis, 3)
            press = np.transpose(press, transpose)
            assert (press.shape == vcoord.shape), (press.shape, vcoord.shape)
            ovar[:] = np.transpose(interpolate(press), options.transpose)
        elif all([_x in source.variables for _x in ("a", "b", "ps", "ps0")]):
            hyam = source.variables["a"][:][np.newaxis, np.newaxis, np.newaxis, :]
            hybm = source.variables["b"][:][np.newaxis, np.newaxis, np.newaxis, :]
            ps = source.variables["ps"][:][:, :, :, np.newaxis]
            p0 = source.variables["p0"][0]
            press = (hyam * p0 + hybm * ps) * 1e-2
            transpose = [0, 1, 2]
            transpose.insert(options.vertical_axis, 3)
            press = np.transpose(press, transpose)
            assert (press.shape == vcoord.shape)
            ovar[:] = np.transpose(interpolate(press), options.transpose)
        elif "press" in source.dimensions:
            press = np.empty_like(vcoord)
            press[:] = source.variables["press"][:][np.newaxis, np.newaxis, np.newaxis, :]
            ovar[:] = np.transpose(interpolate(press), options.transpose)
        else:
            raise NotImplementedError


def convert(options):
    with nc.Dataset(options.inpfile) as sou, nc.Dataset(options.outfile, "w", format='NETCDF4_CLASSIC') as tar:
        altitudes = get_altitudes(options)

        # copy attributes
        for att_name in [_x for _x in sou.ncattrs() if _x != "_FillValue"]:
            setattr(tar, att_name, getattr(sou, att_name))

        history = datetime.datetime.now().isoformat() + ": " + " ".join(sys.argv)
        if hasattr(tar, "histor"):
            history += ";\n" + tar.history
        tar.history = history
        tar.date_modified = datetime.datetime.now().isoformat()
        tar.juregrid3d_version = juregrid3d.__version__
        tar.juregrid3d_revision = juregrid3d.__revision__

        if options.level_type == "theta":
            vert_name = "theta"
            vert_standard_name = "atmosphere_potential_temperature_coordinate"
            if not options.vert_unit:
                options.vert_unit = "K"
            if not options.vert_var:
                options.vert_var = "THETA"
        elif options.level_type == "pressure":
            vert_name = "pressure"
            vert_standard_name = "atmosphere_pressure_coordinate"
            if not options.vert_unit:
                options.vert_unit = "hPa"
            if not options.vert_var:
                options.vert_var = "PRESS"
        else:
            vert_name = "height"
            vert_standard_name = "atmosphere_altitude_coordinate"
            if not options.vert_unit:
                options.vert_unit = "km"
            if not options.vert_var:
                options.vert_var = "GPH"

        for var_name, var in sou.variables.items():
            if len(var.dimensions) < 4:
                continue
            old_dims = list(var.dimensions)
            new_dims = list(var.dimensions)
            break
        new_dims[options.vertical_axis] = vert_name
        new_dims = np.asarray(new_dims)[options.transpose].tolist()

        # copy dimensions
        add_var = set()
        for dim_name, dim in sou.dimensions.items():
            if dim_name in ["hybrid", "press", "lev", "level"] and "PRESS" not in sou.variables:
                add_var.add("press")
            if dim_name == "height":
                add_var.add("gph")
            if dim_name == old_dims[options.vertical_axis]:
                continue
            if dim.isunlimited():
                tar.createDimension(dim_name, None)
            else:
                tar.createDimension(dim_name, len(dim))

        tar.createDimension(vert_name, len(altitudes))
        tar.createVariable(vert_name, "f4", (vert_name),
                           zlib=1, shuffle=1, fletcher32=1)[:] = altitudes
        tar[vert_name].standard_name = vert_standard_name
        tar[vert_name].units = options.vert_unit

        if "press" in add_var:
            p = tar.createVariable("PRESS", "f4", new_dims,
                                   zlib=1, shuffle=1, fletcher32=1)
            p.units = "hPa"
        if "gph" in add_var:
            p = tar.createVariable("GPH", "f4", new_dims,
                                   zlib=1, shuffle=1, fletcher32=1)
            p.units = sou.variables["height"].units

        needs_conversion = []
        for var_name, var in sou.variables.items():
            if var_name in ["press", "hybrid"]:
                continue
            if var_name == options.vert_var:
                continue
            if len(var.dimensions) < 4:
                if var_name == old_dims[options.vertical_axis]:
                    continue
                if not all(x in tar.dimensions for x in var.dimensions):
                    continue
                print("Copying ", var_name, var.dimensions, var.dtype)
                new_var = tar.createVariable(var_name, var.dtype, var.dimensions,
                                             zlib=1, shuffle=1, fletcher32=1,
                                             fill_value=np.nan if var.dtype is float else None)
                new_var[:] = var[:]
                # fill variable attributes.
                for att_name in var.ncattrs():
                    if att_name != "_FillValue":
                        setattr(new_var, att_name, getattr(var, att_name))
            else:
                print("Creating ", var_name, new_dims)
                if var_name in tar.variables:
                    continue
                new_var = tar.createVariable(var_name, var.dtype, new_dims,
                                             zlib=1, shuffle=1, fletcher32=1,
                                             fill_value=np.nan)
                for att_name in var.ncattrs():
                    if att_name != "_FillValue":
                        setattr(new_var, att_name, getattr(var, att_name))

                needs_conversion.append(var_name)

        print("Converting hybrid/press to altitude...")
        reinterpolate(options, sou, tar, needs_conversion, add_var)
    print("done.")


def main():
    """
    interpolate_model.py

    Interpolates a Model Atmosphere in netCDF format and hybrid coordinates to a
    regular altitude grid.

    Usage: interpolate_ecmwf.py [options] <netCDF ECMWF file> <target file>

    """

    # build option list
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('inpfile', metavar='<input netCDF file>', type=str)
    parser.add_argument('outfile', metavar='<output netCDF file>', type=str)
    parser.add_argument('--levels', '-l', type=str, metavar="alt1,alt2,alt3",
                        help="Provide an comma-separated altitude list in km, e.g., 0.5,1,1.5")
    parser.add_argument("--vertical_axis", "-v", type=int, default=3,
                        help="Define vertical axis in input")
    parser.add_argument("--transpose", "-t", type=int, nargs=4, default=[0, 1, 2, 3],
                        help="Exchange spatial dimension axes of the output variable")
    parser.add_argument(
        '--level-type', metavar='type', type=str, default="gph",
        choices=["gph", "pressure", "theta"],
        help="Choose level type to interpolate to (gph, pressure, theta)")
    parser.add_argument('--vert-var', action='store', type=str, default=None)
    parser.add_argument('--vert-unit', action='store', type=str, default=None)
    args = parser.parse_args()

    if not os.path.exists(args.inpfile):
        print("ERROR: No file to convert found!")
        exit()
    if os.path.exists(args.outfile):
        print("ERROR: target file exists already!")
        exit()

    convert(args)


if __name__ == "__main__":
    main()
