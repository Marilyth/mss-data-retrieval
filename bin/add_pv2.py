"""
Copyright (C) 2012 by Forschungszentrum Juelich GmbH
Author(s): Joern Ungermann

Please see docstring of main().
"""
#metpy?
from __future__ import print_function

import datetime
import itertools
import optparse
import os
import sys
from metpy.calc import *
from metpy.units import units
from metpy.calc import thickness_hydrostatic
from metpy.constants import earth_gravity, dry_air_gas_constant
import xarray as xr

import netCDF4
import numpy as np
from scipy import interpolate
import tqdm

# Default arguments for creating NetCDF variables
NETCDF4_PARAMS = {
    "zlib": 1, "shuffle": 1, "fletcher32": 1, "fill_value": np.nan}

# Handle variables which are differently named in varying models here
NAMES = {
    "WACCM": {"TEMP": "T",
              "GPH": ("Z3", 1000.)},
    "MSSZ": {"TEMP": "TEMP",
             "GPH": ("GPH", 1)},
    "MSSL": {"TEMP": "TEMP",
             "GPH": ("GPH", 9.81 * 1000)},
    "MSSP": {"TEMP": "TEMP",
             "GPH": ("GPH", 9.81 * 1000)},
    "MSST": {"TEMP": "TEMP",
             "GPH": ("GPH", 9.81 * 1000)},
    "GWFC": {"TEMP": "TEMP",
             "GPH": ("GPH", 1000.)},
    "FNL": {"TEMP": "TEMP",
            "GPH": ("GPH", 1000.)},
    "MACC": {"TEMP": "temp",
             "GPH": ("GPH", 1000.)},
    "ECMWFP": {"TEMP": "TEMP",
               "GPH": ("GPH", 9.81 * 1000.)},
    "ECMWFZ": {"TEMP": "TEMP",
               "GPH": ("GPH", 9.81 * 1000.)},
    "ECMWFT": {"TEMP": "TEMP",
               "GPH": ("GPH", 9.81 * 1000.)},
    "ECMWFH": {"TEMP": "TEMP",
               "GPH": ("GPH", 9.81 * 1000.)},
}

# Enter suitable dimensions for different models here
DIMS = {
    "WACCM": {"FULL": ("time", "lev", "lat", "lon"),
              "HORIZONTAL": ("time", "lat", "lon")},
    "MSSZ": {"FULL": ("time", "height", "lat", "lon"),
             "HORIZONTAL": ("time", "lat", "lon")},
    "MSSL": {"FULL": ("time", "level", "lat", "lon"),
             "HORIZONTAL": ("time", "lat", "lon")},
    "MSSP": {"FULL": ("time", "pressure", "lat", "lon"),
             "HORIZONTAL": ("time", "lat", "lon")},
    "MSST": {"FULL": ("time", "theta", "lat", "lon"),
             "HORIZONTAL": ("time", "lat", "lon")},
    "GWFC": {"FULL": ("time", "altitude", "latitude", "longitude"),
             "HORIZONTAL": ("time", "latitude", "longitude")},
    "FNL": {"FULL": ("time", "press", "lat", "lon"),
            "HORIZONTAL": ("time", "lat", "lon")},
    "MACC": {"FULL": ("time", "lev", "lat", "lon"),
             "HORIZONTAL": ("time", "lat", "lon")},
    "ECMWFP": {"FULL": ("time", "lon", "lat", "press"),
               "HORIZONTAL": ("time", "lon", "lat")},
    "ECMWFH": {"FULL": ("time", "lon", "lat", "hybrid"),
               "HORIZONTAL": ("time", "lon", "lat")},
    "ECMWFT": {"FULL": ("time", "lon", "lat", "theta"),
               "HORIZONTAL": ("time", "lon", "lat")},
    "ECMWFZ": {"FULL": ("time", "lon", "lat", "height"),
               "HORIZONTAL": ("time", "lon", "lat")},
}
assert NAMES.keys() == DIMS.keys()

VARIABLES = {
    "PRESSURE": ("FULL", "hPa", "air_pressure", "Pressure"),
    "THETA": ("FULL", "K", "air_potential_temperature", "Potential Temperature"),
    "PV": ("FULL", "m^2 K s^-1 kg^-1 10E-6", "ertel_potential_vorticity", "Potential Vorticity"),
    "MOD_PV": ("FULL", "m^2 K s^-1 kg^-1 10E-6", "", "Modified Potential Vorticity"),
    "EQLAT": ("FULL", "degree N", "equivalent_latitude", "Equivalent Latitude"),
    "GPH": ("FULL", "km", "geopotential_height", "Geopotential Altitude"),
    "T_BACKGROUND": ("FULL", "K", "air_temperature_background", "Background Temperature"),
    "T_RESIDUAL": ("FULL", "K", "air_temperature_residual", "Temperature Residual"),
    "W": ("FULL", "m/s", "upward_air_velocity", "Vertical Wind Velocity"),
    "N2": ("FULL", "s^-2", "square_of_brunt_vaisala_frequency_in_air", "N^2"),
    "DIVERGENCE": ("FULL", "s^-1", "divergence_of_wind", "Wind Divergence"),
    "SURFACE_UV": ("HORIZONTAL", "m s^-1", "", "Horizontal Wind Speed at "),
    "SURFACE_PV": ("HORIZONTAL", "m^2 K s^-1 kg^-1 10E-6", "", "Potential Vorticity at "),
    "ALPHA_VEL": ("HORIZONTAL", "m s^- ", "",
                  "Average horizontal wind speed between 100 and 400 hPa (threshold 30 m/s)"),
    "DELTA_V_REL": ("HORIZONTAL", "unitless", "",
                    "relative horizontal wind shear between 200 and 500 hPa (threshold 0.4)"),
    "TROPOPAUSE": ("HORIZONTAL", "km", "tropopause_altitude",
                   "vertical location of first WMO thermal tropopause"),
    "TROPOPAUSE_PRESSURE": ("HORIZONTAL", "Pa", "tropopause_air_pressure",
                            "vertical location of first WMO thermal tropopause"),
    "TROPOPAUSE_THETA": ("HORIZONTAL", "K", "tropopause_air_potential_temperature",
                         "vertical location of first WMO thermal tropopause"),
    "TROPOPAUSE_SECOND": ("HORIZONTAL", "km", "secondary_tropopause_altitude",
                          "vertical location of second WMO thermal tropopause"),
    "TROPOPAUSE_SECOND_PRESSURE": ("HORIZONTAL", "Pa", "secondary_tropopause_air_pressure",
                                   "vertical location of second WMO thermal tropopause"),
    "TROPOPAUSE_SECOND_THETA": ("HORIZONTAL", "K", "secondary_tropopause_air_potential_temperature",
                                "vertical location of second WMO thermal tropopause"),
    "N2_MAX_ABOVE_TROPO": ("HORIZONTAL", "s^-2", "max_of_square_of_brunt_vaisala_frequency_above_tropopause_in_air", ""),
    "N2_MEAN_ABOVE_TROPO": ("HORIZONTAL", "s^-2", "mean_of_square_of_brunt_vaisala_frequency_above_tropopause_in_air", ""),
}


def smooth_polynomially(xs, ys, n_points, degree):
    """
    Smoothens a vector by polynomial fit.

    Args:
        xs: vector of coordinates
        ys: vector of values
        n_points: number of points to use in local polynomial fit
        degree: degree of polynomial to fit with
    Result:
        vector of smoothened values
    """
    n_points = min(n_points, len(xs))
    result = ys.copy()
    n_points_h = n_points / 2
    for i in range(len(ys)):
        imin = max(0, i - n_points_h)
        imax = min(len(ys), imin + n_points)
        imin = imax - n_points
        assert imin >= 0
        rads = ys[imin:imax]
        sel = np.isfinite(rads)
        if sel.sum() >= degree * 2:
            p = np.polyfit(xs[imin:imax][sel], rads[sel], degree)
            result[i] = np.polyval(p, xs[i])
    return result

def get_create_variable(ncin, model, name):
    """
    Either retrieves a variable from NetCDF or creates it,
    in case it is not yet present.
    """
    if name not in ncin.variables:
        if name in VARIABLES:
            dim, units, standard_name, long_name = VARIABLES[name]
        else:
            fields = name.split("_")
            assert fields[1] == "SURFACE"
            dim, units, long_name = VARIABLES["_".join(fields[1:4:2])]
            long_name += fields[2]
        var_id = ncin.createVariable(name, "f4", DIMS[model][dim],
                                     **NETCDF4_PARAMS)
        var_id.units = units
        var_id.long_name = long_name
        if standard_name:
            var_id.standard_name = standard_name
    return ncin.variables[name]


def find_tropopause(alts, temps):
    """
    Identifies position of thermal tropopauses in given altitude/temperature
    profile. Has some issues with inversions, which is circumventyed partly by
    setting seek to False, which is not strictly necessary by WMO definition.

    The thermal definition of the tropopause, WMO, 1957:

    (a) The first tropopause is defined as the lowest level at which the lapse
    rate decreases to 2 degree C/km or less, provided also the average lapse rate
    between this level and all higher levels within 2 km does not exceed 2 degree C/km.

    (b) If above the first tropopause the average lapse rate between any level
    and all higher levels within 1 km exceeds 3 degree C/km, then a second tropopause
    is defined by the same criterion as under (a). This tropopause may be either
    within or above the 1 km layer.
    """
    dtdz_wmo = -2
    zmin = 5
    zmax = 22
    alts = np.asarray(alts)
    temps = np.asarray(temps)
    valid = (~(np.isnan(alts) | np.isnan(temps))) & (alts > 2.0) & (alts < 30.0)
    alts, temps = alts[valid], temps[valid]
    if len(alts) < 3:
        return []
    if alts[0] > alts[1]:  # check for proper order and reverse if necessary
        alts = alts[::-1]
        temps = temps[::-1]

    result = []
    # This differentiation is sufficient as we are looking at average lapse rate
    # with respect to higher levels anyway, so using a more accurate left/right
    # differentiation does not really improve things here.
    lapse_rate = (temps[1:] - temps[:-1]) / (alts[1:] - alts[:-1])
    lapse_alts = (alts[1:] + alts[:-1]) / 2.
    seek = True
    for j in range(1, len(lapse_rate)):
        if not seek and lapse_rate[j] < -3:
            ks = [k for k in range(len(temps)) if lapse_alts[j] <= alts[k] <= lapse_alts[j] + 1.]
            # This way of calculating the average lapse rate is optimal. Don't
            # try to improve. Integrate t'/(z1-z0) numerically (no trapez! do it
            # stupid way) with infinitesimal h. Differentiate numerically using
            # same h. Simplify. Voila. As h can be assumed as small as possible,
            # this is accurate.
            if len(ks) > 1:
                k, ks = ks[0], ks[1:]
                avg_lapse = (temps[ks] - temps[k]) / (alts[ks] - alts[k])
                if all(avg_lapse < -3):
                    seek = True
            else:
                seek = True

        if seek and lapse_rate[j - 1] <= dtdz_wmo < lapse_rate[j] \
                and zmin < lapse_alts[j] < zmax:
            alt = np.interp([dtdz_wmo],
                            lapse_rate[j - 1:j + 1], lapse_alts[j - 1:j + 1])[0]

            ks = [_k for _k in range(len(temps)) if alt <= alts[_k] <= alt + 2.]
            if len(ks) > 1:
                k, ks = ks[0], ks[1:]
                avg_lapse = (temps[ks] - temps[k]) / (alts[ks] - alts[k])
                if all(avg_lapse > dtdz_wmo):
                    result.append(alt)
                    seek = False
            else:
                result.append(alt)
                seek = False
    return result


def swap_axes_write(model, variable):
    """
    Swaps the axis of the variable according to model for writing.
    """
    if model.startswith("ECMWF"):
        return variable[:, ::-1, :, :].swapaxes(1, 3)
    return variable


def swap_axes_read(model, variable):
    """
    Swaps the axis of the variable according to model for reading.
    """
    if model.startswith("ECMWF"):
        result = variable.swapaxes(1, 3)[:, ::-1, :, :]
    else:
        result = variable
    if isinstance(result, np.ma.core.MaskedArray):
        result[result.mask] = np.nan
    return result


def swap_axes_hor(model, variable):
    if model.startswith("ECMWF"):
        return variable.swapaxes(1, 2)
    return variable


def parse_args(args):
    oppa = optparse.OptionParser(usage="""
    add_pv.py

    Adds PV and ancillary quantities to 4D model data given as NetCDF.
    Supported model types are ECMWFP (ECMWF on pressure levels), ECMWFZ
    (JURASSIC ECMWF format on altitude levels), FNL, WACCM.

    Usage: add_pv.py [options] <model type> <netCDF file>

    Example:
    add_pv.py ECMWFP ecmwfr_ana_ml_06072912.nc
    """)

    oppa.add_option('--pressure', '', action='store_true',
                    help="Add PRESSURE field")
    oppa.add_option('--theta', '', action='store_true',
                    help="Add THETA potential temperature field")
    oppa.add_option('--n2', '', action='store_true',
                    help="Add N2 static stability.")
    oppa.add_option('--pv', '', action='store_true',
                    help="Add PV potential vorticity.")
    oppa.add_option('--tropopause', '', action='store_true',
                    help="Add first and second tropopause")
    oppa.add_option('--eqlat', '', action='store_true',
                    help="Add equivalent latitude")
    oppa.add_option('--surface_pressure', '', action='store', type=str,
                    help="Add PV and UV on given hPa surfaces, e.g., 200:300:400.")
    oppa.add_option('--surface_theta', '', action='store', type=str,
                    help="Add PV and UV on given theta surfaces, e.g., 200:300:400.")
    oppa.add_option('--jetstream', '', action='store_true',
                    help="Add scalar fields for jetstream identification")
    oppa.add_option('--temp_background', '', action='store_true',
                    help="Add background Temperature to model data")
    oppa.add_option('--vertical_wind', '-w', action='store_true',
                    help="Add vertical wind (requires OMEGA field)")
    oppa.add_option('--divergence', '', action='store_true',
                    help="Add wind divergence")
    opt, arg = oppa.parse_args(args)

    if len(arg) != 2:
        print(oppa.get_usage())
        exit(1)
    if arg[0] not in NAMES:
        print("Unsupported model type:", arg[0])
        print("Supported types are", list(NAMES.keys()))
        exit(1)
    if not os.path.exists(arg[1]):
        print("Cannot find model data at", arg[1])
        exit(1)
    return opt, arg[0], arg[1]


def get_pressure(ncin, model):
    """
    This function gets a pressure field from a model.
    """
    if "PRESSURE" in ncin.variables:
        press = swap_axes_read(model, ncin.variables["PRESSURE"][:])
    else:
        temp = get_temperature(ncin, model)
        if model == "WACCM":
            hyam = ncin.variables["hyam"][:][np.newaxis, :, np.newaxis, np.newaxis]
            hybm = ncin.variables["hybm"][:][np.newaxis, :, np.newaxis, np.newaxis]
            ps = ncin.variables["PS"][:][:, np.newaxis, :, :]

            press = (hyam * 100000. + hybm * ps) / 100.
            del hyam, hybm, ps
        if model == "MACC":
            hyam = ncin.variables["a"][:][np.newaxis, :, np.newaxis, np.newaxis]
            hybm = ncin.variables["b"][:][np.newaxis, :, np.newaxis, np.newaxis]
            ps = ncin.variables["ps"][:][:, np.newaxis, :, :] / 100.
            p0 = ncin.variables["p0"][0] / 100.
            press = hyam * p0 + hybm * ps
            del hyam, hybm

            if "GPH" not in ncin.variables:
                from pyjurassic.core import construct_gph
                gph = np.zeros(press.shape)
                ps = ps.astype(float)
                for iti, ilo, ila in tqdm.tqdm(
                        itertools.product(range(gph.shape[0]), range(gph.shape[3]), range(gph.shape[2])),
                        total=gph.shape[0] * gph.shape[3] * gph.shape[2], ascii=True):
                    gph[iti, :, ila, ilo] = construct_gph(
                        temp[iti, :, ila, ilo].astype(float),
                        press[iti, :, ila, ilo].astype(float),
                        1013.25)

                get_create_variable(ncin, model, "GPH")[:] = swap_axes_write(model, gph)
            del ps, p0
        elif model in ["ECMWFZ", "ECMWFH", "ECMWFT"]:
            press = swap_axes_read(model, ncin.variables["PRESS"][:])
        elif model == "ECMWFP":
            press = np.zeros(temp.shape)
            press.swapaxes(1, 3)[:] = ncin.variables["press"][:][::-1]
        elif model in ["MSSL", "MSST", "MSSZ"]:
            press = ncin.variables["PRESS"][:]
            assert ncin.variables["PRESS"].units in ["Pa", "hPa"]
            if ncin.variables["PRESS"].units == "Pa":
                press /= 100.
        elif model == "MSSP":
            press = np.zeros(temp.shape)
            press.swapaxes(1, 3)[:] = ncin.variables["pressure"][:]
            assert ncin.variables["pressure"].units in ["Pa", "hPa"]
            if ncin.variables["pressure"].units == "Pa":
                press /= 100.
        else:
            assert model == "FNL"
            press = np.zeros(temp.shape)
            press.swapaxes(1, 3)[:] = ncin.variables["press"][:]

    return press


def get_gph(ncin, model):
    if NAMES[model]["GPH"][0] in ncin.variables:
        gph = swap_axes_read(
            model, ncin.variables[NAMES[model]["GPH"][0]][:] / NAMES[model]["GPH"][1])
    else:
        temp = get_temperature(ncin, model)
        gph = np.zeros_like(temp)
        var = None
        if "altitude" in ncin.variables and len(ncin.variables["altitude"].shape) == 1:
            var = "altitude"
        if "height" in ncin.variables and len(ncin.variables["height"].shape) == 1:
            var = "height"
        assert var is not None
        fac = 1
        if ncin.variables[var].units == "m":
            fac = 1. / 1000.
        if model.startswith("ECMWF"):
            gph[:] = ncin.variables[var][:][np.newaxis, ::-1, np.newaxis, np.newaxis] * fac
        else:
            gph[:] = ncin.variables[var][:][np.newaxis, :, np.newaxis, np.newaxis] * fac
    return gph


def get_temperature(ncin, model):
    try:
        temp = swap_axes_read(model, ncin.variables[NAMES[model]["TEMP"]][:])
    except KeyError:
        temp = swap_axes_read(model, ncin.variables["TEMP_RESIDUAL"][:] +
                              ncin.variables["TEMP_BACKGROUND"][:])
    return temp


def add_pressure(ncin, model):
    """
    This function adds pressure to a model, if not already present as
    "PRESSURE". This is mostly here to more efficiently deal with data present
    on model levels in a way consistent with data presented on other grids.
    """
    if "PRESSURE" not in ncin.variables:
        print("Adding PRESSURE...")
        press = get_pressure(ncin, model)
        get_create_variable(ncin, model, "PRESSURE")[:] = swap_axes_write(model, press)


def get_theta(ncin, model):
    """
    This function computes potential temperature from a model if not already present
    """
    if model in ["MSST"]:
        temp = get_temperature(ncin, model)
        theta = np.zeros_like(temp)
        theta[:] = ncin.variables["theta"][:][np.newaxis, :, np.newaxis, np.newaxis]
    elif "THETA" in ncin.variables:
        theta = swap_axes_read(model, ncin.variables["THETA"][:])
    else:
        temp = get_temperature(ncin, model)
        press = get_pressure(ncin, model)
        theta = temp * (1000. / press) ** (287.04 / 1004.64) # **(2 / 7)
    return theta


def add_theta(ncin, model):
    """
    This function add potential temperature to a model, if not already present.
    """
    if "THETA" not in ncin.variables and model not in ["MSST"]:
        print("Adding THETA...")
        theta = get_theta(ncin, model)
        get_create_variable(ncin, model, "THETA")[:] = swap_axes_write(model, theta)


def add_n2(ncin, model):
    """
    This function adds static stability to a model, if not already present.
    """
    if "N2" not in ncin.variables:
        print("Adding N2...")
        theta = get_theta(ncin, model)
        gph = get_gph(ncin, model)
        n_sq = np.empty_like(theta)
        n_sq[:] = np.nan
        n_sq[:, 1:, :, :] = 1e-3 * (9.81 / theta[:, 1:, :, :]) * \
            np.diff(theta, axis=1) / np.diff(gph, axis=1)
        n_sq[:, 0, :, :] = n_sq[:, 1, :, :]
        get_create_variable(ncin, model, "N2")[:] = swap_axes_write(model, n_sq)
    else:
        n_sq = swap_axes_read(model, ncin.variables["N2"][:])
    if "N2_MAX_ABOVE_TROPO" not in ncin.variables and "TROPOPAUSE" in ncin.variables:
        tropo = swap_axes_hor(model, ncin.variables["TROPOPAUSE"][:])
        max_n_sq = np.empty((gph.shape[0], gph.shape[2], gph.shape[3]))
        max_n_sq[:] = np.nan
        mean_n_sq = max_n_sq.copy()
        for iti, ilo, ila in tqdm.tqdm(
                itertools.product(range(gph.shape[0]), range(gph.shape[3]), range(gph.shape[2])),
                total=gph.shape[0] * gph.shape[3] * gph.shape[2], ascii=True):
            sel = ((tropo[iti, ila, ilo] < gph[iti, :, ila, ilo]) &
                   (gph[iti, :, ila, ilo] < (tropo[iti, ila, ilo] + 2)))
            if sel.sum() > 0:
                n_sq_sel = n_sq[iti, sel, ila, ilo]
                max_n_sq[iti, ila, ilo] = n_sq_sel.max()
                mean_n_sq[iti, ila, ilo] = n_sq_sel.mean()
        get_create_variable(ncin, model, "N2_MAX_ABOVE_TROPO")[:] = swap_axes_hor(model, max_n_sq)
        get_create_variable(ncin, model, "N2_MEAN_ABOVE_TROPO")[:] = swap_axes_hor(model, mean_n_sq)


def add_pv(ncin, model):
    """
    This function calculates isentropic PV and adds it as a variable to a model.
    """
    print("Adding PV...")

    def differentiate(f, dh, axis, coord=False):
        """
        if coord is set, dh is not the grid distance, but the coordinates themselves.
        """
        if axis != 3:
            f = f.swapaxes(axis, 3)
            if len(dh.shape) == 4:
                dh = dh.swapaxes(axis, 3)

        df = np.empty(f.shape)
        df[..., 1:-1] = f[..., 2:] - f[..., :-2]
        if len(dh.shape) == 0 or (len(dh.shape) == 4 and dh.shape[-1] == 1):
            assert not coord
            df[..., 1:-1] /= 2 * dh
        else:
            if coord:
                df[..., 1:-1] /= dh[..., 2:] - dh[..., :-2]
            else:
                df[..., 1:-1] /= 2 * dh[..., 1:-1]
        for i, i1, i2 in [(0, 1, 0), (-1, -1, -2)]:
            if coord:
                df[..., i] = (f[..., i1] - f[..., i2]) / (dh[..., i1] - dh[..., i2])
            else:
                if len(dh.shape) == 0:
                    df[..., i] = (f[..., i1] - f[..., i2]) / dh
                else:
                    df[..., i] = (f[..., i1] - f[..., i2]) / dh[..., i]
        if axis != 3:
            df = df.swapaxes(axis, 3)
        return df

    rearth = 6.356766e6                  # Radius of the Earth (meters)
    angrot = 7.29212e-5                  # Omega of Earth (rad/s)

    latc = ncin.variables["lat"][:]
    lonc = ncin.variables["lon"][:]

    u, v = [swap_axes_read(model, ncin.variables[_x][:]) for _x in ["U", "V"]]
    press = get_pressure(ncin, model)
    theta = get_theta(ncin, model)
    temp = get_temperature(ncin, model)

    # Define finely spaced pressure grid to calculate potential vorticity.
    # p_* variables are with respect to this pressure grid.
    newshape = temp.shape

    coslat = np.cos(np.deg2rad(latc))[np.newaxis, np.newaxis, :, np.newaxis]

    dx = np.deg2rad(lonc[1] - lonc[0]) * rearth * coslat
    dy = np.deg2rad(latc[1] - latc[0]) * rearth

    # PV and vertical theta gradient in pressure

    vort = np.zeros(newshape)

    dvdx = differentiate(v, dx, 3)
    dudy = differentiate(u * coslat, dy * coslat, 2)
    dthetadx = differentiate(theta, dx, 3)
    dthetady = differentiate(theta, dy, 2)

    dvdtheta = differentiate(v, theta, 1, coord=True)
    dvdtheta[np.isnan(dvdtheta)] = 0
    dudtheta = differentiate(u, theta, 1, coord=True)
    dudtheta[np.isnan(dudtheta)] = 0
    dthetadp = differentiate(theta, press, 1, coord=True)

    vort[:] = 2.0 * angrot * np.sin(np.deg2rad(latc))[np.newaxis, np.newaxis, :, np.newaxis]
    vort += dvdx - dthetadx * dvdtheta
    vort -= dudy - dthetady * dudtheta
    vort[:, :, abs(latc) == 90, :] = np.nan

    pv = vort * dthetadp * (-0.01 * 9.81 * 1e6)  # hPa->Pa, g, ?

    get_create_variable(ncin, model, "PV")[:] = swap_axes_write(model, pv)


def add_mod_pv(ncin, model, theta0, epsilon):
    """
    This function adds modified potential vorticity to the model.

    See Mueller and Guenther: A generalized form of lait's modified potential vorticity
    """
    print("Adding MOD_PV...")
    theta = get_theta(ncin, model)
    pv = swap_axes_read(model, ncin.variables["PV"][:])
    mod_pv = pv * ((theta / theta0) ** (-epsilon))

    get_create_variable(ncin, model, "MOD_PV")[:] = swap_axes_write(model, mod_pv)


def add_eqlat(ncin, model):
    print("Adding EQLAT...")
    pv = swap_axes_read(model, ncin.variables["PV"][:])
    theta = get_theta(ncin, model)
    eqlat = np.zeros(pv.shape)

    latc = ncin.variables["lat"][:]
    lonc = ncin.variables["lon"][:]
    if min(latc) > -75 or max(latc) < 75:
        print("WARNING:")
        print("  Not enough latitudes present for this to be a global set.")
        print("  EQLAT may not be meaningful.")

    lats = np.zeros(len(latc) + 1)
    lats[:-1] = latc
    lats[1:] += latc
    lats[1:-1] /= 2
    lats = np.deg2rad(lats)

    area = np.absolute(np.sin(lats[:-1]) - np.sin(lats[1:])) / (2 * len(lonc))
    assert area[0] > 0
    if latc[0] > latc[1]:
        baseareas = (np.sin(np.deg2rad(latc[0])) -
                     np.sin(np.deg2rad(latc))) / 2.
    else:
        baseareas = (np.sin(np.deg2rad(latc[-1])) -
                     np.sin(np.deg2rad(latc)))[::-1] / 2.
        latc = latc[::-1]
    assert(baseareas[1] > baseareas[0])

    if model not in ["MSST"]:
        thetagrid = np.hstack([np.arange(250., 400., 2),
                               np.arange(400., 500., 5.),
                               np.arange(500., 750., 10.),
                               np.arange(750., 1000., 25.),
                               np.arange(1000., 3000., 100.)])
        log_thetagrid = np.log(thetagrid)

        newshape = list(pv.shape)
        newshape[1] = len(thetagrid)

        p_theta = np.zeros(newshape)
        p_theta.swapaxes(1, 3)[:] = thetagrid

        # convert u, v, theta to pressure grid
        theta_pv = np.zeros(newshape)
        lp = np.log(theta[0, :, 0, 0])
        reverse = False
        if lp[0] > lp[-1]:
            theta = theta[:, ::-1]
            pv = pv[:, ::-1]
            reverse = True
        for iti, ilo, ila in tqdm.tqdm(
                itertools.product(range(newshape[0]), range(newshape[3]), range(newshape[2])),
                total=newshape[0] * newshape[3] * newshape[2], ascii=True,
                desc="Interpolation to theta levels:"):
            lp = np.log(theta[iti, :, ila, ilo])
            theta_pv[iti, :, ila, ilo] = np.interp(
                log_thetagrid, lp, pv[iti, :, ila, ilo],
                left=np.nan, right=np.nan)
        print()
    else:
        theta_pv = pv
        newshape = list(pv.shape)
        thetagrid = theta[0, :, 0, 0]

    theta_eqlat = np.zeros(newshape)
    for iti in range(newshape[0]):
        for lev in tqdm.tqdm(range(newshape[1]), desc="Integration", ascii=True):
            areas = np.zeros(len(latc) + 1)
            pv_limits = np.zeros(len(area))
            loc_thpv = theta_pv[iti, lev, :, :]
            loc_lat = np.zeros(loc_thpv.shape, dtype="i8")
            loc_lat.swapaxes(0, 1)[:] = range(len(latc))
            loc_lat = loc_lat.reshape(-1)
            thpv_list = loc_thpv.reshape(-1)
            notnanpv = ~(np.isnan(thpv_list))
            if len(thpv_list[notnanpv]) == 0:
                theta_eqlat[iti, lev, :, :] = np.nan
                continue
            missing_area = area[loc_lat[np.isnan(thpv_list)]].sum()
            areas = baseareas.copy()
            missing_fac = (areas[-1] - missing_area) / areas[-1]
            if missing_fac < 0.99:
                areas *= missing_fac
                print("\nWARNING")
                print("    'Fixing' area due to nan in PV at theta ", thetagrid[lev], end=' ')
                print("by a factor of ", missing_fac)

            minpv, maxpv = thpv_list[notnanpv].min(), thpv_list[notnanpv].max()

            thpv_list = sorted(zip(-thpv_list[notnanpv], loc_lat[notnanpv]))

            aind_lat = np.asarray([x[1] for x in thpv_list], dtype="i8")
            apv = np.asarray([x[0] for x in thpv_list])[:-1]
            cum_areas = np.cumsum(area[aind_lat])[1:]
            if len(cum_areas) >= 2:
                pv_limits = np.interp(areas, cum_areas, apv)

                pv_limits[0], pv_limits[-1] = -maxpv, -minpv
                loc_eqlat = np.interp(-loc_thpv, pv_limits, latc)
                theta_eqlat[iti, lev, :, :] = loc_eqlat
            else:
                print("\nWARNING")
                print("    Filling one level to NaN due to missing PV values")
                theta_eqlat[iti, lev, :, :] = np.nan
    print()

    if model not in ["MSST"]:
        # convert pv back to model grid
        for iti, ilo, ila in tqdm.tqdm(
                itertools.product(range(eqlat.shape[0]), range(eqlat.shape[3]), range(eqlat.shape[2])),
                total=eqlat.shape[0] * eqlat.shape[3] * eqlat.shape[2], ascii=True,
                desc="Interpolation back to model levels:"):
            lp = np.log(theta[iti, :, ila, ilo])
            eqlat[iti, :, ila, ilo] = np.interp(
                lp, log_thetagrid, theta_eqlat[iti, :, ila, ilo],
                left=np.nan, right=np.nan)
        if reverse:
            eqlat = eqlat[:, ::-1]
        print()
    else:
        eqlat = theta_eqlat
    get_create_variable(ncin, model, "EQLAT")[:] = swap_axes_write(model, eqlat)


def add_divergence(ncin, model):
    """
    This function calculates isentropic PV and adds it as a variable to a model.
    """
    print("Adding divergence...")

    def differentiate(f, dh, axis, coord=False):
        """
        if coord is set, dh is not the grid distance, but the coordinates themselves.
        """
        if axis != 3:
            f = f.swapaxes(axis, 3)
            if len(dh.shape) == 4:
                dh = dh.swapaxes(axis, 3)

        df = np.empty(f.shape)
        df[..., 1:-1] = f[..., 2:] - f[..., :-2]
        if len(dh.shape) == 0 or (len(dh.shape) == 4 and dh.shape[-1] == 1):
            assert not coord
            df[..., 1:-1] /= 2 * dh
        else:
            if coord:
                df[..., 1:-1] /= dh[..., 2:] - dh[..., :-2]
            else:
                df[..., 1:-1] /= 2 * dh[..., 1:-1]
        for i, i1, i2 in [(0, 1, 0), (-1, -1, -2)]:
            if coord:
                df[..., i] = (f[..., i1] - f[..., i2]) / (dh[..., i1] - dh[..., i2])
            else:
                if len(dh.shape) == 0:
                    df[..., i] = (f[..., i1] - f[..., i2]) / dh
                else:
                    df[..., i] = (f[..., i1] - f[..., i2]) / dh[..., i]
        if axis != 3:
            df = df.swapaxes(axis, 3)
        return df

    rearth = 6.356766e6                  # Radius of the Earth (meters)

    latc = ncin.variables["lat"][:]
    lonc = ncin.variables["lon"][:]

    u, v = [swap_axes_read(model, ncin.variables[_x][:]) for _x in ["U", "V"]]
    press = get_pressure(ncin, model)
    temp = get_temperature(ncin, model)

    # Define finely spaced pressure grid to calculate potential vorticity.
    # p_* variables are with respect to this pressure grid.
    newshape = temp.shape

    coslat = np.cos(np.deg2rad(latc))[np.newaxis, np.newaxis, :, np.newaxis]

    dx = np.deg2rad(lonc[1] - lonc[0]) * rearth * coslat
    dy = np.deg2rad(latc[1] - latc[0]) * rearth

    # PV and vertical theta gradient in pressure

    div = np.zeros(newshape)

    dudx = differentiate(u, dx, 3)
    dvdy = differentiate(v * coslat, dy * coslat, 2)
    dpdx = differentiate(press, dx, 3)
    dpdy = differentiate(press, dy, 2)

    dvdp = differentiate(v, press, 1, coord=True)
    dvdp[np.isnan(dvdp)] = 0
    dudp = differentiate(u, press, 1, coord=True)
    dudp[np.isnan(dudp)] = 0

    div[:] = dudx - dpdx * dudp
    div[:] += dvdy - dpdy * dvdp
    div[:, :, abs(latc) == 90, :] = np.nan

    get_create_variable(ncin, model, "DIVERGENCE")[:] = swap_axes_write(model, div)


def add_surface(ncin, model, typ, levels):
    """
    This function takes PV and hor. Wind from a model and adds a variable where
    these entities are interpolated on the given horizontal hPa planes.
    """

    if levels is None:
        return
    for p in [int(x) for x in levels.split(":")]:
        print("Adding PV, UV on", typ, "level", p)
        pv = swap_axes_read(model, ncin.variables["PV"][:])
        if typ == "PRESSURE":
            vert = get_pressure(ncin, model)
        elif typ == "THETA":
            vert = get_theta(ncin, model)
        else:
            vert = swap_axes_read(model, ncin.variables[typ][:])
        u = swap_axes_read(model, ncin.variables["U"][:])
        v = swap_axes_read(model, ncin.variables["V"][:])
        pv_surf = np.zeros((pv.shape[0], pv.shape[2], pv.shape[3]))
        uv_surf = np.zeros(pv_surf.shape)
        uv = np.sqrt(u ** 2 + v ** 2)

        if vert[0, 0, 0, 0] < vert[0, -1, 0, 0]:
            order = 1
        else:
            order = -1

        for iti, ilo, ila in tqdm.tqdm(
                itertools.product(range(pv.shape[0]), range(pv.shape[3]), range(pv.shape[2])),
                total=pv.shape[0] * pv.shape[3] * pv.shape[2], ascii=True,
                desc="Interpolation to {} level {}".format(typ, p)):
            uv_surf[iti, ila, ilo] = np.interp(
                [p], vert[iti, ::order, ila, ilo], uv[iti, ::order, ila, ilo],
                left=np.nan, right=np.nan)
            pv_surf[iti, ila, ilo] = np.interp(
                [p], vert[iti, ::order, ila, ilo], pv[iti, ::order, ila, ilo],
                left=np.nan, right=np.nan)

        get_create_variable(ncin, model, "%s_SURFACE_%04d_UV" % (typ, p))[:] = \
            swap_axes_hor(model, uv_surf)
        get_create_variable(ncin, model, "%s_SURFACE_%04d_PV" % (typ, p))[:] = \
            swap_axes_hor(model, pv_surf)


def add_tropopauses(ncin, model):
    """
    Adds first and second thermal WMO tropopause to model. Fill value is -999.
    """
    print("Adding first and second tropopause")

    temp = get_temperature(ncin, model)
    press = np.log(get_pressure(ncin, model))
    gph = get_gph(ncin, model)
    theta = get_theta(ncin, model)

    if gph[0, 1, 0, 0] < gph[0, 0, 0, 0]:
        gph = gph[:, ::-1, :, :]
        press = press[:, ::-1, :, :]
        temp = temp[:, ::-1, :, :]
        theta = theta[:, ::-1, :, :]

    valid = np.isfinite(gph[0, :, 0, 0])
    assert gph[0, valid, 0, 0][1] > gph[0, valid, 0, 0][0]
    assert press[0, valid, 0, 0][1] < press[0, valid, 0, 0][0]

    above_tropo1 = np.empty((gph.shape[0], gph.shape[2], gph.shape[3]))
    above_tropo1[:] = np.nan
    above_tropo2 = above_tropo1.copy()
    above_tropo1_press = above_tropo1.copy()
    above_tropo2_press = above_tropo1.copy()
    above_tropo1_theta = above_tropo1.copy()
    above_tropo2_theta = above_tropo1.copy()
    for iti, ilo, ila in tqdm.tqdm(
            itertools.product(range(gph.shape[0]), range(gph.shape[3]), range(gph.shape[2])),
            total=gph.shape[0] * gph.shape[3] * gph.shape[2], ascii=True):
        tropopauses = find_tropopause(gph[iti, :, ila, ilo], temp[iti, :, ila, ilo])
        tropopauses = [x for x in tropopauses if 5 < x < 22]
        if len(tropopauses) > 0:
            above_tropo1[iti, ila, ilo] = min(tropopauses)
            above_tropo1_press[iti, ila, ilo] = np.interp(
                above_tropo1[iti, ila, ilo], gph[iti, :, ila, ilo], press[iti, :, ila, ilo])
            above_tropo1_theta[iti, ila, ilo] = np.interp(
                above_tropo1[iti, ila, ilo], gph[iti, :, ila, ilo], theta[iti, :, ila, ilo])
            second = [x for x in tropopauses if x > above_tropo1[iti, ila, ilo]]
            if len(second) > 0:
                above_tropo2[iti, ila, ilo] = min(second)
            above_tropo2_press[iti, ila, ilo] = np.interp(
                above_tropo2[iti, ila, ilo], gph[iti, :, ila, ilo], press[iti, :, ila, ilo])
            above_tropo2_theta[iti, ila, ilo] = np.interp(
                above_tropo2[iti, ila, ilo], gph[iti, :, ila, ilo], theta[iti, :, ila, ilo])

    above_tropo1_press = np.exp(above_tropo1_press)
    above_tropo2_press = np.exp(above_tropo2_press)

    get_create_variable(ncin, model, "TROPOPAUSE")[:] = swap_axes_hor(model, above_tropo1)
    get_create_variable(ncin, model, "TROPOPAUSE_SECOND")[:] = swap_axes_hor(model, above_tropo2)
    get_create_variable(ncin, model, "TROPOPAUSE_PRESSURE")[:] = swap_axes_hor(model, above_tropo1_press) * 100
    get_create_variable(ncin, model, "TROPOPAUSE_SECOND_PRESSURE")[:] = swap_axes_hor(model, above_tropo2_press) * 100
    get_create_variable(ncin, model, "TROPOPAUSE_THETA")[:] = swap_axes_hor(model, above_tropo1_theta)
    get_create_variable(ncin, model, "TROPOPAUSE_SECOND_THETA")[:] = swap_axes_hor(model, above_tropo2_theta)


def add_jetstream(press, u, v):
    """
    Adds alpha vel and deltavrel accoring to
    Koch06 - AN EVENT-BASED JET-STREAM CLIMATOLOGY AND TYPOLOGY
    """

    alphavel = np.zeros((press.shape[0], press.shape[2], press.shape[3]))
    deltavrel = np.zeros(alphavel.shape)

    uv = np.sqrt(u ** 2 + v ** 2)

    pmin, pmax = 100., 400.
    for iti, ilo, ila in tqdm.tqdm(
            itertools.product(range(press.shape[0]), range(press.shape[3]), range(press.shape[2])),
            total=press.shape[0] * press.shape[3] * press.shape[2], ascii=True,
            desc="Adding jetstream:"):
        uv_border = np.interp(
            [pmin, pmax, 200, 500], press[iti, :, ila, ilo], uv[iti, :, ila, ilo],
            left=np.nan, right=np.nan)
        press_idc = [k for k in range(press.shape[1])
                     if pmin < press[iti, k, ila, ilo] < pmax]
        press_prof = np.asarray(
            [pmin] + press[iti, press_idc, ila, ilo] + [pmax])
        uv_prof = np.asarray(
            [uv_border[0]] + uv[iti, press_idc, ila, ilo] + [uv_border[1]])

        alphavel[iti, ila, ilo] = np.trapz(uv_prof, x=press_prof) / (pmax - pmin)
        deltavrel[iti, ila, ilo] = (uv_border[2] - uv_border[3]) / uv_border[2]

    get_create_variable(ncin, model, "ALPHA_VEL")[:] = swap_axes_hor(model, alphavel)
    get_create_variable(ncin, model, "DELTA_V_REL")[:] = swap_axes_hor(model, deltavrel)


def add_temp_background(ncin, model):
    """
    Adds temperature background (wavenumbers 0-12). Fill value is -999.
    """
    if model in ["ECMWFZ", "MSSZ"]:
        print("Adding temperature background")
        print("detrending with zonal wavenumber 18.")
        temp = get_temperature(ncin, model)
        temp = np.ma.masked_invalid(temp)
        mask = temp.mask.copy()
        assert ncin.variables["height"].units in ["km", "m"]
        scale = 1
        if ncin.variables["height"].units == "m":
            scale = 1000.
        dalt = abs(ncin.variables["height"][1] - ncin.variables["height"][0]) / scale
        npoints_alt = int(np.round(5.5 / dalt))
        degree_alt = 4
        while npoints_alt < (2 * degree_alt):
            degree_alt -= 1
        print("polynomial smoothing over height with {} points, {} degree.".format(
            npoints_alt, degree_alt))
        dlat = abs(ncin.variables["lat"][1] - ncin.variables["lat"][0])
        npoints_lat = int(np.round(5.1 / dlat))
        degree_lat = 4
        while npoints_lat < (2 * degree_lat):
            degree_lat -= 1
        print("polynomial smoothing over latitude with {} points, {} degree.".format(
            npoints_lat, degree_lat))
        xlon = np.linspace(0, 360, temp.shape[3] + 1)[:-1]
        alt = ncin.variables["height"][:][::-1] / scale
        zminidx = len(alt) - 1
        zmaxidx = 0
        assert alt[zminidx] < alt[zmaxidx]
        while alt[zminidx] < 2:
            zminidx -= 1
        while temp.mask[:, zmaxidx, :, :].sum() > 0:
            zmaxidx += 1
        print(zminidx, zmaxidx)
        for it in range(temp.shape[0]):
            for iz in tqdm.tqdm(range(zmaxidx, zminidx + 1), ascii=True):
                if temp.mask[it, iz, :, :].sum() == 0:
                    continue
                for ilat in range(temp.shape[2]):
                    sel = ~temp.mask[it, iz, ilat, :]
                    sel = np.append(np.append(sel, sel), sel)
                    temp_new = np.append(np.append(temp[it, iz, ilat, :], temp[it, iz, ilat, :]), temp[it, iz, ilat, :])
                    xlon_new = np.append(np.append(xlon - 360., xlon), xlon + 360.)
                    if sel.sum() > 0:
                        f = interpolate.interp1d(xlon_new[sel], temp_new[sel])
                        temp[it, iz, ilat, :] = f(xlon)
                    else:
                        temp[it, iz, ilat, :] = 0

        temp_rfft = np.fft.rfft(temp[:, :, :, :], axis=3)
        temp_rfft[:, :, :, 19:] = 0.
        xlat = np.arange(temp_rfft.shape[2])
        xalt = np.arange(temp_rfft.shape[1])
        for it in range(temp_rfft.shape[0]):
            for wn in tqdm.tqdm(range(19), desc="Smoothing wavenumbers", ascii=True):
                for iz in range(zmaxidx, zminidx + 1):
                    temp_rfft[it, iz, :, wn] = \
                        smooth_polynomially(
                            xlat, temp_rfft[it, iz, :, wn].real, npoints_lat, degree_lat) + \
                        1j * smooth_polynomially(
                            xlat, temp_rfft[it, iz, :, wn].imag, npoints_lat, degree_lat)
                for iy in range(temp_rfft.shape[2]):
                    temp_rfft[it, zmaxidx:zminidx + 1, iy, wn] = \
                        smooth_polynomially(
                            xalt, temp_rfft[it, zmaxidx:zminidx + 1, iy, wn].real, npoints_alt, degree_alt) + \
                        1j * smooth_polynomially(
                            xalt, temp_rfft[it, zmaxidx:zminidx + 1, iy, wn].imag, npoints_alt, degree_alt)
        T_b = np.ones(temp.shape) * np.nan
        T_b[:, :, :, :] = np.fft.irfft(temp_rfft, axis=3, n=temp.shape[3])
        temp_rfft = np.ones((temp.shape[0], temp.shape[1], temp.shape[2], 19)) * np.nan
        temp[mask] = np.nan
        T_b[mask] = np.nan
        T_b = np.ma.masked_invalid(T_b)
        temp = np.ma.masked_invalid(temp)
        get_create_variable(ncin, model, "T_BACKGROUND")[:] = swap_axes_write(model, T_b)
        get_create_variable(ncin, model, "T_RESIDUAL")[:] = swap_axes_write(model, temp - T_b)


def add_vertical_wind(ncin, model):
    """
    Adds vertical wind w. Fill value is -999.
    """
    if model == "ECMWFZ":
        print("Adding vertical wind")
        temp = get_temperature(ncin, model)
        pres = get_pressure(ncin, model) * 100.
        omega = swap_axes_read(model, ncin.variables["OMEGA"][:])
        w = - (omega * 287.058 * temp) / (pres * 9.81)
        get_create_variable(ncin, model, "W")[:] = swap_axes_write(model, w)


def add_metpy(option, filename):
    """
    Adds the variables possible through metpy, this is actually only theta and pv...
    Functions for other variables exist, e.g. n2, but don't exactly work
    """
    with xr.open_dataset(filename) as xin:
        if option.theta or option.pv:
            theta = potential_temperature(xin["PRESS"], xin["TEMP"])
            temperature_from_potential_temperature
            xin["THETA"] = theta
        if option.pv:
            xin["THETA"] = xin["THETA"].metpy.assign_crs(grid_mapping_name='latitude_longitude',
                                                         earth_radius=6.356766e6)
            pv = potential_vorticity_baroclinic(xin["THETA"], xin["PRESS"], xin["U"], xin["V"])
            xin["PV"] = pv
        # if option.n2:
            # This is not working.. geopotential height is not recognised as a length
            # n2 = brunt_vaisala_frequency_squared(xin["GPH"], theta)
            # xin["N2"] = n2


def add_rest(option, model, filename):
    """
    Adds the variables not possible through metpy
    """
    # Open NetCDF file as passed from command line
    with netCDF4.Dataset(filename, "r+") as ncin:

        history = datetime.datetime.now().isoformat() + ":" + " ".join(sys.argv)
        if hasattr(ncin, "history"):
            history += "\n" + ncin.history
        ncin.history = history
        ncin.date_modified = datetime.datetime.now().isoformat()

        if option.eqlat:
            add_eqlat(ncin, model)

        add_surface(ncin, model, "PRESSURE", option.surface_pressure)
        add_surface(ncin, model, "THETA", option.surface_theta)

        if option.n2:
            add_n2(ncin, model)

        if option.tropopause:
            add_tropopauses(ncin, model)

        if option.jetstream:
            add_jetstream(ncin, model)

        if option.temp_background:
            add_temp_background(ncin, model)

        if option.vertical_wind:
            add_vertical_wind(ncin, model)


def main():
    option, model, filename = parse_args(sys.argv[1:])
    add_metpy(option, filename)
    add_rest(option, model, filename)


if __name__ == "__main__":
    main()
