from metpy.interpolate import interpolate_1d
import xarray as xr
import numpy as np
import sys


def interpolate_vertical(ml_file, inter_file, to_interpolate, new_vertical_axis):
    """
    Linearly interpolate the to_interpolate parameter to the levels of new_vertical_axis
    e.g. interpolate GPH (to_interpolate) to levels [340, 360] (levels) of vertical axis THETA (new_vertical_axis)
    """
    with xr.load_dataset(inter_file) as interpolated:
        with xr.open_dataset(ml_file) as ml:
            for variable in to_interpolate:
                x = np.array(ml[new_vertical_axis].data)
                y = np.array(ml[variable].data)
                interpolated_data = interpolate_1d(interpolated["lev"].data, x, y, axis=1)
                attributes = ml[variable].attrs

                interpolated[variable] = interpolated["pres"].copy(data=interpolated_data)
                interpolated[variable].attrs = ml[variable].attrs
                interpolated.to_netcdf(inter_file)


ml = sys.argv[1]
inter_file = sys.argv[2]
to_interpolate = sys.argv[3].split(",")
vertical_axis = sys.argv[4]
interpolate_vertical(ml, inter_file, to_interpolate, vertical_axis)
