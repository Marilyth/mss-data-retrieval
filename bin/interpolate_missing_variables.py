from metpy.interpolate import interpolate_1d
import xarray as xr
import numpy as np
import sys


def interpolate_vertical(ml_file, inter_file, to_interpolate, levels, new_vertical_axis):
    """
    Linearly interpolate the to_interpolate parameter to the levels of new_vertical_axis
    e.g. interpolate GPH (to_interpolate) to levels [340, 360] (levels) of vertical axis THETA (new_vertical_axis)
    """
    with xr.open_dataset(ml_file) as ml:
        x = np.array(ml[new_vertical_axis].data)
        y = np.array(ml[to_interpolate].data)
        interpolated_data = interpolate_1d(levels, x, y, axis=1)
        attributes = ml[to_interpolate].attrs

    with xr.load_dataset(inter_file) as interpolated:
        interpolated[to_interpolate] = interpolated["q"].copy(data=interpolated_data)
        interpolated[to_interpolate].attrs = attributes
        interpolated.to_netcdf(inter_file)


ml = sys.argv[1]
inter_file = sys.argv[2]
to_interpolate = sys.argv[3]
levels = [float(x) for x in sys.argv[4].split(",")]
vertical_axis = sys.argv[5]
interpolate_vertical(ml, inter_file, to_interpolate, levels, vertical_axis)
