import xarray as xr
import numpy as np


def inverse(xin, base_variable, base_value, variable, time, lat, lon):
    """
    Return the value for variable at time/lat/lon/base_variable through interpolation
    """
    slice = xin[base_variable][time, :, lat, lon]
    index = np.abs(slice - base_value).argmin()
    print(index.shape)
    return index.data[0]


def interpolate(xin, base_variable, levels, variables=None):
    """
    1. Calculates the ml where the base_variable reaches the desired levels for each time/lat/lon
    2. Interpolate all variables to the calculated ml for their respective time/lat/lon
    3. Delete all unused ml afterwards
    """
    dataset = xin[base_variable]
    print(dataset.shape)
    for time in range(dataset.shape[0]):
        for x in range(dataset.shape[-1]):
            for y in range(dataset.shape[-2]):
                for level in levels:
                    pass


def main():
    filename = "/home/mayb/Desktop/MSS/data-retrieval/analysis-pipeline/mss/2010-04-01T12.an.ml.nc"
    remodel_to = "THETA"
    levels = [340, 360]
    with xr.open_dataset(filename) as xin:
        test = inverse(xin, remodel_to, 340, "PRESS", 0, 0, 0)
        interpolate(xin, remodel_to, levels)


if __name__ == "__main__":
    main()