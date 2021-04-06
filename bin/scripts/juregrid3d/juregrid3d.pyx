# distutils: language = c++
#cython: boundscheck=False
#cython: wraparound=False

cimport libcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport pair
import numpy as np
cimport numpy as cnp
import cython

cdef extern from "juregrid3d.hpp":
    cdef cppclass Interpolate3D:
        Interpolate3D(int, double*, double*, double*, int, double*, const string&) except +
        vector[pair[int, double]] interpolate(double, double, double) except +
    cdef void compute_weight_matrix_c(
            double*, size_t, size_t, size_t, size_t,
            size_t,
            double*, size_t,
            float*, int*, int*, size_t, size_t, libcpp.bool) except +


def compute_weight_matrix(coords, idx, levels, allow_extrapolate=False):
    import scipy.sparse
    cdef double[:, :, :, ::1] coords_arr = coords
    cdef double[::1] levels_arr = levels
    cdef int cols
    cdef int rows

    cols = coords.shape[0] * coords.shape[1] * coords.shape[2] * coords.shape[3]
    rows = cols / coords.shape[idx] * len(levels)

    indptr = np.zeros(rows + 1, dtype=np.int32)
    indices = np.zeros(rows * 2, dtype=np.int32)
    data = np.zeros(rows * 2, dtype=np.float32)
    cdef float[::1] data_arr = data
    cdef int[::1] indices_arr = indices
    cdef int[::1] indptr_arr = indptr
    cdef libcpp.bool cpp_allow_extrapolate = allow_extrapolate
    print(allow_extrapolate)
    compute_weight_matrix_c(
        &coords_arr[0, 0, 0, 0],
        coords.shape[0], coords.shape[1], coords.shape[2], coords.shape[3],
        idx,
        &levels_arr[0], len(levels_arr),
        &data_arr[0], &indices_arr[0], &indptr_arr[0], rows, cols, cpp_allow_extrapolate)
    return scipy.sparse.csr_matrix(
            (data, indices, indptr), shape=(rows, cols))


class RectilinearInterpolate4D:
    """
    Class to re-interpolate rectilinear 4-D data along one chosen axis
    (call multiple times for a full regridding across multiple axis).

    Very fast, can deal with NaN within the coordinate array as long as they are
    contiguously at the upper or lower boundaries. Does not extrapolate by default,
    if extrapolation is enabled (not recommended) data is linearily extrapolated from
    the neighbouring two valid data points.
    """
    def __init__(self, coords, idx, levels, allow_extrapolate=False):
        """
        Initializes interpolation. Effectively precomputes the weights and factors.

        If the coords array contains NaN at the upper or lower borders, these values
        are ignored. NaN in the middle of the coords array likely leads to nonsensical
        results.

        Parameters
        ----------
        coords (array)           : 4-D array of current coordinates of the axis to reinterpolate.
        idx (int)                : Index of the axis to interpolate in.
        levels (array)           : New values of reinterpolated axis
        allow_extrapolate (bool) : Indicate whether extrapolation is allowed
        """
        self.levels = levels
        self.idx = idx
        self.A = compute_weight_matrix(coords, idx, levels, allow_extrapolate)
        self.oldshape = coords.shape
        self.newshape = list(coords.shape)
        self.newshape[idx] = len(levels)
        self.newshape = tuple(self.newshape)

    def __call__(self, data):
        """
        Computes interpolation.

        Parameters
        ----------
        data (array): 4-D array of the data to be interpolated.

        Returns
        -------
        array: The reinterpolated 4-D array.
        """
        if data.shape != self.oldshape:
            raise RuntimeError("Shape mismatch: {} != {}".format(data.shape, self.oldshape))
        return np.ma.masked_invalid(
            self.A.dot(data.reshape(-1)).reshape(self.newshape))


cdef class PyInterpolate3D:
    """
    Class to regrid irregular 3-D data. The global data is separated in individual
    patches across the globe like a football, and within each, the data is interpolated
    using the chosen method with a local stereographic projection. The patches are small
    enough such that boundary problems are not apparent.
    """

    cdef Interpolate3D *thisptr

    cdef double[::1] alts_arr
    cdef double[::1] lons_arr
    cdef double[::1] lats_arr

    def __cinit__(self, alts, lons, lats, interp="NEAREST", **kwargs):
        """
        Constructors takes the coordinates of the irregular data: alts, lons, and lats.
        'interp' defines the algorithm of how to interpolate with the 3-D Triangulation.
        """
        cdef double[::1] levels_view

        self.alts_arr, self.lons_arr, self.lats_arr = [
            np.ascontiguousarray(x, dtype=np.double) for x in (alts, lons, lats)]
        assert len(self.alts_arr) == len(self.lons_arr)
        assert len(self.alts_arr) == len(self.lats_arr)
        assert all([x.ndim == 1 for x in (self.alts_arr, self.lons_arr, self.lats_arr)])
        assert interp in ["NATURALNEIGHBOUR", "BARYCENTRIC", "DISTANCE", "NEAREST"], interp

        levels_view = np.ascontiguousarray([kwargs.get("scale", 110.)] +
                list(kwargs.get("levels", [])), dtype=np.double)
        self.thisptr = new Interpolate3D(
            len(self.lons_arr), &self.alts_arr[0], &self.lons_arr[0], &self.lats_arr[0],
            len(levels_view), &levels_view[0], interp)

    def __dealloc__(self):
        del self.thisptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef get_weights(self, double z, double lon, double lat):
        cdef vector[pair[int, double]] indice_weights
        cdef int i
        cdef cnp.ndarray indice
        cdef cnp.ndarray[double] weights

        indice_weights = self.thisptr.interpolate(z, lon, lat)
        indice = np.empty(len(indice_weights), dtype=np.int)
        weights = np.empty(len(indice_weights), dtype=np.double)

        for i in range(indice_weights.size()):
            indice[i] = indice_weights[i].first
            weights[i] = indice_weights[i].second
        return indice, weights

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double interpolate(self, double z, double lon, double lat, double[::1] data):
        cdef double result = 0
        cdef vector[pair[int, double]] indice_weights
        cdef int i, length

        indice_weights = self.thisptr.interpolate(z, lon, lat)
        length = indice_weights.size()

        if length > 0:
            for i in range(length):
                result += data[indice_weights[i].first] * indice_weights[i].second
            return result
        else:
            return np.nan

    def regrid(self, alts, lons, lats, data):
        alts = np.ascontiguousarray(alts, dtype=np.double)
        lons = np.ascontiguousarray(lons, dtype=np.double)
        lats = np.ascontiguousarray(lats, dtype=np.double)

        oldshape = alts.shape
        if len(oldshape) > 1:
            alts = alts.reshape(-1)
            lons = lons.reshape(-1)
            lats = lats.reshape(-1)

        assert (len(alts) == len(lats)) and (len(alts) == len(lons)), (len(alts), len(lons), len(lats))

        if type(data) is list:
            result = self._mass_regrid(alts, lons, lats, data)
            if len(oldshape) > 1:
                result = [x.reshape(oldshape) for x in result]
        else:
            result = self._regrid(alts, lons, lats, data)
            if len(oldshape) > 1:
                result = result.reshape(oldshape)
        return result

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef cnp.ndarray _regrid(self, double[::1] alts, double[::1] lons, double[::1] lats, data):
        cdef double[::1] data_view
        cdef int ilon, ilat, ialt
        cdef cnp.ndarray[double] result

        data_view = np.ascontiguousarray(data, dtype=np.double)
        result = np.empty_like(alts, dtype=np.double)
        assert self.alts_arr.shape[0] == data_view.shape[0]

        for ipt in range(alts.shape[0]):
            result[ipt] = self.interpolate(alts[ipt], lons[ipt], lats[ipt], data_view)
        return result

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef list _mass_regrid(self, double[::1] alts, double[::1] lons, double[::1] lats, datas):
        cdef int ipt, idata, i, ndata = len(datas)
        cdef size_t length
        cdef double[::1] view
        cdef vector[double*] results_ptr
        cdef vector[double*] datas_ptr
        cdef double nan = np.nan

        assert all([self.alts_arr.shape[0] == data.shape[0] for data in datas])

        results = [np.zeros_like(alts) for _ in range(len(datas))]

        for idata in range(ndata):
            view = datas[idata]
            datas_ptr.push_back(&view[0]);
            view = results[idata]
            results_ptr.push_back(&view[0]);

        for ipt in range(alts.shape[0]):
            indice_weights = self.thisptr.interpolate(alts[ipt], lons[ipt], lats[ipt])

            length = indice_weights.size()
            if length > 0:
                for idata in range(ndata):
                    for i in range(length):
                        results_ptr[idata][ipt] += datas_ptr[idata][indice_weights[i].first] * indice_weights[i].second
            else:
                for idata in range(ndata):
                    results_ptr[idata][ipt] = nan

        return results


def regrid(xalts, xlons, xlats, alts, lons, lats, data, **kwargs):
    """
    Regrids an irregularly spaced 3-D dataset onto a regular 3-D grid by means of natural
    neighbour interpolation.

    Parameters
    ----------
    xalts, xlons, xlats : 1-D arrays with the regular grid coordinates
    alts, lons, lats, data: n-D arrays with irregular data and associated coordinates
    algo : string, optional
        Defines the interpolation class to employ

    Returns
    -------
    array : n-D array with regridded data
    """
    intp = PyInterpolate3D(alts, lons, lats, **kwargs)
    return intp.regrid(xalts, xlons, xlats, data)
