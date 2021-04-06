import numpy as np
import juregrid3d


def test_compute_matrix_0():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(3)[:, np.newaxis, np.newaxis, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    A = juregrid3d.compute_weight_matrix(coords, 0, levels)
    a = A.dot(coords.reshape(-1)).reshape(2, 4, 5, 6)
    assert np.allclose(a, levels[:, np.newaxis, np.newaxis, np.newaxis])


def test_compute_matrix_1():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(4)[np.newaxis, :, np.newaxis, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    A = juregrid3d.compute_weight_matrix(coords, 1, levels)
    a = A.dot(coords.reshape(-1)).reshape(3, 2, 5, 6)
    assert np.allclose(a, levels[np.newaxis, :, np.newaxis, np.newaxis])


def test_compute_matrix_2():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(5)[np.newaxis, np.newaxis, :, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    A = juregrid3d.compute_weight_matrix(coords, 2, levels)
    a = A.dot(coords.reshape(-1)).reshape(3, 4, 2, 6)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, :, np.newaxis])


def test_compute_matrix_3():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(6)[np.newaxis, np.newaxis, np.newaxis, :]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    A = juregrid3d.compute_weight_matrix(coords, 3, levels)
    a = A.dot(coords.reshape(-1)).reshape(3, 4, 5, 2)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, np.newaxis, :])


def test_compute_matrix_one_value():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(5)[np.newaxis, np.newaxis, :, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(1)
    levels[:] = 1.1
    A = juregrid3d.compute_weight_matrix(coords, 2, levels)
    a = A.dot(coords.reshape(-1)).reshape(3, 4, 1, 6)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, :, np.newaxis])


def test_compute_matrix_many_values():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(5)[np.newaxis, np.newaxis, :, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.linspace(1, 4, 200)
    A = juregrid3d.compute_weight_matrix(coords, 2, levels)
    a = A.dot(coords.reshape(-1)).reshape(3, 4, 200, 6)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, :, np.newaxis])


def test_compute_matrix_nan_bottom():
    coords = np.zeros((3, 4, 5, 6))
    coords_1d = np.asarray([0, 1.3, 2.2, 3, 3.8, 5])
    coords[:] = coords_1d[np.newaxis, np.newaxis, np.newaxis, :]
    values = np.zeros((3, 4, 5, 6))
    values[:] = np.arange(6)[np.newaxis, np.newaxis, np.newaxis, :]
    coords[0, 0, 0, 0] = np.nan
    coords[0, 0, 1, 0] = np.nan
    coords[0, 0, 1, 1] = np.nan
    coords[0, 0, 2, 5] = np.nan
    coords[0, 0, 3, 4] = np.nan
    coords[0, 0, 3, 5] = np.nan
    coords[0, 0, 4, 0] = np.nan
    coords[0, 0, 4, 5] = np.nan
    levels = np.zeros(3)
    levels[:] = (1.1, 2.1, 4.1)
    A = juregrid3d.compute_weight_matrix(coords, 3, levels, False)
    a = A.dot(values.reshape(-1)).reshape(3, 4, 5, 3)
    ref = np.empty((3, 4, 5, 3))
    ref[:] = np.interp(levels, coords_1d, np.arange(6))
    ref[0, 0, 0, 0] = np.nan
    ref[0, 0, 1, 0] = np.nan
    ref[0, 0, 1, 1] = np.nan
    ref[0, 0, 2, 2] = np.nan
    ref[0, 0, 3, 2] = np.nan
    ref[0, 0, 4, 0] = np.nan
    ref[0, 0, 4, 2] = np.nan
    assert np.allclose(a, ref, equal_nan=True)
    A = juregrid3d.compute_weight_matrix(coords, 3, levels, True)
    a = A.dot(values.reshape(-1)).reshape(3, 4, 5, 3)
    ref = np.empty((3, 4, 5, 3))
    ref[:] = np.interp(levels, coords_1d, np.arange(6))
    ref[0, 0, 0, 0] = 7 / 9
    ref[0, 0, 1, 0] = 5 / 8
    ref[0, 0, 1, 1] = 15 / 8
    ref[0, 0, 2, 2] = 35 / 8
    ref[0, 0, 3, 2] = 35 / 8
    ref[0, 0, 4, 0] = 7 / 9
    ref[0, 0, 4, 2] = 35 / 8
    assert np.allclose(a, ref, equal_nan=True)


def test_frontend_0():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(3)[:, np.newaxis, np.newaxis, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    method = juregrid3d.RectilinearInterpolate4D(coords, 0, levels)
    a = method(coords)
    assert np.allclose(a, levels[:, np.newaxis, np.newaxis, np.newaxis])


def test_frontend_1():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(4)[np.newaxis, :, np.newaxis, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    method = juregrid3d.RectilinearInterpolate4D(coords, 1, levels)
    a = method(coords)
    assert np.allclose(a, levels[np.newaxis, :, np.newaxis, np.newaxis])


def test_frontend_2():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(5)[np.newaxis, np.newaxis, :, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    method = juregrid3d.RectilinearInterpolate4D(coords, 2, levels)
    a = method(coords)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, :, np.newaxis])


def test_frontend_3():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(6)[np.newaxis, np.newaxis, np.newaxis, :]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(2)
    levels[:] = (1.1, 2.1)
    method = juregrid3d.RectilinearInterpolate4D(coords, 3, levels)
    a = method(coords)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, np.newaxis, :])


def test_frontend_single_value():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(5)[np.newaxis, np.newaxis, :, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.zeros(1)
    levels[:] = 1.1
    method = juregrid3d.RectilinearInterpolate4D(coords, 2, levels)
    a = method(coords)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, :, np.newaxis])


def test_frontend_many_values():
    coords = np.zeros((3, 4, 5, 6))
    coords[:] = np.arange(5)[np.newaxis, np.newaxis, :, np.newaxis]
    coords[:] += np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6) / (3 * 4 * 5 * 6.)
    levels = np.linspace(1, 4, 200)
    method = juregrid3d.RectilinearInterpolate4D(coords, 2, levels)
    a = method(coords)
    assert np.allclose(a, levels[np.newaxis, np.newaxis, :, np.newaxis])
