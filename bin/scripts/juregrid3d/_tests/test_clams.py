from __future__ import print_function

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import juregrid3d
import netCDF4 as nc
import time


def _random():
    N = 200000

    alts = np.random.rand(N) * 30
    lons = np.random.rand(N) * 360 - 180
    lats = np.random.rand(N) * 180 - 90

    data = np.sin(4 * np.deg2rad(10 + lons)) * 10 * np.cos(3 * np.deg2rad(lats)) * 10 + alts

    lons_p = np.arange(-180, 180, 1, dtype=np.double)
    lats_p = np.arange(-90, 90, 1, dtype=np.double)

    print(alts.shape, lons.shape, lats.shape)
    intp = juregrid3d.PyInterpolate3D(alts, lons, lats, algo="LEVELSFAST", levels=np.arange(0, 30, 2.))
#    intp = juregrid3d.PyInterpolate3D(alts, lons, lats, algo="PATCHFAST")
    print(1)
    start = time.time()
    regrid3 = intp.regrid([10], lons_p, lats_p, data)[0, :, :]
    stop = time.time()
    print("time:", stop - start)
    plt.pcolor(lons_p, lats_p, regrid3.T)
    plt.colorbar()
    plt.show()


def _test2():
    ncf = nc.Dataset("/share/gloria/1/tacts/model/clams/T11_int/init_T11_12092612.nc")
    alts, lons, lats, data, levels = [ncf.variables[x][:] for x in ["ZETA", "LON", "LAT", "HONO2", "ZETA_GRID"]]

    INTERPOLATIONS = [
        "BASIC", "BASICFAST", "PATCH", "PATCHFAST", "LEVELS", "LEVELSFAST"]
    INTERPOLATIONS = ["LEVELSSPHERE"]
    INTERPOLATIONS = [
        "PATCHDIST", "PLEVELSDIST"]

    point = 87.159090909090907, 0, -55
    print(point)
    for scale in [1]:
        for i, algo in enumerate(INTERPOLATIONS):
            intp = juregrid3d.PyInterpolate3D(
                alts, lons, lats, algo=algo, levels=levels, scale=scale)
            result = intp.interpolate(point[0], point[1], point[2])
            print(algo, scale)
            print(result)
            for i, w in zip(result[0], result[1]):
                print(i, alts[i], lons[i], lats[i], w)
            print()
    point = -87.159090909090907, 0, -55
    print(point)
    for scale in [1]:
        for i, algo in enumerate(INTERPOLATIONS):
            intp = juregrid3d.PyInterpolate3D(
                alts, lons, lats, algo=algo, levels=levels, scale=scale)
            result = intp.interpolate(point[0], point[1], point[2])
            print(algo, scale)
            print(result)
            for i, w in zip(result[0], result[1]):
                print(i, alts[i], lons[i], lats[i], w)
            print()


def _main():
    ncf = nc.Dataset("/share/gloria/1/tacts/model/clams/T11_int/init_T11_12092612.nc")
    alts, lons, lats, data, levels = [ncf.variables[x][:] for x in ["ZETA", "LON", "LAT", "HONO2", "ZETA_GRID"]]

    target = 380
    sel = alts > 0  # (alts > 360) & (alts < 400)
    alts, lons, lats, data = [x[sel] for x in (alts, lons, lats, data)]

    print(sel.sum())
    INTERPOLATIONS = [
        "LEVELSFAST", "LEVELSNEAR", "LEVELSDIST",
        "PLEVELSFAST", "PLEVELSNEAR", "PLEVELSDIST", "PLEVELS",
        "PATCH", "PATCHFAST", "PATCHNEAR", "PATCHDIST",
        "SPHERE", "SPHEREFAST", "SPHERENEAR", "SPHEREDIST",
        "BASIC", "BASICFAST", "BASICNEAR", "BASICDIST"]

    alts_p = np.linspace(alts.min(), alts.max(), 100)
    lons_p = np.arange(-180, 181, 1, dtype=np.double)
    lats_p = np.arange(-90, 91, 1, dtype=np.double)
#    for scale in [0.1, 1, 10, 100, 1000, 10000]:
    for algo in INTERPOLATIONS:
        for scale in [0.1, 1, 10, 100, 1000, 10000]:
            print(algo, scale)
            try:
                start = time.time()
                intp = juregrid3d.PyInterpolate3D(
                    alts, lons, lats, algo=algo, levels=levels, scale=scale)
                setup = time.time()
                gridded1 = intp.regrid([target], lons_p, lats_p, data)
                stop1 = time.time()
                gridded2 = intp.regrid(alts_p, [0], lats_p, data)
                stop2 = time.time()
                gridded3 = intp.regrid(alts_p, lons_p, [-70], data)
                stop3 = time.time()
                print(algo, scale, setup - start, stop1 - setup, stop2 - stop1, stop3 - stop2, stop3 - start)
                inv = np.ma.masked_invalid(gridded2[:, 0, :]).mask
                x, y = np.meshgrid(lats_p, alts_p)
                print(list(zip(x[inv], y[inv])))
                plt.clf()
                plt.pcolormesh(lons_p, lats_p, np.ma.masked_invalid(gridded1[0, :, :]).T,
                               vmin=0, vmax=1.0e-8, cmap=plt.cm.spectral)
                plt.colorbar()
                plt.savefig("Z_{}_{:08.2f}.png".format(algo, scale))
                plt.clf()
                plt.pcolormesh(lats_p, alts_p, np.ma.masked_invalid(gridded2[:, 0, :]),
                               vmin=0, vmax=1.6e-8, cmap=plt.cm.spectral)
                plt.colorbar()
                plt.savefig("X_{}_{:08.2f}.png".format(algo, scale))
                plt.clf()
                plt.pcolormesh(lons_p, alts_p, np.ma.masked_invalid(gridded3[:, :, 0]),
                               vmin=0, vmax=1.4e-8, cmap=plt.cm.spectral)
                plt.colorbar()
                plt.savefig("Y_{}_{:08.2f}.png".format(algo, scale))
            except Exception as e:
                print(e)

            if algo.count("LEVELS") > 0:
                break


def _main2():
    ncf = nc.Dataset("/share/gloria/1/tacts/model/clams/T11_int/init_T11_12092612.nc")
    alts, lons, lats, data, levels = [ncf.variables[x][:] for x in ["ZETA", "LON", "LAT", "HONO2", "ZETA_GRID"]]

    target = 380

    TRIAGS = ["PATCH", "BASIC", "SPHERE", "LEVELS"]
    TRIAGS = ["LEVELS"]
    INTERPS = ["BARYCENTRIC", "DISTANCE", "NEAREST", "NATURALNEIGHBOUR"]
    INTERPS = ["DISTANCE2", "NEAREST2"]
    SCALES = [1, 10, 40, 100, 1000]

    alts_p = np.linspace(alts.min(), alts.max(), 100)
    lons_p = np.arange(-180, 181, 1, dtype=np.double)
    lats_p = np.arange(-90, 91, 1, dtype=np.double)

    for interp in INTERPS:
        algo = "PLEVELS"
        scale = 10
        start = time.time()
        intp = juregrid3d.PyInterpolate3D(
            alts, lons, lats, triag=algo, interp=interp, levels=levels, scale=scale)
        setup = time.time()
        gridded1_ref = intp.regrid([target], lons_p, lats_p, data)
        stop1 = time.time()
        gridded2_ref = intp.regrid(alts_p, [0], lats_p, data)
        stop2 = time.time()
        gridded3_ref = intp.regrid(alts_p, lons_p, [-70], data)
        stop3 = time.time()
        print("MEAS", algo, interp, scale, setup - start, stop1 - setup, stop2 - stop1, stop3 - stop2, stop3 - start)
        plt.clf()
        plt.pcolormesh(lons_p, lats_p, np.ma.masked_invalid(gridded1_ref[0, :, :]).T,
                       vmin=0, vmax=1.0e-8, cmap=plt.cm.spectral)
        plt.colorbar()
        plt.savefig("Z_{}_{}_{:08.2f}.png".format(algo, interp, scale))
        plt.clf()
        plt.pcolormesh(lats_p, alts_p, np.ma.masked_invalid(gridded2_ref[:, 0, :]),
                       vmin=0, vmax=1.6e-8, cmap=plt.cm.spectral)
        plt.colorbar()
        plt.savefig("X_{}_{}_{:08.2f}.png".format(algo, interp, scale))
        plt.clf()
        plt.pcolormesh(lons_p, alts_p, np.ma.masked_invalid(gridded3_ref[:, :, 0]),
                       vmin=0, vmax=1.4e-8, cmap=plt.cm.spectral)
        plt.colorbar()
        plt.savefig("Y_{}_{}_{:08.2f}.png".format(algo, interp, scale))
        for algo in TRIAGS:
            for scale in SCALES:

                print(algo, scale)
                try:
                    start = time.time()
                    intp = juregrid3d.PyInterpolate3D(
                        alts, lons, lats, triag=algo, interp=interp, levels=levels, scale=scale)
                    setup = time.time()
                    gridded1 = intp.regrid([target], lons_p, lats_p, data)
                    stop1 = time.time()
                    gridded2 = intp.regrid(alts_p, [0], lats_p, data)
                    stop2 = time.time()
                    gridded3 = intp.regrid(alts_p, lons_p, [-70], data)
                    stop3 = time.time()
                    inv = np.ma.masked_invalid(gridded2[:, 0, :]).mask
                    x, y = np.meshgrid(lats_p, alts_p)
                    print(list(zip(x[inv], y[inv])))
                    plt.clf()
                    plt.pcolormesh(lons_p, lats_p, np.ma.masked_invalid(gridded1[0, :, :]).T,
                                   vmin=0, vmax=1.0e-8, cmap=plt.cm.spectral)
                    plt.colorbar()
                    plt.savefig("Z_{}_{}_{:08.2f}.png".format(algo, interp, scale))
                    plt.clf()
                    plt.pcolormesh(lats_p, alts_p, np.ma.masked_invalid(gridded2[:, 0, :]),
                                   vmin=0, vmax=1.6e-8, cmap=plt.cm.spectral)
                    plt.colorbar()
                    plt.savefig("X_{}_{}_{:08.2f}.png".format(algo, interp, scale))
                    plt.clf()
                    plt.pcolormesh(lons_p, alts_p, np.ma.masked_invalid(gridded3[:, :, 0]),
                                   vmin=0, vmax=1.4e-8, cmap=plt.cm.spectral)
                    plt.colorbar()
                    plt.savefig("Y_{}_{}_{:08.2f}.png".format(algo, interp, scale))
                    plt.clf()

                    plt.pcolormesh(lons_p, lats_p, np.ma.masked_invalid((gridded1 - gridded1_ref)[0, :, :]).T,
                                   vmin=-4e-9, vmax=4e-9,
                                   cmap=plt.cm.RdBu)
                    plt.colorbar()
                    plt.savefig("Z_DELTA_{}_{}_{:08.2f}.png".format(algo, interp, scale))
                    plt.clf()
                    plt.pcolormesh(lats_p, alts_p, np.ma.masked_invalid((gridded2 - gridded2_ref)[:, 0, :]),
                                   vmin=-4e-9, vmax=4e-9,
                                   cmap=plt.cm.RdBu)
                    plt.colorbar()
                    plt.savefig("X_DELTA_{}_{}_{:08.2f}.png".format(algo, interp, scale))
                    plt.clf()
                    plt.pcolormesh(lons_p, alts_p, np.ma.masked_invalid((gridded3 - gridded3_ref)[:, :, 0]),
                                   vmin=-4e-9, vmax=4e-9,
                                   cmap=plt.cm.RdBu)
                    plt.colorbar()
                    plt.savefig("Y_DELTA_{}_{}_{:08.2f}.png".format(algo, interp, scale))
                except Exception as e:
                    print(e)

                if algo.count("LEVELS") > 0:
                    break


def _basic():
    N = 20000

    alts = np.random.rand(N) * 30
    lons = np.random.rand(N) * 360 - 180
    lats = np.random.rand(N) * 180 - 90

    INTERPOLATIONS = [
        "NATURALNEIGHBOUR", "DISTANCE", "BARYCENTRIC", "NEAREST"]
    POINTS = [
        (10, 0, 0),
        (1.5, 0, 0),
        (0.1, 0, 0),
        # (-1, 0, 0),
        (28.1, 0, 0),
        (29.9, 0, 0),
        # (30.9, 0, 0),
        (10, 180, 0),
        (10, -180, 0),
        (10, 90, 0),
        (10, -90, 0),
        (10, 190, 0),
        (10, -190, 0),
        (10, 0, 95),
        (10, 0, -95),
    ]

    levels = np.arange(0., 31., 2.)
    scale = 110
    for algo in INTERPOLATIONS:
        intp = juregrid3d.PyInterpolate3D(
            alts, lons, lats,
            interp=algo,
            levels=levels, scale=scale)
        for point in POINTS:
            print(algo, point)
            result = intp.get_weights(point[0], point[1], point[2])
            print(result)
            for i, w in zip(result[0], result[1]):
                print(alts[i], lons[i], lats[i], w)
            print()


def _basic2():
    N = 20000

    alts = np.random.rand(N) * 300 + 000
    lons = np.random.rand(N) * 360 - 180
    lats = np.random.rand(N) * 180 - 90

    levels = []  # np.asarray([100, 110, 130, 160, 200, 250, 310, 380, 460], dtype=np.double)
    print("01")
    POINTS = [
        (100, 0, 0),
        (150, 0, 0),
        (0.1, 0, 0),
        (-1, 0, 0),
        (28.1, 0, 0),
        (29.9, 0, 0),
        (31, 0, 0),
        (10, 180, 0),
        (10, -180, 0),
        (10, 90, 0),
        (10, -90, 0),
        (10, 190, 0),
        (10, -190, 0),
        (10, 0, 95),
        (10, 0, -95),
    ]

    scale = 100
    intp = juregrid3d.PyInterpolate3D(
        alts, lons, lats,
        interp="NEAREST", levels=levels, scale=scale)
    for point in POINTS:
        print(point)
        result = intp.get_weights(point[0], point[1], point[2])
        print(result)
        for i, w in zip(result[0], result[1]):
            print(alts[i], lons[i], lats[i], w)
        print()


def _basic3():
    alts = np.asarray([0., 0, 0, 1])
    lons = np.asarray([0., 1, 0, 0])
    lats = np.asarray([0., 0, 1, 0])

    POINTS = [
        (0.1, 0.1, 0.1),
        (1, 0, 0),
        (0, 0, -2),
        (1, 2, 1),
        (0, 0.6, -0.1)
    ]

    scale = 110
    print(alts.shape)
    intp = juregrid3d.PyInterpolate3D(
        alts, lons, lats,
        interp="NATURALNEIGHBOUR", scale=scale)
    for point in POINTS:
        print(point)
        result = intp.get_weights(point[0], point[1], point[2])
        print(result)
        for i, w in zip(result[0], result[1]):
            print(alts[i], lons[i], lats[i], w)
        print()

# _basic3()
# _test2()
# _main2()
# _random()
