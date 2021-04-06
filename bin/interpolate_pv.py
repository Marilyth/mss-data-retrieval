import matplotlib.pyplot as plt
import numpy as np
import tqdm
import itertools
import netCDF4 as nc
import sys
import cython


comp = cython.inline("""
def comp(double[::1] pvlevels,
         float[:, :, :, ::1] pv,
         float[:, :, :, ::1] gph,
         float[:, :, :, ::1] pre,
         float[:, :, :, ::1] the,
         float[:, :, :, ::1] gph2,
         float[:, :, :, ::1] pre2,
         float[:, :, :, ::1] the2
         ):
    cdef:
        int it, ipv, ila, ilo, idx

    for it in range(pv.shape[0]):
        for ila in range(pv.shape[2]):
            for ilo in range(pv.shape[3]):
                for ipv in range(pvlevels.shape[0]):
                    for idx in range(1, pv.shape[1]):
                        if gph[it, idx, ila, ilo] < 4 * 9810:
                            continue
                        if pv[it, idx - 1, ila, ilo] < pvlevels[ipv] < pv[it, idx, ila, ilo]:
                            w1 = (pv[it, idx - 1, ila, ilo] - pvlevels[ipv]) / (pv[it, idx - 1, ila, ilo] - pv[it, idx, ila, ilo])
                            gph2[it, ipv, ila, ilo] = gph[it, idx - 1, ila, ilo] * (1 - w1) + gph[it, idx, ila, ilo] * w1
                            pre2[it, ipv, ila, ilo] = pre[it, idx - 1, ila, ilo] * (1 - w1) + pre[it, idx, ila, ilo] * w1
                            the2[it, ipv, ila, ilo] = the[it, idx - 1, ila, ilo] * (1 - w1) + the[it, idx, ila, ilo] * w1
                            break
""")["comp"]

ifn, ofn = sys.argv[1:3]

pv_levels = [2., 4., 6.]

with nc.Dataset(ifn) as inf:
    pv = abs(inf.variables["PV"][:])
    lon = inf.variables["lon"][:]
    lat = inf.variables["lat"][:]
    shape = pv.shape
    gph = inf.variables["GPH"][:]
    pre = inf.variables["PRESS"][:]
    the = inf.variables["THETA"][:]
    newshape = list(shape)
    newshape[1] = len(pv_levels)
    gph2 = np.full(newshape, np.nan, dtype=np.float32)
    pre2 = np.full(newshape, np.nan, dtype=np.float32)
    the2 = np.full(newshape, np.nan, dtype=np.float32)

    comp(np.asarray(pv_levels), pv, gph, pre, the, gph2, pre2, the2)

    with nc.Dataset(ofn, "w", format="NETCDF4_CLASSIC") as onf:
        for dim in ("time", "lat", "lon"):
            ncdim = onf.createDimension(dim, len(inf.dimensions[dim]))
            ncvar = onf.createVariable(dim, np.float64, (dim,))
            ncvar.standard_name = inf.variables[dim].standard_name
            ncvar.units = inf.variables[dim].units
            ncvar[:] = inf.variables[dim][:]
        ncdim = onf.createDimension("pv", len(pv_levels))
        ncvar = onf.createVariable('pv', np.float64, ('pv',))
        ncvar.standard_name = "atmosphere_ertel_potential_vorticity_coordinate"
        ncvar.units = "PVU"
        ncvar[:] = pv_levels

        for var, data in (("GPH", gph2), ("PRESS", pre2), ("THETA", the2)):
            ncvar = onf.createVariable(
                var, "f4", ('time', 'pv', 'lat', 'lon',),
                fill_value=np.nan, zlib=1, shuffle=1, fletcher32=1)
            ncvar.standard_name = inf.variables[var].standard_name
            ncvar.units = inf.variables[var].units
            ncvar[:] = data
    # lo, la = [inf.variables[x][:] for x in ("lon", "lat")]
    # plt.subplot(3, 1, 1)
    # plt.pcolormesh(lo, la, gph2[0,0] / 9.81 / 1000, vmax=22)
    # plt.colorbar()
    # plt.subplot(3, 1, 2)
    # plt.pcolormesh(lo, lat, gph2[0,1] / 9.81 / 1000, vmax=22)
    # plt.colorbar()
    # plt.subplot(3, 1, 3)
    # plt.pcolormesh(lo, la, gph2[0,2] / 9.81 / 1000, vmax=22)
    # plt.colorbar()
    # plt.show()
