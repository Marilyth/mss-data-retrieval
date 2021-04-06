import netCDF4 as nc
import sys

new_units = sys.argv[1]

for fn in sys.argv[2:]:
    with nc.Dataset(fn, "r+") as nf:
        time = nf.variables["time"]
        vals, units = time[:], time.units
        dts = nc.num2date(vals, units)
        new_vals = nc.date2num(dts, new_units)
        time[:] = new_vals
        time.units = new_units
