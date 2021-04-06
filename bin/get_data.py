import sys
import os
from datetime import datetime, date, timedelta
from dateutil import parser
import subprocess

os.environ["OMP_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = os.environ["OMP_NUM_THREADS"]
os.environ["NUMEXPR_NUM_THREADS"] = os.environ["OMP_NUM_THREADS"]
os.environ["OPENBLAS_NUM_THREADS"] = os.environ["OMP_NUM_THREADS"]
os.environ["VECLIB_MAXIMUM_THREADS"] = os.environ["OMP_NUM_THREADS"]

WORK = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))
GRIB2NCDF = WORK + "/bin/scripts/grib2ncdf.py"
INTPMOD = WORK + "/bin/scripts/interpolate_model.py"
ADDPV = WORK + "/bin/add_pv.py"
CONVW = WORK + "/bin/convert_omega_to_w.py"
INTPPV = WORK + "/bin/interpolate_pv.py"

if sys.argv[1] == "today":
    DATE = datetime.combine(date.today(), datetime.min.time())
    DATE += timedelta(hours=int(sys.argv[2]))
elif sys.argv[1] == "yesterday":
    DATE = datetime.combine(date.today() - timedelta(days=1), datetime.min.time())
    DATE += timedelta(hours=int(sys.argv[2]))
else:
    DATE = parser.parse(" ".join(sys.argv[1:]))

print("Initiating download")

BASE = f"{DATE.strftime('%Y-%m-%dT%H%M%S')}.an"
GRIB = WORK + f"/grib/{BASE}.grib"
mlfile = WORK + f"/mss/{BASE}.ml.nc"
plfile = WORK + f"/mss/{BASE}.pl.nc"
alfile = WORK + f"/mss/{BASE}.al.nc"
tlfile = WORK + f"/mss/{BASE}.tl.nc"
pvfile = WORK + f"/mss/{BASE}.pv.nc"
btfile = WORK + f"/mss/{BASE}.bt.nc"
sfcfile = WORK + f"/mss/{BASE}.sfc.nc"
s3dfile = WORK + f"/mss/{BASE}.s3d.nc"
tmpfile = WORK + f"/mss/.{BASE}.tmp.nc"

command = []
if not os.path.isfile(WORK + f"/grib/${BASE}.ml.grib"):
    command.append("python " + WORK + f"/bin/download_an_ml.py {DATE.strftime('%Y-%m-%d %H:%M:%S')}")

if not os.path.isfile(WORK + f"/grib/${BASE}.sfc.grib"):
    command.append("python " + WORK + f"/bin/download_an_sfc.py {DATE.strftime('%Y-%m-%d %H:%M:%S')}")
subprocess.run(" & ".join(command), shell=True, capture_output=True)

if not os.path.isfile(WORK + f"/grib/${BASE}.ml.grib"):
   print("FATAL `date` Model level file is missing")
   exit()

if not os.path.isfile(WORK + f"/grib/${BASE}.sfc.grib"):
   print("FATAL `date` Surface file is missing")
   exit()
"""
cat grib/${BASE}.ml.grib grib/${BASE}.sfc.grib > ${GRIB}

for fn in ./$mlfile ./$pvfile ./$tmpfile ./$plfile ./$tlfile ./$alfile ./$btfile ./$sfcfile ./$s3dfile; do
    if [ -f $fn ]; then
        rm $fn
    fi
done

export gph_levels=0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9,9.25,9.5,9.75,10,10.25,10.5,10.75,11,11.25,11.5,11.75,12,12.25,12.5,12.75,13,13.25,13.5,13.75,14,14.25,14.5,14.75,15,15.25,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,19.75,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60
init_date=$(date +%Y-%m-%dT%H:%M:%S)
if [[ "$init_date" > "$DATE" ]] 
then 
    init_date="${DATE}T${TIME}" 
fi
export time_units="hours since ${init_date}"

# convert grib to netCDF
python $GRIB2NCDF ${GRIB} --output $mlfile --time-units "$time_units" --pressure-units Pa --gph-units "m^2s^-2" --level-type level --format NETCDF3_64BIT_OFFSET
rm ${GRIB}

# Add ancillay information
python $ADDPV MSSL $mlfile --pv --theta --tropopause --n2 --eqlat

# separate sfc from ml variables
ncks -7 -L 4 -C -O -x -vlevel,N2,CLWC,U,Q,TEMP,PRESS,GPH,CC,W,V,CIWC,THETA,PV,MOD_PV,O3,DIVERGENCE,EQLAT $mlfile $sfcfile
ncatted -O -a standard_name,MSL,o,c,air_pressure_at_sea_level $sfcfile
ncks -6 -C -O -vtime,level,lon,lat,N2,CLWC,U,Q,TEMP,PRESS,GPH,CC,W,V,CIWC,THETA,PV,MOD_PV,O3,DIVERGENCE,EQLAT $mlfile $tmpfile
mv $tmpfile $mlfile

# interpolate to different grids
python $INTPMOD -v 1 $mlfile $plfile --level-type pressure --vert-unit hPa --levels 850,500,400,300,200,150,120,100,80,65,50,40,30,20,10,5,1
python $INTPMOD -v 1 $mlfile $tlfile --level-type theta --vert-unit K --levels 340,360,370,380,390,400,410,420
python $INTPPV $mlfile $pvfile
ncks -6 -C -O -vtime,level,lon,lat,N2,U,TEMP,PRESS,GPH,W,V,THETA,PV $mlfile $tmpfile
python $INTPMOD -v 1 $tmpfile $alfile --level-type gph --vert-unit km --levels $gph_levels
rm $tmpfile

ncks -6 -O -d level,0,0 -d level,16,28,4 -d level,32,124,2 $mlfile $tmpfile
rm $mlfile
nccopy -7 -s -d7 $tmpfile $mlfile
rm $tmpfile"""
