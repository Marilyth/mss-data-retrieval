#!/bin/bash
source /home/mayb/miniconda3/bin/activate 3dgrid
set -x

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=${OMP_NUM_THREADS}
export NUMEXPR_NUM_THREADS=${OMP_NUM_THREADS}
export OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS}
export VECLIB_MAXIMUM_THREADS=${OMP_NUM_THREADS}

export WORK=/home/mayb/Desktop/MSS/github_data-retrieval/mss-data-retrieval
export MSS=/home/mayb/mss/testdata
export GRIB2NCDF=$WORK/bin/grib2ncdf.sh
export ADDPV=$WORK/bin/add_pv.py
export INTP=$WORK/bin/interpolate_missing_variables.py

cd $WORK

export DATE=$1
export TIME=$2

echo Waiting for $3

if [ x${DATE} == "xtoday" ]; then
    export DATE=`date +%Y-%m-%d`
elif [ x${DATE} == "xyesterday" ]; then
    export DATE=`date +%Y-%m-%d -d yesterday`
fi

if [ x$3 != "x" ]; then
    sleep $3
fi

echo Initiating download

export BASE=${DATE}T${TIME}.an
export GRIB=grib/${BASE}.grib
export mlfile=mss/${BASE}.ml.nc
export plfile=mss/${BASE}.pl.nc
export alfile=mss/${BASE}.al.nc
export tlfile=mss/${BASE}.tl.nc
export pvfile=mss/${BASE}.pv.nc
export btfile=mss/${BASE}.bt.nc
export sfcfile=mss/${BASE}.sfc.nc
export s3dfile=mss/${BASE}.s3d.nc
export tmpfile=mss/.${BASE}.tmp.nc


if [ ! -f grib/${BASE}.ml.grib ]; then
    python bin/download_an_all.py $DATE $TIME ml &
fi
if [ ! -f grib/${BASE}.sfc.grib ]; then
    python bin/download_an_all.py $DATE $TIME sfc &
fi
if [ ! -f grib/${BASE}.pv.grib ]; then
    python bin/download_an_all.py $DATE $TIME pv &
fi
if [ ! -f grib/${BASE}.tl.grib ]; then
    python bin/download_an_all.py $DATE $TIME tl &
fi
wait

if [ ! -f grib/${BASE}.ml.grib ]; then
   echo	FATAL `date` Model level file is missing
   exit
fi
if [ ! -f grib/${BASE}.sfc.grib ]; then
   echo	FATAL `date` Surface file is missing
   exit
fi

cat grib/${BASE}.ml.grib grib/${BASE}.sfc.grib > ${GRIB}

for fn in ./$mlfile ./$pvfile ./$tmpfile ./$plfile ./$tlfile ./$alfile ./$sfcfile; do
    if [ -f $fn ]; then
        rm $fn
    fi
done

export gph_levels=0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,3100,3150,3200,3250,3300,3350,3400,3450,3500,3550,3600,3650,3700,3750,3800,3850,3900,3950,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,5900,6000
init_date=$(date +%Y-%m-%dT%H:%M:%S)
if [[ "$init_date" > "$DATE" ]] 
then 
    init_date="${DATE}T${TIME}"
fi
export time_units="hours since ${init_date}"

# convert grib to netCDF
#python $GRIB2NCDF ${GRIB} --output $mlfile --time-units "$time_units" --pressure-units Pa --gph-units "m^2s^-2" --level-type level --format NETCDF3_64BIT_OFFSET
$GRIB2NCDF input=${GRIB} output=$mlfile time_units="$time_units" pressure_units=Pa gph_units="m^2s^-2" level_type=level
cdo -f nc4 copy grib/${BASE}.tl.grib $tlfile
ncrename -d .lev,level -v .lev,level $tlfile
ncrename -v .pres,PRESS $tlfile
ncatted -a units,time,o,c,$time_units $tlfile
ncrename -v .d,DIVERGENCE $tlfile
ncrename -v .q,Q $tlfile
ncrename -v .v,V $tlfile
ncrename -v .o3,O3 $tlfile
ncrename -v .w,W $tlfile
ncrename -v .u,U $tlfile
ncrename -v .pv,PV $tlfile
cdo -f nc4 copy grib/${BASE}.pv.grib $pvfile
ncrename -d .lev,level -v .lev,level $pvfile
ncrename -v .pres,PRESS $pvfile
ncatted -a units,time,o,c,$time_units $pvfile
ncrename -v .d,DIVERGENCE $pvfile
ncrename -v .q,Q $pvfile
ncrename -v .v,V $pvfile
ncrename -v .o3,O3 $pvfile
ncrename -v .w,W $pvfile
ncrename -v .u,U $pvfile
ncrename -v .pt,THETA $pvfile

# Add ancillay information
python $ADDPV MSSL $mlfile --pv --theta --tropopause --n2 #--eqlat nan values cause issues for now due to no 180° coverage

# separate sfc from ml variables
ncks -7 -L 4 -C -O -x -vlevel,N2,clwc,U,Q,TEMP,PRESS,GPH,cc,W,V,ciwc,THETA,PV,MOD_PV,O3,DIVERGENCE $mlfile $sfcfile
ncatted -O -a standard_name,msl,o,c,air_pressure_at_sea_level $sfcfile
ncks -6 -C -O -vtime,level,lon,lat,N2,clwc,U,Q,TEMP,PRESS,GPH,cc,W,V,ciwc,THETA,PV,MOD_PV,O3,DIVERGENCE,hyai,hyam,hybi,hybm,sp,lnsp $mlfile $tmpfile
mv $tmpfile $mlfile

# interpolate to different grids
# pressure levels
cdo ml2pl,85000,50000,40000,30000,20000,15000,12000,10000,8000,6500,5000,4000,3000,2000,1000,500,100 $mlfile $plfile
ncatted -O -a standard_name,plev,o,c,atmosphere_pressure_coordinate $plfile
ncks -C -O -x -v lev,sp,lnsp $plfile $plfile
mv $plfile ${MSS}/EUR_LL015.an.pl.nc

# theta levels
python $INTP $mlfile $tlfile GPH THETA
python $INTP $mlfile $tlfile N2 THETA
python $INTP $mlfile $tlfile TEMP THETA
ncatted -O -a standard_name,level,o,c,atmosphere_potential_temperature_coordinate $tlfile
mv $tlfile ${MSS}/EUR_LL015.an.tl.nc

# potential vorticity levels
python $INTP $mlfile $pvfile GPH PV
python $INTP $mlfile $pvfile N2 PV
python $INTP $mlfile $pvfile TEMP PV
ncatted -O -a standard_name,level,o,c,atmosphere_ertel_potential_vorticity_coordinate $pvfile
mv $pvfile ${MSS}/EUR_LL015.an.pv.nc

# altitude levels
ncks -6 -C -O -vtime,level,lon,lat,N2,U,TEMP,PRESS,GPH,W,V,THETA,PV,hyai,hyam,hybi,hybm,lnsp $mlfile $tmpfile
cdo ml2hl,$gph_levels $tmpfile $alfile
ncatted -O -a standard_name,height,o,c,atmosphere_altitude_coordinate $alfile
ncks -C -O -x -v lev,sp,lnsp $alfile $alfile
mv $alfile ${MSS}/EUR_LL015.an.al.nc
rm $tmpfile

ncks -6 -O -d level,0,0 -d level,16,28,4 -d level,32,124,2 $mlfile $tmpfile
rm $mlfile
nccopy -7 -s -d7 $tmpfile $mlfile
rm $tmpfile
ncatted -O -a standard_name,level,o,c,atmosphere_hybrid_sigma_pressure_coordinate $mlfile
ncks -C -O -x -v lev,sp,lnsp,nhyi,nhym,hyai,hyam,hybi,hybm $mlfile $mlfile
mv $mlfile ${MSS}/EUR_LL015.an.ml.nc
mv $sfcfile ${MSS}/EUR_LL015.an.sfc.nc
