#!/bin/bash
source /home/mayb/miniconda3/bin/activate 3dgrid
set -x

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=${OMP_NUM_THREADS}
export NUMEXPR_NUM_THREADS=${OMP_NUM_THREADS}
export OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS}
export VECLIB_MAXIMUM_THREADS=${OMP_NUM_THREADS}

export WORK=/home/mayb/Desktop/MSS/github_data-retrieval/mss-data-retrieval
export GRIB2NCDF=$WORK/bin/grib2ncdf.sh
export INTPMOD=$WORK/bin/scripts/interpolate_model.py
export ADDPV=$WORK/bin/add_pv2.py
export CONVW=$WORK/bin/convert_omega_to_w.py
export INTPPV=$WORK/bin/interpolate_pv.py

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

# Add ancillay information
python $ADDPV MSSL $mlfile --pv --theta --tropopause --n2 #--eqlat nan values cause issues for now

# separate sfc from ml variables
ncks -7 -L 4 -C -O -x -vlevel,N2,clwc,U,Q,TEMP,PRESS,GPH,cc,W,V,ciwc,THETA,PV,MOD_PV,O3,DIVERGENCE $mlfile $sfcfile
ncatted -O -a standard_name,msl,o,c,air_pressure_at_sea_level $sfcfile
ncks -6 -C -O -vtime,level,lon,lat,N2,clwc,U,Q,TEMP,PRESS,GPH,cc,W,V,ciwc,THETA,PV,MOD_PV,O3,DIVERGENCE,hyai,hyam,hybi,hybm,sp,lnsp $mlfile $tmpfile
mv $tmpfile $mlfile

# interpolate to different grids
#python $INTPMOD -v 1 $mlfile $plfile --level-type pressure --vert-unit hPa --levels 850,500,400,300,200,150,120,100,80,65,50,40,30,20,10,5,1
cdo ml2pl,85000,50000,40000,30000,20000,15000,12000,10000,8000,6500,5000,4000,3000,2000,1000,500,100 $mlfile $plfile
#python $INTPMOD -v 1 $mlfile $tlfile --level-type theta --vert-unit K --levels 340,360,370,380,390,400,410,420
python $INTPPV $mlfile $pvfile
ncks -6 -C -O -vtime,level,lon,lat,N2,U,TEMP,PRESS,GPH,W,V,THETA,PV,hyai,hyam,hybi,hybm,lnsp $mlfile $tmpfile
#python $INTPMOD -v 1 $tmpfile $alfile --level-type gph --vert-unit km --levels $gph_levels
cdo ml2hl,$gph_levels $tmpfile $alfile
rm $tmpfile

ncks -6 -O -d level,0,0 -d level,16,28,4 -d level,32,124,2 $mlfile $tmpfile
rm $mlfile
nccopy -7 -s -d7 $tmpfile $mlfile
rm $tmpfile
