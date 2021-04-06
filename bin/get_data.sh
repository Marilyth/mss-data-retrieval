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
export ADDPV=$WORK/bin/add_pv.py
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
    python bin/download_an_ml.py $DATE $TIME &
fi
if [ ! -f grib/${BASE}.sfc.grib ]; then
    python bin/download_an_sfc.py $DATE $TIME &
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
#python $GRIB2NCDF ${GRIB} --output $mlfile --time-units "$time_units" --pressure-units Pa --gph-units "m^2s^-2" --level-type level --format NETCDF3_64BIT_OFFSET
$GRIB2NCDF input=${GRIB} output=$mlfile time_units="$time_units" pressure_units=Pa gph_units="m^2s^-2" level_type=level

# Add ancillay information
python $ADDPV MSSL $mlfile --pv --theta --tropopause --n2 --eqlat

# separate sfc from ml variables
ncks -7 -L 4 -C -O -x -.vlevel,.N2,.CLWC,.U,.Q,.TEMP,.PRESS,.GPH,.CC,.W,.V,.CIWC,.THETA,.PV,.MOD_PV,.O3,.DIVERGENCE,.EQLAT $mlfile $sfcfile
ncatted -O -a standard_name,MSL,o,c,air_pressure_at_sea_level $sfcfile
ncks -6 -C -O -.vtime,.level,.lon,.lat,.N2,.CLWC,.U,.Q,.TEMP,.PRESS,.GPH,.CC,.W,.V,.CIWC,.THETA,.PV,.MOD_PV,.O3,.DIVERGENCE,.EQLAT $mlfile $tmpfile
mv $tmpfile $mlfile

# interpolate to different grids
python $INTPMOD -v 1 $mlfile $plfile --level-type pressure --vert-unit hPa --levels 850,500,400,300,200,150,120,100,80,65,50,40,30,20,10,5,1 # cdo ml2pl
python $INTPMOD -v 1 $mlfile $tlfile --level-type theta --vert-unit K --levels 340,360,370,380,390,400,410,420
python $INTPPV $mlfile $pvfile
ncks -6 -C -O -vtime,level,lon,lat,N2,U,TEMP,PRESS,GPH,W,V,THETA,PV $mlfile $tmpfile
python $INTPMOD -v 1 $tmpfile $alfile --level-type gph --vert-unit km --levels $gph_levels
rm $tmpfile

ncks -6 -O -d level,0,0 -d level,16,28,4 -d level,32,124,2 $mlfile $tmpfile
#rm $mlfile
nccopy -7 -s -d7 $tmpfile $mlfile
rm $tmpfile
