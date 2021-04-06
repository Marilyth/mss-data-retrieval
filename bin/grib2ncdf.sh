#!/bin/bash

# Set all grib2ncdf.py global attributes
ECMWF_META=(
    "class" "type" "stream" "typeOfProcessedData"
    "experimentVersionNumber" "gridType" "g2grid")
declare -A VARNAME_TRANSLATE=( ["T"]="TEMP" ["D"]="DIVERGENCE" )
declare -A STANDARD_NAME_TRANSLATE=(
    ["10 metre U wind component"]="surface_eastward_wind"
    ["10 metre V wind component"]="surface_northward_wind"
    ["Boundary layer height"]="atmosphere_boundary_layer_thickness"
    ["Low cloud cover"]="low_cloud_area_fraction"
    ["Medium cloud cover"]="medium_cloud_area_fraction"
    ["High cloud cover"]="high_cloud_area_fraction"
    ["Sea-ice cover"]="sea_ice_area_fraction"
    ["Specific cloud ice water content"]="specific_cloud_ice_water_content"
    ["Specific cloud liquid water content"]="specific_cloud_liquid_water_content"
    ["Fraction of cloud cover"]="cloud_area_fraction_in_atmosphere_layer"
)
declare -A GPH_FAC=( ["m"]=1. ["km"]=1000.0 ["m^2s^-2"]=9.80665 )
declare -A PRESSURE_FAC=( ["hPa"]=100.0 ["Pa"]=1.0 )
declare -A NETCDF_PARAMS=( ["fill_value"]="np.nan" )

tmpfile="./tmp.nc"
tmpgrib="./tmp.grib"

# Dynamically assign named command line arguments
for ARGUMENT in "$@"
do
    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=) 
    declare "$KEY=$VALUE"     
done

add_pressure_gph () {
    # gheight is in meters, convert
    cdo gheight $input gph.grib
    cdo expr,"gh=gh*${GPH_FAC[$gph_units]}" gph.grib gph2.grib
    mv gph2.grib gph.grib
    
    # pressure is in Pa, convert
    cdo pressure_fl $input pressure.grib
    cdo expr,"var1=var1*${PRESSURE_FAC[$pressure_units]}" pressure.grib pressure2.grib
    mv pressure2.grib pressure.grib

    cat $input pressure.grib gph.grib > $tmpgrib
    rm gph.grib pressure.grib
}

echo "Constructing gph and pressure"
#add_pressure_gph

echo "Creating ncdf, Renaming variables, dropping unused"
cdo -f nc4 copy $tmpgrib $tmpfile
rm $tmpgrib
ncrename -d .lev_2,level -v .lev_2,level $tmpfile
ncrename -v .gh,GPH $tmpfile
ncrename -v .var1,PRESS $tmpfile
ncatted -a units,GPH,o,c,$gph_units $tmpfile
ncatted -a units,PRESS,o,c,$pressure_units $tmpfile
ncatted -a units,time,o,c,$time_units $tmpfile
ncrename -v .d,DIVERGENCE $tmpfile
ncrename -v .t,TEMP $tmpfile
ncrename -v .q,Q $tmpfile
ncrename -v .v,V $tmpfile
ncrename -v .o3,O3 $tmpfile
ncrename -v .w,W $tmpfile
ncrename -v .u,U $tmpfile
#ncwa -a nhyi $tmpfile tmp2.nc
#ncwa -a nhym tmp2.nc $tmpfile
#ncks -C -O -x -v hyai $tmpfile tmp2.nc
#ncks -C -O -x -v hybi tmp2.nc $tmpfile
#ncks -C -O -x -v hyam $tmpfile tmp2.nc
#ncks -C -O -x -v hybm tmp2.nc $tmpfile
rm tmp2.nc
mv $tmpfile $output
