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
tmpgrib2="./tmp2.grib"

# Dynamically assign named command line arguments
for ARGUMENT in "$@"
do
    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=) 
    declare "$KEY=$VALUE"     
done

cdo -f nc4 copy $input $tmpfile

add_pressure_gph () {
    # add surface pressure, as lnsp is somehow not working for ERA5 and gheight
    ncap2 -s 'sp=exp(lnsp);sp@units="Pa";sp@standard_name="surface_air_pressure";sp@code=134;sp@table=128' $tmpfile sp.nc
    cdo gheight sp.nc gph.nc
    # gheight is in meters, convert
    ncap2 -s "zh=zh*${GPH_FAC[$gph_units]}" gph.nc gph2.nc
    mv gph2.nc gph.nc
    
    cdo pressure_fl sp.nc pressure.nc
    # pressure is in Pa, convert
    ncap2 -s "pressure=pressure*${PRESSURE_FAC[$pressure_units]}" pressure.nc pressure2.nc
    mv pressure2.nc pressure.nc

    cdo merge sp.nc gph.nc pressure.nc merged.nc
    mv merged.nc $tmpfile
    rm gph.nc pressure.nc sp.nc
}

echo "Constructing gph and pressure"
add_pressure_gph

echo "Renaming variables, dropping unused"
ncatted -a units,zh,o,c,$gph_units $tmpfile
ncatted -a units,pressure,o,c,$pressure_units $tmpfile
ncatted -a units,time,o,c,"${time_units}" $tmpfile
mv $tmpfile $output
