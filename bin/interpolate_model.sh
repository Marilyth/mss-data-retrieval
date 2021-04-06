# Dynamically assign named command line arguments
for ARGUMENT in "$@"
do
    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=) 
    declare "$KEY=$VALUE"     
done

# Interpolate linear the desired levels for each lat/lon combination and save in variable
# cdo remapbic/remapbil
ncap2 -s 'defdim("level2");level2[level2]={$levels};level2@units=$vert_unit;' $in $out

