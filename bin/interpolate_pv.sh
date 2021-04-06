#!/bin/bash

# Example pressure
# Inverse e.g. press(ml) and temp(ml) to temp(press)
# Afterwards interpolate temp(press) to wanted press values
ncks -v TEMP,PRESS /home/mayb/Desktop/MSS/github_data-retrieval/mss-data-retrieval/test_env/original.ml.nc /home/mayb/Desktop/MSS/github_data-retrieval/mss-data-retrieval/test_env/tmp.ml.nc
ncap2 -s "plev={850, 500};TEMP2=array(0,0,/$time,$plev,$lat,$lon/);"
