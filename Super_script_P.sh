#!/bin/bash
rm exe.x

gfortran -O2 -o exe.x covid.f95 Subroutines_Covid.f95 dranxor2.f Random2.f

mkdir out
mkdir err
city=("Results")
#alphas=("0")
mkdir ${city}
mkdir ${city}/"sum"
mkdir ${city}/"T"
realizations=1 #Number of realizations to execute
for  ((jjcontrol=1;jjcontrol<=$realizations;jjcontrol++)); do
		time echo $jjcontrol | ./exe.x #For Local
done
