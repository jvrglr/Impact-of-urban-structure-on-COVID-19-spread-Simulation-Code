#!/bin/bash
rm exe.x
rm exeNoT.x

gfortran -O2 -o exe.x covid.f95 Subroutines_Covid.f95 dranxor2.f Random2.f
gfortran -O2 -o exeNoT.x covid_no_traj.f95 Subroutines_Covid.f95 dranxor2.f Random2.f
mkdir out
mkdir err
city=("Results")
#alphas=("0")
mkdir ${city}
mkdir ${city}/"sum"
mkdir ${city}/"T"

for jjcontrol in {0..0}; do
		time echo $jjcontrol | ./exe.x #For Local
done
