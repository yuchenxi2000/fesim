#!/bin/sh
cd "$(dirname "$0")"
module add intel/2018.0
module add openmpi/3.0.0-intel2018.0 
mpic++ -O2 -o qmat -std=c++11 -lmkl_rt main.cpp 

