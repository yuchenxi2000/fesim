#!/bin/sh
cd "$(dirname "$0")"
module add mpi/intel/oneapi/2021.7.1
mpicxx -O2 -o qmat -std=c++11 -lmkl_rt main.cpp 
