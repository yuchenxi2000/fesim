#/bin/sh
module load intel/2018.1
cd "$(dirname "$0")"
icpc -std=c++11 -fopenmp -O2 -lmkl_rt -o main main.cpp

