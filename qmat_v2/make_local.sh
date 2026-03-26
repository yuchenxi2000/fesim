#!/bin/sh
cd "$(dirname "$0")"
mpic++ -O2 -o qmat -std=c++11 -DACCELERATE_NEW_LAPACK -DACCELERATE_LAPACK_ILP64 -framework Accelerate main.cpp 
