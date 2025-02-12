#!/bin/bash

WORKING_DIR=$(pwd)

cd ../../../../mesh

rm domain2D.msh
./test.py naca.dat 0.4 $1
gmsh NACA_2408.geo -2 -o domain2D.msh


cd ${WORKING_DIR}
module load gcc-glibc dealii

mpirun -n 2 main
mv forces_vs_time.csv output_${1}.csv


