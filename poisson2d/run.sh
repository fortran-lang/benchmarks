#!/bin/bash

# This may not apply to your distribution
# change or comment out at will, but remark
# the stack-size.
ulimit -s unlimited

OPTS="-O3 -march=native"
FC=gfortran
CC=gcc

$CC $OPTS -c c_common.c
for c in 1 2 3 4_1d 4_2d
do
    echo "Compiling and running: c_${c}"
    $CC $OPTS -o c_${c} c_${c}.c c_common.o

    # Now run for
    {
	for N in 100 200 300 ; do
	    ./c_${c} $N 10000
	done
	for N in 400 500 600; do
	    ./c_${c} $N 4000
	done
    } > c_${c}.data
done

$FC $OPTS -c fortran_common.f90
for f in 2 3 4
do
    echo "Compiling and running: fortran_${f}"
    $FC $OPTS -o fortran_${f} fortran_${f}.f90 fortran_common.o

    # Now run for
    {
	for N in 100 200 300 ; do
	    ./fortran_${f} $N 10000
	done
	for N in 400 500 600; do
	    ./fortran_${f} $N 4000
	done
    } > fortran_${f}.data
done
