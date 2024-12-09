#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Usage: ./run.sh <q> <p>"
    exit 1
fi

q=$1
p=$2
mpicc -o start start.c -lm
# mpiexec --use-hwthread-cpus -n $((2**p)) ./start $q $p
mpiexec --oversubscribe -n $((2**p)) ./start $q $p