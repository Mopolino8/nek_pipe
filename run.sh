#!/bin/bash -l
#
# usage: run.sh casename np

echo $1 > SESSION.NAME
echo $PWD/ >> SESSION.NAME

rm $1'.sch'

mpirun -n $2 ./nek5000
