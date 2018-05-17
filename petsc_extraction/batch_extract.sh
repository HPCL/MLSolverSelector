#!/bin/bash

# Directory containing the matrices we want to test
MATRIX_DIR=./data

# The file extension for the matricies in that dir 
MATRIX_EX=pbin

# the output file name 
OUTPUTFILE=matrix_data.data

# How do we want to run it  
PROCS=1

for i in $MATRIX_DIR/*.$MATRIX_EX; do
  echo RUNNING Matrix $i 
  #./matrix_binary.a $i $OUTPUTFILE
  mpirun -n $PROCS ./matrix_binary.a $i $OUTPUTFILE
done



