#!/bin/bash

MATRIX_DIR=./data
MATRIX_EX=pbin
OUTPUTFILE=matrix_data.data

for i in $MATRIX_DIR/*.$MATRIX_EX; do
  echo RUNNING Matrix $i 
  ./matrix_binary.a $i $OUTPUTFILE
done


