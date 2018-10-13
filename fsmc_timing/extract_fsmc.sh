
#!/bin/bash

#Run a single matrix with <matrix> <output> <edges> <interior> 
single_run () {
 ./$FSET.a $1 $2 $3 $4  
} 

# Run all matrices in fsmc directory <output> <edges> <interior> 
loop_de_loop() { 
for i in fsmc/*.petsc; do
  echo RUNNING Matrix $i 
  single_run $i $1 $2 $3 
done
}

# Run a set of tests testing sample size with input FSET 
run_tests() {

  export FSET=$1
  mkdir -p $2/${FSET}  
  
  #Fixed Value Tests 
  loop_de_loop $2/${FSET}/10_10 10 10 
  loop_de_loop $2/${FSET}/10_20 10 20
  loop_de_loop $2/${FSET}/10_50 10 50
  loop_de_loop $2/${FSET}/10_75 10 75
  loop_de_loop $2/${FSET}/10_100 10 100
  loop_de_loop $2/${FSET}/10_250 10 250
  loop_de_loop $2/${FSET}/10_500 10 500

  # Percentage tests 
  loop_de_loop $2/${FSET}/p_10_10 10 .10 
  loop_de_loop $2/${FSET}/p_10_20 10 .20
  loop_de_loop $2/${FSET}/p_10_50 10 .50
  loop_de_loop $2/${FSET}/_p_10_75 10 .75
}


do_all_tests() {
   run_tests RSKANIKA $1
   run_tests RS1 $1
   run_tests RS2 $1
   run_tests FULL $1
}

do_all_tests ./testResults  



