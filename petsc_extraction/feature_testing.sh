
#!/bin/bash

###Ussage ./batch_extract matrix-dir outputfile matrix-ex edge-points interior-points rand solve

### Download the matrices from the FSM collection and convert them to petsc binary. This downloads all the
### matrices in the exact_mat_features.csv file if they don't already exist in dir

setup () {
export FSET=$1
export PROCS=1
export EXACTFEATURES="./exact_mat_features.csv" 
export WORKING_DIR="./WORKING/"
export MATDIR=$WORKING_DIR/matrices/
export RESDIR=$WORKING_DIR/$FSET/results/

mkdir -p $WORKING_DIR
mkdir -p $MATDIR
mkdir -p $RESDIR
make $FSET
}

# This function downloads the matrices listed in the first column of EXACTFEATURES file
# and downloads them from the FSMC and converts them to binary, installing in MATDIR, 
# then deleting all the .tar.gz files. 
download_matrices() {
  python ../petsc_mm_converter/petsc_mm_convert.py $EXACTFEATURES $MATDIR
}

# Run a single matrix with <matrix> <output> <edges> <interior> <solve> 
single_run () {
  mpirun -n $PROCS ./$FSET.a $1 $2 $3 $4 $5 
} 

# Run all matrices with extension .pbin in MatDIR with <output> <edges> <interior> <solve>
loop_de_loop() { 
for i in $MATDIR/*.pbin; do
  echo RUNNING Matrix $i 
  single_run $i $1 $2 $3 $4 
done
}

# Solves all the matrices in the matrix dir, uses a cheap 10,10 sample set to save time. 
solve_all_matrices() {
  loop_de_loop $WORKING_DIR/solver_timing 10 10 1 
}

# Debug features for 1 matrix. Basically, this runs a test for a single matrix, calculates 
# all the features for the matrix and compares to the exact solution. The output is a table 
# showing exact, calcuated and error. 
debug_feature_extraction () {
  if [ $FSET != full ] ; then
      setup full 
  fi
  mkdir -p $WORKING_DIR/DEBUG/
  rm -f $WORKING/DEBUG/features_temp
  single_run $1 $WORKING_DIR/DEBUG/features_temp -1 -1 0 

  python ./process_data.py Verify $EXACTFEATURES $WORKING_DIR/DEBUG/features_temp 0.1 

}


# Extracts featues from all matrices using the full rows (expensive) sample
run_all_rows() {
  loop_de_loop $WORKING_DIR/features_all_interior -1 -1 0 
}

# Extract features from all matrices using $1 $2 edge,interior sample set
run_all_rows_fixed() {
  loop_de_loop $RESDIR/features_$1\_$2 $1 $2 0
}

# Run the process_data script to check the accuracy of all features against Exact features
check_exact_calcs_match_exact_repo() {
  python process_data.py Verify $EXACTFEATURES $WORKING_DIR/features_all_interior 0.1 
}

# Run a few differnt sample sets. 
test_feature_accuracy() {
    run_all_rows_fixed 10 10    # Sample set S(10,10)
    run_all_rows_fixed 10 20    # Sample set S(10,20)
    run_all_rows_fixed 10 100   # Sample set S(10,100)
    run_all_rows_fixed 10 500   # Sample set S(10,500)
}


##### Plotting functions. 


# Plot the accuracy of the features against the correct exact values in EXACTFEATURES for a list 
# of features
#
# Needs tests: 
#   run_all_rows_fixed 10 10
#   run_all_rows_fixed 10 20
#   run_all_rows_fixed 10 100
#   run_all_rows_fixed 10 500
#   
plot_feature_accuracy() {
  echo $WORKING_DIR
  python ./process_data.py Error $WORKING_DIR/features_all_interior $RESDIR/features\_ 3 10 10 10 100 10 500 nnz AvgNonzerosPerRow AbsoluteNonZeroSum Trace AbsoluteTrace DiagonalMean 
}


# Plot the extraction times and ratio for a list of features
# Needs tests :
#   solve_all_matrices 
#   run_all_rows_fixed 10 10
#   run_all_rows_fixed 10 20
#   run_all_rows_fixed 10 100
#   run_all_rows_fixed 10 500
plot_feature_extraction_time () { 
  python process_data.py Extraction $WORKING_DIR/features_all_interior $WORKING_DIR/solver_timing_timing $RESDIR/features_ 3 10 10 10 100 10 500 
}

# Compare the extraction times for the different features sets :
#
# Needs tests:  
#   run_all_rows_fixed 10 10
#   run_all_rows_fixed 10 20
#   run_all_rows_fixed 10 100
#   run_all_rows_fixed 10 500
#
#   for each of the three feature sets. ( call setup <FSET> to switch feature set for testing ) 
compare_extraction_fset () {
  python process_data.py Compare $WORKING_DIR/features_all_interior $WORKING_DIR 3 full/results/features_ RS1/results/features_ RS2/results/features_ 3 10 10 10 100 10 500 
}


do_all_tests() {

  setup full 
  test_feature_accuracy      ## Get the sample solutions

  setup RS1                  ## setup env for RS1
  test_feature_accuracy      ## Get the sample solutions

  setup RS2
  test_feature_accuracy      ## same for RS2

  setup full                 ## switch back to the full ( might not be needed ) 
}

generate_all_figures() {
  plot_feature_accuracy
  plot_feature_extraction_time
  compare_extraction_fset 
}
  
