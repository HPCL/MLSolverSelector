

SS=mpicxx -g -O3 -fPIC -std=c++11 

PETSC_DIR ?= /home/boneill/Documents/ML_NEAMS/petsc/petsc_local/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables


SS_FLAGS = $(PETSC_CCPPFLAGS)  
SS_LFLAGS =$(PETSC_LIB)

all: matrix_binary matrix_free 

.PHONY: matrix_binary 
matrix_binary: petsc_feature_extract.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

.PHONY: matrix_free
matrix_free: petsc_feature_extract_mf.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

clean: 
	rm -f matrix_free.a matrix_binary.a
	

