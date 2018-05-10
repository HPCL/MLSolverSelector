

SS=mpicxx -g -O3 -fPIC -std=c++11 

PETSC_DIR ?= /home/boneill/Documents/ML_NEAMS/petsc/petsc_local/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables


SS_FLAGS = $(PETSC_CCPPFLAGS)  
SS_LFLAGS =$(PETSC_LIB)

all: matrix_market 

.PHONY: matrix_market 
matrix_market: MM2Petsc.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

.PHONY: clean
clean: 
	rm -f matrix_market.a 
