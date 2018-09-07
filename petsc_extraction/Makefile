

SS=mpicxx -g -O3 -std=c++11 

PETSC_DIR ?= /home/boneill/Documents/ML_NEAMS/petsc/petsc_local/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables


SS_FLAGS = $(PETSC_CCPPFLAGS)  
SS_LFLAGS =$(PETSC_LIB)

all: full RS1 RS2 RSKANIKA matrix_free_full matrix_free_RS1 matrix_free_RS2

.PHONY: full, RS1, RS2
full: petsc_feature_extract.C
	$(SS) -DFULL=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
RS1: petsc_feature_extract.C
	$(SS) -DRS1=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
RS2: petsc_feature_extract.C
	$(SS) -DRS2=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
RSKANIKA: petsc_feature_extract.C
	$(SS) -DRSKANIKA=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

.PHONY: snes_ex5
snes_ex5:
	$(SS) -DFULL $(LOCINCS) $(SS_FLAGS) -o $@.a snes_ex5_extraction.C $(SS_LFLAGS)

.PHONY: matrix_free
matrix_free_RS1: petsc_feature_extract_mf.C 
	$(SS) -DRS1=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
matrix_free_RS2: petsc_feature_extract_mf.C 
	$(SS) -DRS2=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
matrix_free_full: petsc_feature_extract_mf.C 
	$(SS) -DFULL=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
clean: 
	rm -f full.a RS1.a RS2.a matrix_free_RS1.a matrix_free_RS2.a matrix_free_full.a	

