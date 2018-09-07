
SS=mpicxx -g -O3 -std=c++11 

PETSC_DIR ?= /home/boneill/Documents/ML_NEAMS/petsc/petsc_local/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables


SS_FLAGS = $(PETSC_CCPPFLAGS)  
SS_LFLAGS =$(PETSC_LIB)

all: full RS1 RS2 RSKANIKA 

.PHONY: full, RS1, RS2
full: extract_fsmc.C
	$(SS) -DFULL=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
RS1: extract_fsmc.C
	$(SS) -DRS1=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
RS2: extract_fsmc.C
	$(SS) -DRS2=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   
RSKANIKA: extract_fsmc.C
	$(SS) -DRSKANIKA=1 $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

clean: 
	rm -f full.a RS1.a RS2.a RSKANINA.a 	

