include ../../Makefile.inc

LOCINCS = -I ../ 

all: petscex2 bfm petscex1

.PHONY: petscf 
petscf: petsc_feature_extract.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

.PHONY: petscmf 
petscmf: petsc_feature_extract_mf.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

.PHONY: petscex2 
petscex2: petscex2.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

.PHONY: petscex1 
petscex1: petscex1.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

.PHONY: bfm  
bfm: build_from_file.C 
	$(SS) $(LOCINCS)  $(SS_FLAGS) -o $@.a $< $(SS_LFLAGS)   

clean: 
	rm -f *.a 	

