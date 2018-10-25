
PETSC_DIR ?=/opt/moose/petsc/mpich_petsc-3.7.5/gcc-opt-superlu

include $(PETSC_DIR)/lib/petsc/conf/variables
SS_DEP_FLAGS += $(PETSC_CCPPFLAGS)  
SS_DEP_LFLAGS +=$(PETSC_LIB)
PETSCSCAN= petsc_scan.C 
PETSCSCANO= petsc_scan.a 

MM2PETSC= MM2Petsc.C
MM2PETSCO= MM2Petsc.a


all: $(PETSCSCANO) $(MM2PETSCO)

%.a : %.C 
	$(CC) $(PETSC_CCPPFLAGS) -o $@ $^ $(PETSC_LIB)

.PHONY : clean
clean:
	rm -f $(PETSCSCAN0) $(MM2PETSCO)

