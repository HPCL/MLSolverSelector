
include ../../Makefile.inc
PETSCSCAN= petsc_scan.C 
PETSCSCANO= petsc_scan.a 

MM2PETSC= MM2Petsc.C
MM2PETSCO= MM2Petsc.a


all: $(PETSCSCANO) $(MM2PETSCO)

%.a : %.C 
	$(CC) $(SS_FLAGS) -o $@ $^ $(SS_LFLAGS)

.PHONY : clean
clean:
	rm -f $(PETSCSCAN0) $(MM2PETSCO)

