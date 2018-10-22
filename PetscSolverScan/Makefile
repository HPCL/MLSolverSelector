
PETSC_DIR ?=/opt/moose/petsc/mpich_petsc-3.7.5/gcc-opt-superlu

include $(PETSC_DIR)/lib/petsc/conf/variables


all: $(PETSCSCANO)

%.a : %.C 
	$(CC) $(PETSC_CCPPFLAGS) -o $@ $^ $(PETSC_LIB)

.PHONY : clean
clean:
	rm -f *.a PetscPreconOptions.db PetscSolverOptions.db petscoutput1


