
# Set the compiler and the flags. Note that c++11 support is required. 
SS=mpicxx 
SS_FLAGS = -g -fopenmp -std=c++11 

# The directory this makefile is is (automatic)!
SS_SRC_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
SS_LIB=sslite

# Directory to install the library and dependencies 
SS_INSTALL_DIR = /usr/local/solverselecter-lite/

# Set up PETSC  
PETSC_DIR?=/opt/moose/petsc/mpich_petsc-3.7.5/gcc-opt-superlu
include $(PETSC_DIR)/lib/petsc/conf/variables
SS_DEP_FLAGS += -DWITH_PETSCUI=1 $(PETSC_CCPPFLAGS)  
SS_DEP_LFLAGS +=$(PETSC_LIB)

#### This is the commands to call to compile the library
SS_INCLUDE = -I $(SS_INSTALL_DIR)/include $(SS_DEP_FLAGS) 
SS_LIBS = -L$(SS_INSTALL_DIR)/lib -l$(SS_LIB) $(SS_DEP_LFLAGS)


