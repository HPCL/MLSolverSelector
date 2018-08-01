
# Set the compiler and the flags. Note that c++11 support is required. 
SS=mpicxx 
SS_FLAGS = -g -std=c++11 

# The directory this makefile is is (automatic)!
SS_SRC_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
SS_LIB=SolverSelecter

# Directory to install the library and dependencies 
SS_INSTALL_DIR = /usr/local/solverselecter/

# Do we want to use the SQLITE 3  
WITH_SQLITE3=1
SQLITE3_DIR=$(SS_INSTALL_DIR)/sqlite3
ifeq ($(WITH_SQLITE3),1)
SS_DEP_FLAGS += -DWITH_SQLITE3=1 -I $(SQLITE3_DIR)/include 
SS_DEP_LFLAGS += -L$(SQLITE3_DIR)/lib -lsqlite3
endif

# Do we want to use WAFFLES 
WITH_WAFFLES=1
WAFFLES_DIR ?=$(SS_INSTALL_DIR)/waffles
ifeq ($(WITH_WAFFLES),1)
SS_DEP_FLAGS += -DWITH_WAFFLES=1 -I $(WAFFLES_DIR)/include 
SS_DEP_LFLAGS += -L$(WAFFLES_DIR)/lib -lGClasses -lpthread
endif

# Do we want to use Petsc 
WITH_PETSC=1
PETSC_DIR=/opt/moose/petsc/mpich_petsc-3.7.5/gcc-opt-superlu
ifeq ($(WITH_PETSC),1)
include $(PETSC_DIR)/lib/petsc/conf/variables
SS_DEP_FLAGS += -DWITH_PETSCUI=1 $(PETSC_CCPPFLAGS)  
SS_DEP_LFLAGS +=$(PETSC_LIB)
endif

# Do we want to use boost file system. 
WITH_BOOST_FILESYSTEM=1
BOOST_DIR ?= /opt/moose/boost_1_63_0/
ifeq ($(WITH_BOOST_FILESYSTEM),1)
SS_DEP_FLAGS += -DWITH_BOOST_FILESYSTEM=1   
SS_DEP_LFLAGS += -lboost_filesystem -lboost_system 
endif

#### This is the commands to call to compile the library
SS_INCLUDE = -I $(SS_INSTALL_DIR)/include $(SS_DEP_FLAGS) 
SS_LIBS = -L$(SS_INSTALL_DIR)/lib -l$(SS_LIB) $(SS_DEP_LFLAGS)


