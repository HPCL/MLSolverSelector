# 1. Set PETSC_DIR (export PETSC_DIR= in env should be enough. ) 
# 2. run make install.deps
# 3. run make install 


SRC_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
INSTALL_DIR=$(SRC_DIR)/build/solverselecter
EXTRA_DIR=$(SRC_DIR)/build/temp/
EXT_BUILD=$(SRC_DIR)/build/
SSLIBNAME='solverselecter'

SS=mpicxx -g -O3 -std=c++11 
SS_FLAGS= 
SS_LFLAGS= -L $(INSTALL_DIR)/lib -l$(SSLIBNAME) 

### The names of the subdirs to search for *.C files to compile. 
MODNAMES = base features measurements database machine_learning interfaces

##### The makefile is set up as follows. Additional library
# dependencies are set with the WITH_XXX variables. The make 
# file attempts to build every *.C file in the src directories
# of the $(MODFILES) subdirectories. Implimentations can be turned
# on and off with the WITH_XXX variables. Currently, the WITH_XXX variables
# are used as PRAGMAs inside the code to decide when to add objects. 

#### There are some makefile targets  download.zoltan download.sqlite3 and download.waffles
# and install.zoltan install.sqlite3 and install.waffles that may or may not work to install
# each of these dependencies. 

# Do we have zoltan installed. This is used for matrix free (not working yet) 
WITH_ZOLTAN=1
ZOLTAN_DIR=$(EXT_BUILD)/zoltan/
ZOLTAN_DOWNLOAD=http://www.cs.sandia.gov/~kddevin/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz
ifeq ($(WITH_ZOLTAN),1)
SS_FLAGS += -DWITH_ZOLTAN=1 -I $(ZOLTAN_DIR)/include 
SS_LFLAGS += -L $(ZOLTAN_DIR)/lib -lzoltan
endif

# Do we want to build the sqlite3 database implimentation 
WITH_SQLITE3=1
SQLITE3_DIR=$(EXT_BUILD)/sqlite/
SQLITE3_DOWNLOAD=https://sqlite.org/2017/sqlite-autoconf-3210000.tar.gz
ifeq ($(WITH_SQLITE3),1)
SS_FLAGS += -DWITH_SQLITE3=1 -I $(SQLITE3_DIR)/include 
SS_LFLAGS += -L $(SQLITE3_DIR)/lib -lsqlite3
endif

# Do we want to install the waffles ML implimentation 
WITH_WAFFLES=1
WAFFLES_DIR=$(EXT_BUILD)/waffles/
WAFFLES_DOWNLOAD=https://github.com/mikegashler/waffles/archive/1.0.0.tar.gz
ifeq ($(WITH_WAFFLES),1)
SS_FLAGS += -DWITH_WAFFLES=1 -I $(WAFFLES_DIR)/include 
SS_LFLAGS += -L $(WAFFLES_DIR)/lib -lGClasses -lpthread
endif

# Do we want to build the Petsc Interface for the solver
WITH_PETSCUI=1
PETSC_DIR ?= /home/boneill/Documents/ML_NEAMS/petsc/petsc_local/petsc
ifeq ($(WITH_PETSCUI),1)
include ${PETSC_DIR}/lib/petsc/conf/variables
SS_FLAGS += -DWITH_PETSCUI=1 $(PETSC_CCPPFLAGS)  
SS_LFLAGS +=$(PETSC_LIB)
endif

# Do we want to build the DummyUI Interface for the solver 
WITH_DUMMYUI=1
ifeq ($(WITH_DUMMYUI),1)
SS_FLAGS+= -DWITH_DUMMYUI=1 
endif

### Get the source directory and names of all the subdirs
SS_FLAGS += $(addprefix -I $(SRC_DIR), $(addsuffix /include, $(MODNAMES)))


