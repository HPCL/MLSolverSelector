

SS=mpicxx -g -O3 -fPIC -std=c++11 
TOP := $(dir $(lastword $(MAKEFILE_LIST)))
EXT_BUILD = $(TOP)/dependencies/
EXTRA_DIR = $(TOP)/dependencies/source/
SS_FLAGS += -I $(TOP)/base/

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
