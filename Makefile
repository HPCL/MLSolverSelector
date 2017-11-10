# 1. Set PETSC_DIR (export PETSC_DIR= in env should be enough. ) 
# 2. run make install.deps
# 3. run make install 


SRC_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
INSTALL_DIR=$(SRC_DIR)/build/solverselecter
EXTRA_DIR=$(SRC_DIR)/build/temp/
EXT_BUILD=$(SRC_DIR)/build/
SSLIBNAME=solverselecter

SS=mpicxx -g -std=c++11 
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

### Get the source directory and names of all the subdirs
SS_FLAGS += $(addprefix -I $(SRC_DIR), $(addsuffix /include, $(MODNAMES)))

$(foreach name, clean lib install , $(eval $(name)_targets:=$(addsuffix .$(name), $(MODNAMES))))
$(foreach name, clean lib install , $(eval .PHONY : $(name)_targets))

############################################################################################################
################################### BUILD THE LIBRARY ######################################################

.PHONY: help build install install.deps clean  
help: 
	@echo " The makefile has targets: \n" \
				" 1. help: this message \n" \
				" 2. build: build the library \n" \
				" 3. install: build and install the library \n" \
				" 4: install.deps downloads and installs zoltan, sqlite3 and waffles in the XXX_DIR directories \n" \
				" 5: clean: removes all the object files and the install directory \n" \



build: $(lib_targets)

install: preinstall $(install_targets) 
	@$(foreach ob,  $(wildcard $(INSTALL_DIR)/include/*.h), $(file >>$(INSTALL_DIR)/include/api.h, #include "$(notdir $(ob))") ) 

install.deps: download install_deps  

clean: $(clean_targets) 
	@rm -f $(INSTALL_DIR)/include/*.h
	@rm -f $(INSTALL_DIR)/lib/lib$(SSLIBNAME).a


###########################################################################################################
###########################################################################################################

#### Build the individual object files 
%.o : %.C
	$(SS) -c $(SS_FLAGS) $< -o $@



%.ss : %.C 
	$(SS) -I $(INSTALL_DIR)/include/ $(SS_FLAGS) -o $@ $< $(SS_LFLAGS) 

petsc_ex: $(SRC_DIR)/examples/petsc_example/ex1.ss

dummy_ex: $(SRC_DIR)/examples/dummy/dummy.ss
	

$(foreach i, $(MODNAMES), $(eval $(i)_list:=$(patsubst %.C, %.o, $(wildcard $(SRC_DIR)$(i)/src/*.C))))  
$(foreach i, $(MODNAMES), $(eval .PHONY : $(i).object) ${\n} $(eval $(i).object : $($(i)_list)))

$(clean_targets): 
	@rm -f $(wildcard $(SRC_DIR)$(basename $@)/src/*.o)

$(lib_targets): %.lib : %.object
	@ar -rcs lib$(SSLIBNAME).a $(wildcard $(SRC_DIR)$(basename $@)/src/*.o)	

$(install_targets): %.install : %.lib
	@install $(SRC_DIR)lib$(SSLIBNAME).a $(INSTALL_DIR)/lib/
	@install $(SRC_DIR)$(basename $@)/include/*.h $(INSTALL_DIR)/include/
	@for f in $(wildcard $(SRC_DIR)$(basename $@)/include/*.T); do cp $$f $(INSTALL_DIR)/include/ ; done; 

preinstall:
	@mkdir -p $(INSTALL_DIR)/include
	@mkdir -p $(INSTALL_DIR)/lib
	@rm -f $(INSTALL_DIR)/include/api.h
	
install_deps: install.waffles install.sqlite3 install.zoltan

download: download.waffles download.sqlite3 download.zoltan

download.waffles: 
	mkdir -p $(EXTRA_DIR)
	cd $(EXTRA_DIR) ; if [ ! -f $(notdir $(WAFFLES_DOWNLOAD)) ]; then wget $(WAFFLES_DOWNLOAD); tar -xf $(notdir $(WAFFLES_DOWNLOAD)) ; fi ;
	
download.sqlite3: 
	mkdir -p $(EXTRA_DIR)
	cd $(EXTRA_DIR) ; if [ ! -f $(notdir $(SQLITE3_DOWNLOAD)) ] ; then wget $(SQLITE3_DOWNLOAD); tar -xf $(notdir $(SQLITE3_DOWNLOAD)) ; fi;

download.zoltan:
	mkdir -p $(EXTRA_DIR)
	$(info You should prob go to the zoltan download page register, just to be nice )
	cd $(EXTRA_DIR) ; if [ ! -f $(notdir $(ZOLTAN_DOWNLOAD)) ] ; then wget $(ZOLTAN_DOWNLOAD) ; tar -xf $(notdir $(ZOLTAN_DOWNLOAD)) ; fi;

install.sqlite3:
	cd $(EXTRA_DIR)/$(basename $(basename $(notdir $(SQLITE3_DOWNLOAD)))) ; \
	./configure --prefix=$(abspath $(SQLITE3_DIR)) ; make ; make install ;	 

install.waffles:
	mkdir -p $(WAFFLES_DIR)/lib
	mkdir -p $(WAFFLES_DIR)/include/GClasses
	cd $(EXTRA_DIR)/waffles-1.0.0/src/GClasses ; sed -e s/-Werror//g -i.backup Makefile ; make opt ;
	install $(EXTRA_DIR)/waffles-1.0.0/lib/libGClasses.a $(WAFFLES_DIR)/lib
	install $(EXTRA_DIR)/waffles-1.0.0/src/GClasses/*.h $(WAFFLES_DIR)/include/GClasses/ 

install.zoltan:
	cd $(EXTRA_DIR)/Zoltan_v3.83/ ; mkdir -p build ; cd build ; \
	../configure --prefix=$(abspath $(ZOLTAN_DIR)) ; make ; make install ;
	
 			
