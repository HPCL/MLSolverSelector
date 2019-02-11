
include Makefile.inc

.PHONY: clean
.PHONY: install
.PHONY: feature_set_generation 

OUT=lib$(SS_LIB).a
BASE_SRC=$(wildcard src/*.cpp)
BASE_SRC+=$(wildcard feature_sets/src/*.cpp)

BASE_OBS=$(BASE_SRC:%.cpp=%.o)

BASE_OBS=$(BASE_SRC:%.cpp=%.o)

SS_INCS = -I include/ -I feature_sets/include  

SS_FLAGS +=-Wno-write-strings 

default: $(OUT)

feature_set_generation: 
	cd feature_sets && ./generateExtraction.sh 

$(OUT): feature_set_generation $(BASE_OBS) 
	ar rcs $@ $(BASE_OBS) $(IMPL_OBS) 

.cpp.o:
	$(SS) $(SS_FLAGS) $(SS_INCS) $(SS_DEP_FLAGS) -c $< -o $@

clean:
	rm -f $(BASE_OBS) $(OUT) 


#Install the library, includes and a copy of the makefile.inc file.
#Using that, you should be able to compile any Petsc application 
install:
	mkdir -p $(SS_INSTALL_DIR)/lib
	mkdir -p $(SS_INSTALL_DIR)/include
	install $(OUT) $(SS_INSTALL_DIR)/lib
	install include/*h feature_sets/include/*.h $(SS_INSTALL_DIR)/include
	install Makefile.inc $(SS_INSTALL_DIR) 
	

