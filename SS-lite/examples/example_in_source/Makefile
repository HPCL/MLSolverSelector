


include ../../Makefile.inc

### Get all the src files in the current directory
BASE_SRC=$(wildcard ./*.cpp)

### create list of all source files with a .o extension  
BASE_OBS=$(BASE_SRC:%.cpp=%.ss)

# The default rule is to build all .cpp files in this directory.
default: $(BASE_OBS)

%.ss: %.cpp
	$(SS) $(SS_FLAGS) -I ../../include/ -I ../../feature_sets/include/ $(SS_DEP_FLAGS) $^ $(SS_DEP_LFLAGS) -L../.. -lsslite -o $@

clean:
	rm -f $(BASE_OBS) 
