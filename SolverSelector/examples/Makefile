
SS_INSTALL_DIR=/home/jrobinson/Document/solverselecter
include $(SS_INSTALL_DIR)/Makefile.inc 

### Get all the src files
BASE_SRC=$(wildcard ./*.cpp)

### Get all the object files 
BASE_OBS=$(BASE_SRC:%.cpp=%.ss)

default: $(BASE_OBS)

%.ss: %.cpp
	$(SS) $(SS_FLAGS) $(SS_INCLUDE) -g $^ $(SS_LIBS) -o $@ -I/home/jrobinson/Documents/papi-5.6.0/src/testlib -I/home/jrobinson/Documents/papi-5.6.0/src/ -I/home/jrobinson/Documents/papi-5.6.0/src/components/rapl/tests -L/home/jrobinson/Documents/papi-5.6.0/src/ -lpapi

clean:
	rm -f $(BASE_OBS) 
