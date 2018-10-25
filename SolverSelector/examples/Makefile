
SS_INSTALL_DIR=/usr/local/solverselecter
include $(SS_INSTALL_DIR)/Makefile.inc 

### Get all the src files
BASE_SRC=$(wildcard ./*.cpp)

### Get all the object files 
BASE_OBS=$(BASE_SRC:%.cpp=%.ss)

default: $(BASE_OBS)

%.ss: %.cpp
	$(SS) $(SS_FLAGS) $(SS_INCLUDE) $^ $(SS_LIBS) -o $@

clean:
	rm -f $(BASE_OBS) 
