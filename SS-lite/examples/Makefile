
#### This is an example of a makefile so building an executable using 
#### SS-lite. Because most people will be using this out-of-source, this
#### makefile assumes the library is installed at 

SS_INSTALL_DIR=/usr/local/solverselecter-lite

#### So, if you changed it, make sure that you change it here too. 


include $(SS_INSTALL_DIR)/Makefile.inc 

### Get all the src files in the current directory
BASE_SRC=$(wildcard ./*.cpp)

### create list of all source files with a .o extension  
BASE_OBS=$(BASE_SRC:%.cpp=%.ss)

# The default rule is to build all .cpp files in this directory.
default: $(BASE_OBS)

%.ss: %.cpp
	$(SS) $(SS_FLAGS) $(SS_INCLUDE) $^ $(SS_LIBS) -o $@

clean:
	rm -f $(BASE_OBS) 
