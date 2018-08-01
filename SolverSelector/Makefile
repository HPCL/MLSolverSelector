
include Makefile.inc

OUT=lib$(SS_LIB).a

### Get all the src files
BASE_SRC=$(wildcard base/src/*.cpp)
IMPL_SRC=$(wildcard impls/src/*.cpp)

### Get all the object files 
BASE_OBS=$(BASE_SRC:%.cpp=%.o)
IMPL_OBS=$(IMPL_SRC:%.cpp=%.o)

SS_INCS = -I base/include -I base/template -I impls/include 

default: $(OUT)

$(OUT): $(BASE_OBS) $(IMPL_OBS)
	ar rcs $@ $(BASE_OBS) $(IMPL_OBS) 

.cpp.o:
	$(SS) $(SS_FLAGS) $(SS_INCS) $(SS_DEP_FLAGS) -c $< -o $@

clean:
	rm -f $(BASE_OBS) $(IMPL_OBS) $(OUT) $(BASE_SRC:%.cpp=%.d) $(IMPL_SRC:%.cpp=%.d) Makefile.bak

distclean: check_clean clean 
	rm -r $(SS_INSTALL_DIR)
	
check_clean: 
	@( read -p "This will remove $(SS_INSTALL_DIR) including everything installed in there. Are you sure?!? [y/N]: " sure && case "$$sure" in [yY]) true;; *) false;; esac )

.PHONY: check_clean clean

install:
	mkdir -p $(SS_INSTALL_DIR)/lib
	mkdir -p $(SS_INSTALL_DIR)/include
	
	install $(OUT) $(SS_INSTALL_DIR)/lib
	install base/include/*h base/template/*.tpp impls/include/*.h $(SS_INSTALL_DIR)/include
	install Makefile.inc $(SS_INSTALL_DIR) 
	install Makefile_Moose.inc $(SS_INSTALL_DIR) 
	
###########################################################################################################
###########################################################################################################
##################################### Install the dependencies ##########################################

### install the dependencies
build_deps: sqlite3.build waffles.build 
install_deps: sqlite3.install waffles.install
clean_deps: waffles.clean sqlite3.clean 

sqlite3.build: 
	cd ./3rd_party/ ; \
	tar -xf sqlite-autoconf-3210000.tar.gz ; \
	cd sqlite-autoconf-3210000 ; \
	./configure --prefix=$(SQLITE3_DIR) ; make ; 	 

waffles.build: 
	cd ./3rd_party ; tar -xf 1.0.0.tar.gz; cd waffles-1.0.0/src/GClasses ; sed -e s/-Werror//g -i.backup Makefile ; make opt ; 

sqlite3.install: 
	mkdir -p $(SQLITE3_DIR)
	cd ./3rd_party/sqlite-autoconf-3210000 ; make install ; 

waffles.install:
	mkdir -p $(WAFFLES_DIR)/lib
	mkdir -p $(WAFFLES_DIR)/include/GClasses
	install ./3rd_party/waffles-1.0.0/lib/libGClasses.a $(WAFFLES_DIR)/lib ;\
	install ./3rd_party/waffles-1.0.0/src/GClasses/*.h $(WAFFLES_DIR)/include/GClasses/ ;

# Instead of cleaning, we just delete the entire 3rd party directories ( these are not needed once installed ) 
sqlite3.clean: 
	rm -r  ./3rd_party/sqlite-autoconf-3210000 ;  
waffles.clean: 
	rm -r ./3rd_party/waffles-1.0.0/;  

