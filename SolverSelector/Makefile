
include Makefile.inc

###########################################################################################################
###########################################################################################################
##################################### Install the dependencies ##########################################
.PHONY: all examples 
all: examples 

petsc_example:
	$(MAKE) -C petsc_example 

#install and download dependencies
install.deps: download install_deps  

#download the dependencies
download: download.waffles download.sqlite3 

download.waffles: 
	mkdir -p $(EXTRA_DIR)
	cd $(EXTRA_DIR) ; if [ ! -f $(notdir $(WAFFLES_DOWNLOAD)) ]; then wget $(WAFFLES_DOWNLOAD); tar -xf $(notdir $(WAFFLES_DOWNLOAD)) ; fi ;
	
download.sqlite3: 
	mkdir -p $(EXTRA_DIR)
	cd $(EXTRA_DIR) ; if [ ! -f $(notdir $(SQLITE3_DOWNLOAD)) ] ; then wget $(SQLITE3_DOWNLOAD); tar -xf $(notdir $(SQLITE3_DOWNLOAD)) ; fi;


### install the dependencies
install_deps: install.waffles install.sqlite3 

install.sqlite3:
	cd $(EXTRA_DIR)/$(basename $(basename $(notdir $(SQLITE3_DOWNLOAD)))) ; \
	./configure --prefix=$(abspath $(SQLITE3_DIR)) ; make ; make install ;	 

install.waffles:
	mkdir -p $(WAFFLES_DIR)/lib
	mkdir -p $(WAFFLES_DIR)/include/GClasses
	cd $(EXTRA_DIR)/waffles-1.0.0/src/GClasses ; sed -e s/-Werror//g -i.backup Makefile ; make opt ;
	install $(EXTRA_DIR)/waffles-1.0.0/lib/libGClasses.a $(WAFFLES_DIR)/lib
	install $(EXTRA_DIR)/waffles-1.0.0/src/GClasses/*.h $(WAFFLES_DIR)/include/GClasses/ 

	

