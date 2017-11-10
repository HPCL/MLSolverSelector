
include ../makefile.inc 
$(info sdfsdfsdsdfsdfsdf$(abspath $(INSTALL_DIR)/include))

.PHONY: clean all petscex1 dummy

all: petscex1 petscex2 dummy
petscex1: 
	$(SS) -I $(INSTALL_DIR)/include $(SS_FLAGS) -o petscex1.ss petscex1.C $(SS_LFLAGS) 

petscex2: 
	$(SS) -I $(INSTALL_DIR)/include $(SS_FLAGS) -o petscex2.ss petscex2.C $(SS_LFLAGS) 

dummy: 
	$(SS) -I $(INSTALL_DIR)/include $(SS_FLAGS) -o dummy.ss dummy.C $(SS_LFLAGS) 

clean:
	rm petscex1.ss dummy.ss
	

### Callgrind profiles the example. gprof2dot changes output to a flowchart, which is nice. eog
### is the gnome image viewe
callgrind: 
	rm -f callgrind_out
	valgrind --tool=callgrind --callgrind-out-file=callgrind_out ./petscex2.ss 
	python ../gprof2dot/build/lib/gprof2dot.py -e 0.2 -n 0.2 -f callgrind callgrind_out | dot -Tsvg -o output.svg
	eog output.svg

#mem checker
valgrind:
	valgrind --tool=memcheck --leak-check=full ./petscex2.ss


