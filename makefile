
OBJECTS = mesh.o fem.o

include $(SLEPC_DIR)/conf/slepc_common

.SECONDARY:
.PHONY: clobber

CXXFLAGS += -std=c++11
CFLAGS += -std=c99

%.exe: $(OBJECTS) %.o
	$(CLINKER) $*.o -o $@ $(OBJECTS) $(SLEPC_EPS_LIB)

cgal_experiment: cgal_experiment.cpp
	$(CXX) -o cgal_experiment cgal_experiment.cpp -lCGAL -lCGAL_Core -lgmp -lmpfr -fround-math

clobber:
	rm -rf *.o
