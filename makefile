
OBJECTS = mesh.o fem.o

include $(SLEPC_DIR)/conf/slepc_common

.SECONDARY:
.PHONY: clobber

CXXFLAGS += -std=c++11
CFLAGS += -std=c99

%.exe: $(OBJECTS) %.o
	$(CLINKER) $*.o -o $@ $(OBJECTS) $(SLEPC_EPS_LIB)

clobber:
	rm -rf *.o *.so *.exe
