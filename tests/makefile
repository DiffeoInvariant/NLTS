PETSC_DIR=/usr/local/petsc
PETSC_ARCH=arch-linux-cxx-debug
CXX=clang++
CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -g
LDFLAGS=$(PETSC_WITH_EXTERNAL_LIB) -L../lib -lnlts 
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules

NLTS_INCL=../include/
INCLS=-I$(NLTS_INCL) -I/usr/include/boost $(PETSC_CC_INCLUDES) 

TARGET_IO = test_io

.PHONY: all io txt2petsc in diff embed

all: io txt2petsc in diff embed

embed: test_embedding.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_embed $(INCLS) $(LDFLAGS)

diff: test_diff.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_diff $(INCLS) $(LDFLAGS)
in: test_in.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_input $(INCLS) $(LDFLAGS)

txt2petsc: txt2petsc.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(LDFLAGS)

io: test_io.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_io $(INCLS) $(LDFLAGS)
clean:
	@$(RM) *.o
