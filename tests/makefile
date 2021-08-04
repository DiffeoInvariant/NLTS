
default: all
CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -g $(INCLS)
LDFLAGS=$(PETSC_WITH_EXTERNAL_LIB) -Wl,-rpath,../lib -L../lib -lnlts
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules

NLTS_INCL=../include/
ODEINT_INCL=../packages/odeint/include/

INCLS=-I$(NLTS_INCL) -I$(ODEINT_INCL) -I../packages/boost_1_66_0/include -I../packages/GSL/include $(PETSC_CC_INCLUDES)

TARGET_IO = test_io

.PHONY: all txt2petsc io in diff embed nndb test_bcd

all: txt2petsc io in diff embed nndb test_bcd


embed: test_embedding.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_embed $(INCLS) $(LDFLAGS)
	./test_embed

diff: test_diff.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_diff $(INCLS) $(LDFLAGS)
	./test_diff
in: test_in.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_input $(INCLS) $(LDFLAGS)
	./test_input

txt2petsc: txt2petsc.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(LDFLAGS)
	./txt2petsc

io: test_io.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_io $(INCLS) $(LDFLAGS)
	./test_io

test_bcd: test_bcd.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(LDFLAGS) 
	./test_bcd

clean:
	@$(RM) *.o
