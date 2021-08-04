
CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -g
LDFLAGS=$(PETSC_WITH_EXTERNAL_LIB) -Wl,-rpath,../lib -L../lib -lnlts
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules

NLTS_INCL=../include/
ODEINT_INCL=/home/diffeoinvariant/odeint/include/

INCLS=-I$(NLTS_INCL) -I$(ODEINT_INCL) -I/usr/include/boost -I../packages/GSL/include $(PETSC_CC_INCLUDES)

BOOST_LINK = -Wl,-rpath,../packages/boost_1_66_0/lib -L../packages/boost_1_66_0/lib -lboost_system -lboost_filesystem

TARGET_IO = test_io

.PHONY: all io txt2petsc in diff embed nndb test_bcd

all: io txt2petsc in diff embed nndb test_bcd

nndb: test_nndb.cpp
	$(CXX) $(CXXFLAGS) $^ -o test_nndb $(INCLS) $(LDFLAGS) 

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

test_bcd: test_bcd.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(LDFLAGS) 

clean:
	@$(RM) *.o
