
PETSC_DIR=/usr/local/petsc
PETSC_ARCH=arch-linux-cxx-debug
CXX=clang++
CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -fPIC -g 
LDFLAGS=-shared $(PETSC_WITH_EXTERNAL_LIB) -L/lib/x86_64-linux-gnu -lboost_filesystem
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules

NLTS_INCL=include/
INCLS=-I$(NLTS_INCL) -I/usr/include/boost $(PETSC_CC_INCLUDES) 

NLTS_SRC_DIR=src/
NLTS_SRCS:=$(shell find $(NLTS_SRC_DIR) -name '*.cpp')


TARGET_DIR=lib
TARGET=$(TARGET_DIR)/libnlts.so

NLTS_PROGRAMS_DIR=examples

.PHONY: all $(TARGET) mutual

all: $(TARGET) mutual

mutual: $(NLTS_PROGRAMS_DIR)/mutual.cpp
	$(CXX) $(CXXFLAGS) $^ -o bin/mutual $(INCLS) $(PETSC_WITH_EXTERNAL_LIB) -L./lib -lnlts
	LD_LIBRARY_PATH=lib

$(TARGET): $(NLTS_SRCS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(LDFLAGS)
	LD_LIBRARY_PATH=lib

clean:
	@$(RM) $(TARGET_DIR)
