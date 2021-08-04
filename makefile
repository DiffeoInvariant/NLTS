PETSC_DIR=/usr/local/petsc
PETSC_ARCH=arch-linux-c-opt
CXX=clang++
CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -fPIC -g
LDFLAGS=-shared $(PETSC_WITH_EXTERNAL_LIB) 

include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules

NLTS_INCL=include/
INCLS=-I$(NLTS_INCL) -I./packages/GSL/include $(PETSC_CC_INCLUDES)
NLTS_SRC_DIR=src/
NLTS_SRCS:=$(shell find $(NLTS_SRC_DIR) -name '*.cpp')

TARGET_DIR=lib
TARGET=$(TARGET_DIR)/libnlts.so
NLTS_PROGRAMS_DIR=examples

.PHONY: all $(TARGET) mutual mutual-opt

all: $(TARGET) mutual mutual-opt

mutual-opt:$(NLTS_PROGRAMS_DIR)/mutual_opt.cpp
	$(CXX) $(CXXFLAGS) $^ -o bin/mutual-opt $(INCLS) $(PETSC_WITH_EXTERNAL_LIB) -L./lib -lnlts -Wl,-rpath,$(shell pwd)/lib

mutual: $(NLTS_PROGRAMS_DIR)/mutual.cpp
	$(CXX) $(CXXFLAGS) $^ -o bin/mutual $(INCLS) $(PETSC_WITH_EXTERNAL_LIB) -L./lib -lnlts -Wl,-rpath,$(shell pwd)/lib

$(TARGET): $(NLTS_SRCS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(LDFLAGS)
