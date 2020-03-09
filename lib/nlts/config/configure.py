#!/usr/bin/python3
from __future__ import print_function
import os, sys, re
import argparse
from urllib.request import urlretrieve
'''
Based on, and partly copied from, PETSc's configure.py
'''
PETSC_DOWNLOAD_URL='http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.12.4.tar.gz'
DEFAULT_PETSC_ARCH='arch-libnlts-cxx-debug'
DEFAULT_PETSC_DIR='packages/petsc'
def check_python_version():
    if sys.version_info < (2,6):
        print('************************************************************************')
        print('*      Python version 2.6+ or 3.4+ is required to run ./configure      *')
        print('*         Try: "python2.7 ./configure" or "python3 ./configure"        *')
        print('************************************************************************')
        sys.exit(4)
        
def download_petsc():
    print('************************************************************************')
    print('*              Attempting to download PETSc from                       *')
    print('*http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.12.4.tar.gz *')
    print('************************************************************************')
    os.system('mkdir packages')
    petsc_tarball, headers = urlretrieve(PETSC_DOWNLOAD_URL, 'packages/petsc.tar.gz')
    os.system('tar -xzf packages/petsc.tar.gz')

def configure_and_build_petsc(CC, CXX):
    config_str = f"./configure PETSC_ARCH={DEFAULT_PETSC_ARCH} --with-cc={CC} --with-cxx={CXX} --download-mpich --download-sowing --download-c2html --download-zlib --download-libpng --download-fftw --with-fc=0 COPTFLAGS='-O3 -march=native -mtune=native -fPIC' CXXOPTFLAGS='-O3 -march=native -mtune=native -fPIC'"

    os.system(f"cd {DEFAULT_PETSC_DIR} && {config_str}")
    os.system('make PETSC_DIR={DEFAULT_PETSC_DIR} PETSC_ARCH={DEFAULT_PETSC_ARCH} all')
    os.system('make PETSC_DIR={DEFAULT_PETSC_DIR} PETSC_ARCH={DEFAULT_PETSC_ARCH} check')
    os.system('cd -')
                            
    
def build_option(lst, pre, name, helpstr=""):
    lst.append((pre + name, helpstr))
def build_options_and_downloads(options_pre="--with_", download_pre="--download_"):
    options_list = []
    download_list = []
    build_option(options_list, options_pre, "PETSC_DIR", "Directory with your PETSc installation.")
    build_option(options_list, options_pre, "PETSC_ARCH", "PETSC_ARCH for your PETSc installation.")
    build_option(download_list, download_pre, "petsc", "Download PETSc")
    build_option(options_list, options_pre, "boost_dir", "Directory for your Boost installation.")
    build_option(options_list, options_pre, "boost_lib_dir", "Directory of your boost/filesystem library.")
    #build_option(download_list, download_pre, "boost", "Download Boost")
    build_option(options_list, options_pre, "cc", "Which C compiler to use? Default clang.")
    build_option(options_list, options_pre, "cxx", "Which C++ compiler to use? Default clang++. If you provide a different cxx, it MUST support the C++17 standard.")
    
    return options_list, download_list


def build_parser(opts, download_opts):
    parser = argparse.ArgumentParser()
    for op in opts:
        parser.add_argument(op[0], help=op[1])
    
    for dop in download_opts:
        parser.add_argument(dop[0], help=dop[1], action="store_true")

    return parser

def assign_if(options, arg, name, default):
    options[name] = arg if arg else default
    return options
        
def parse_arguments(parser):
    args = parser.parse_args()
    options = {}

    options = assign_if(options, args.with_PETSC_DIR, 'PETSC_DIR', DEFAULT_PETSC_DIR)
    options = assign_if(options, args.with_PETSC_ARCH, 'PETSC_ARCH', DEFAULT_PETSC_ARCH)
    options = assign_if(options, args.with_boost_dir, 'boost_dir', '/usr/include/boost')
    options = assign_if(options, args.with_boost_lib_dir, 'boost_lib', '/lib/x86_64-linux-gnu')
    options = assign_if(options, args.with_cc, 'CC', 'clang')
    options = assign_if(options, args.with_cxx, 'CXX', 'clang++')
    options['download_petsc'] = args.download_petsc
    
    return options

def handle_dependencies(options):
    if options['download_petsc']:
        download_petsc()
        configure_and_build_petsc(options['CC'], options['CXX'])


def build_makefile_lines(options):
    lines = [f"PETSC_DIR={options['PETSC_DIR']}\n",
             f"PETSC_ARCH={options['PETSC_ARCH']}\n",
             f"CXX={options['CXX']}\n",
             f"CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -fPIC -g\n",
             f"LDFLAGS=-shared $(PETSC_WITH_EXTERNAL_LIB) -L{options['boost_lib']} -lboost_filesystem\n\n",
             "include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables\n",
             "include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules\n\n",
             "NLTS_INCL=include/\n",
             f"INCLS=-I$(NLTS_INCL) -I{options['boost_dir']} $(PETSC_CC_INCLUDES)\n",
             "NLTS_SRC_DIR=src/\n",
             "NLTS_SRCS:=$(shell find $(NLTS_SRC_DIR) -name '*.cpp')\n\n",
             "TARGET_DIR=lib\n",
             "TARGET=$(TARGET_DIR)/libnlts.so\n",
             "NLTS_PROGRAMS_DIR=examples\n\n",
             ".PHONY: all $(TARGET) mutual mutual-opt\n\n",
             "all: $(TARGET) mutual mutual-opt\n\n",
             "mutual-opt:$(NLTS_PROGRAMS_DIR)/mutual_opt.cpp\n",
             "\t$(CXX) $(CXXFLAGS) $^ -o bin/mutual-opt $(INCLS) $(PETSC_WITH_EXTERNAL_LIB) -L./lib -lnlts\n",
             "\tLD_LIBRARY_PATH=lib\n\n",
             "mutual: $(NLTS_PROGRAMS_DIR)/mutual.cpp\n",
             "\t$(CXX) $(CXXFLAGS) $^ -o bin/mutual $(INCLS) $(PETSC_WITH_EXTERNAL_LIB) -L./lib -lnlts\n\n",
             "$(TARGET): $(NLTS_SRCS)\n",
             "\t$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(LDFLAGS)\n",
             "\tLD_LIBRARY_PATH=lib\n\n"]

    return lines

def write_makefile(lines):
    with open('makefile', 'w') as mf:
        for ln in lines:
            mf.write(ln)
    
             


if __name__ == '__main__':
    check_python_version()
    options, download_opts = build_options_and_downloads()
    parser = build_parser(options, download_opts)
    options = parse_arguments(parser)
    handle_dependencies(options)
    makefile_lines = build_makefile_lines(options)
    write_makefile(makefile_lines)
