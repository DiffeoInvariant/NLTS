# NLTS
A C++ library for non-linear time series, built on PETSc.

# Installation
```
git clone git@github.com:DiffeoInvariant/NLTS.git
```
## Dependencies
  You need working installations of PETSc (can be downloaded via `configure`) and boost. 
  You can download PETSc from https://www.mcs.anl.gov/petsc/download/index.html.
  
## Configuration
  Run `./configure -h` to see configure options.
  
  
## Building
```
cd NLTS
make all
```

# TISEAN Programs
NLTS contains several (currently, only one--`mutual`) implementations of programs included in the TISEAN library.
These can be found in `NLTS/bin` (implementations are contained in `examples`), and have very similar syntax to TISEAN programs. Pass the `-h` option to any of these programs
to see help options.

# Examples
Examples and tests can be found in the so-named directories. For example, to run the file `tests/test_diff.cpp`, from the root of the NLTS directory, do
```
cd tests
make all
./txt2petsc
./test_diff
```
The program `txt2petsc` converts the supplied example data into PETSc binary files that can be read in by the other test files,
and only needs to be run once, before running any of the other tests.
