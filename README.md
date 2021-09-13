# NLTS
A C++ library for non-linear time series, built on PETSc.

# Installation
```
git clone git@github.com:DiffeoInvariant/NLTS.git
```
## Dependencies
  You need working installations of PETSc (can be downloaded via `configure`) and boost. 
  You can download PETSc from https://www.mcs.anl.gov/petsc/download/index.html, or have NLTS do it for you by passing the `--download_petsc` flag to `configure`. For plotting to work (and to use the default PETSc downloaded and configured with `--download_petsc`) you will need X windows headers, which you can
  install on ubuntu with `sudo apt install libx11-dev`, or on Arch with `pacman -S libx11-devel`.
  
## Configuration
  Run `./configure -h` to see configure options.
  
  
## Building
```
./configure <options>
make
```
where `<options>` are any options you want to pass to `configure`, which you can view by running `./configure -h`

# TISEAN Programs
NLTS contains several (currently, only one--`mutual`) implementations of programs included in the [https://www.pks.mpg.de/tisean/](TISEAN) library.
These can be found in `NLTS/bin` (implementations are contained in `examples`), and have very similar syntax to TISEAN programs. Pass the `-h` option to any of these programs
to see help options.

NLTS also contains `mutual-opt`, which is like TISEAN's `mutual` except it optimizes the lag parameter `tau` to minimize mutual information.

# Examples
Examples and tests can be found in the so-named directories. For example, to run the file `tests/test_diff.cpp`, from the root of the NLTS directory, do
```
cd tests
make all
./txt2petsc
./test_diff
```
The program `txt2petsc` converts the supplied example data into PETSc binary files that can be read in by the other test files,
and only needs to be run once, before running any of the other tests. If you want to run the box-counting dimension test (`test_bcd`), you will need to run `git clone https://github.com/diffeoinvariant/odeint.git` to install the `odeint` library to generate a test trajectory, and you will also need to install Boost version 1.66.0 in `packages/boost_1_66_0` (no need to build any of the Boost libraries) for access to `boost::multiprecision` to compare actual memory usage with that from a naive algorithm for the box-counting dimension.
