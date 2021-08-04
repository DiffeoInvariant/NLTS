#ifndef NLTS_IO_H
#define NLTS_IO_H
#include <petscsys.h>
#include <petscvec.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>
#include <optional>
namespace nlts
{
  inline namespace io
  {
    /* param 1 (input): filename of a file formatted like

       x0 t0
       x1 t1
       ...
       xn tn
       
       with ti and xi being times and positions
       returns: pair of std::vectors, the first containing the positions, the second 
       containing the times
    */
    extern std::pair<std::vector<double>, std::vector<double>> read_scalar_trajectory(std::string);

    /* param 1 (input): filename of a file formatted like 

       x0 t0
       x1 t1
       ...
       xn tn

       param 2 (output): Pointer to a Vec for the positions xi
       param 3 (output): Pointer to a Vec for the times ti
       returns: an error code
    */
    extern PetscErrorCode VecReadScalarTrajectory(std::string, Vec *, Vec *);

    /* param 1 (input): filename of a file in PETSc binary format
       param 2 (output): pointer to a Vec to read that file into
       returns: an error code
    */
    extern PetscErrorCode VecReadBinary(std::string, Vec *);

    /* param 1 (input): prefix to the filenames that will be written. The 
       written files will be filename_pre + "_x.dat" and filename_pre + "_t.dat"
       param 2 (input): the spatial component of the trajectory
       param 3 (input): the time component of the trajectory
       param 4 (input, optional): create a new directory if the one specified 
       by filename_pre doesn't exist? defaults to true
       returns: an error code
    */
    extern PetscErrorCode VecWriteScalarTrajectory(std::string filename_pre, Vec X, Vec T, bool make_new_directory=true);

    /* for now, NLTS doesn't need any special initialization or finalization
       other than what PETSc provides, so we will provide those functions in 
       the nlts namespace */
    extern PetscErrorCode Initialize(int *argc, char ***argv,
			    const char file[]=NULL,
			    const char help[]=NULL);

    extern PetscErrorCode Finalize();

    /* param 1 (input): the name of the PETSc option
       param 2 (input, optional): prefix of the option, defaults to none
       param 3 (input, optional): pointer to the PETSc options database, if 
       any exists. defaults to NULL, the default options database extracted
       from argv 
       returns: does the option exist in the database?
    */
    extern bool has_petsc_option(std::string name,
				 std::optional<std::string> prepend=std::nullopt,
				 PetscOptions opts_db=NULL);

    /*
      template param T: the type of the PETSc option, e.g. PetscInt, 
      PetscReal, std::string, etc 

      param 1 (input): the name of the PETSc option
      param 2 (input, optional): prefix of the option, defaults to none
      param 3 (input, optional): pointer to the PETSc options database, if 
      any exists. defaults to NULL, the default options database extracted
      from argv 
      returns: does the option exist in the database?
    */
    template<typename T>
    extern std::pair<T, bool> get_petsc_option(std::string name,
					       std::optional<std::string> prepend=std::nullopt,
					       PetscOptions opts_db=NULL);
  }
}
#endif
