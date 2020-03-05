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

    std::pair<std::vector<double>, std::vector<double>> read_scalar_trajectory(std::string);

    PetscErrorCode VecReadScalarTrajectory(std::string, Vec *, Vec *);

    PetscErrorCode VecReadBinary(std::string, Vec *);

    PetscErrorCode VecWriteScalarTrajectory(std::string filename_pre, Vec X, Vec T, bool make_new_directory=true);

    extern PetscErrorCode Initialize(int *argc, char ***argv,
			    const char file[]=NULL,
			    const char help[]=NULL);

    extern PetscErrorCode Finalize();

    extern bool has_petsc_option(std::string name,
				 std::optional<std::string> prepend=std::nullopt,
				 PetscOptions opts_db=NULL);

    template<typename T>
    extern std::pair<T, bool> get_petsc_option(std::string name,
					       std::optional<std::string> prepend=std::nullopt,
					       PetscOptions opts_db=NULL);
  }
}
#endif
