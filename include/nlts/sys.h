#ifndef NLTS_SYS_H
#define NLTS_SYS_H
#include <petscsys.h>
#include <petscvec.h>
#include <vector>
#include <gsl/span>

namespace nlts
{
  extern PetscErrorCode VecSetFromStd(Vec *, const std::vector<PetscScalar>&);

  extern PetscErrorCode VecGetSubRange(Vec, PetscInt, PetscInt, Vec *);

}
#endif
