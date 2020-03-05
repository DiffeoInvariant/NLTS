#ifndef NLTS_SYS_H
#define NLTS_SYS_H
#include <petscsys.h>
#include <petscvec.h>
#include <vector>

namespace nlts
{
  extern PetscErrorCode VecSetFromStd(Vec *, const std::vector<PetscScalar>&);

}
#endif
