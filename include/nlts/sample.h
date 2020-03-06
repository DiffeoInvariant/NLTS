#ifndef NLTS_SAMPLE_H
#define NLTS_SAMPLE_H

#include <petscvec.h>
namespace nlts
{

  extern PetscErrorCode DecimatingDownsample(Vec X, Vec *Y, PetscInt M);

}
#endif
