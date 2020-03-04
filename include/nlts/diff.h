#ifndef NLTS_DIFF_H
#define NLTS_DIFF_H
#include <petscvec.h>

namespace nlts
{

  enum StencilType
  {
   Forward,
   Backward,
   Central
  };
  
  extern PetscErrorCode FiniteDifference(Vec x, Vec dxdt, PetscReal dt, PetscInt order,
					 StencilType direction, bool periodic_domain=false);






}

#endif
