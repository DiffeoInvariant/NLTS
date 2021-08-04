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

  /* param 1 (input): the vector containing the function values to difference
     param 2 (output): the vector containing the differences. MUST be the same
     local size as the first parameter; the first and last order-many elements
     (e.g. first and last 2 for a second-order scheme) are computed with 
     decreasing-order (and possibly a different type than specified)
     schemes as fewer elements become available to difference, so use these
     at your own peril.
     param 3 (input): the denominator of the finite difference (up to a
     scheme-dependent multiplicative constant and power); if x = f(t),
     this is the difference between t points evaluated in the vector x (the 
     first parameter).
     param 4: The order of the difference scheme in dt. For Central differences,
     this must be even, and <= 8; for Forward and Backward differences, this 
     must be <= 6
     param 5: The direction of the difference scheme; for elements within 
     `order`-many of the beginning or end of dxdt however, a scheme with a 
     different direction (and order) than the one specified may be used*/
  extern PetscErrorCode FiniteDifference(Vec x, Vec dxdt, PetscReal dt, PetscInt order,
					 StencilType direction);






}

#endif
