#include <nlts/diff.h>

namespace nlts
{


  static PetscErrorCode BackwardFd(const PetscScalar *x, PetscScalar *dxdt, PetscInt idx, PetscReal dt, PetscInt order)
  {
    PetscFunctionBeginUser;
    if(order == 0 or order > 6){
      SETERRQ(PETSC_COMM_WORLD, 1, "Must supply backward finite difference order greater than 0 but less than 7.");
    }
    switch(order){
    case 1:
      dxdt[idx] = (-x[idx - 1] + x[idx]) / dt;
      break;
    case 2:
      dxdt[idx] = (0.5 * x[idx - 2] - 2.0 * x[idx - 1] + 1.5 * x[idx]) / dt;
      break;
    case 3:
      dxdt[idx] = (-x[idx + 3] / 3.0 + 1.5 * x[idx + 2] - 3.0 * x[idx + 1] + 11.0 * x[idx] / 6.0) / dt;
      break;
    case 4:
      dxdt[idx] = (x[idx - 4]/4.0 - 4 * x[idx - 3] / 3.0 + 3.0 * x[idx - 2] - 4.0 * x[idx - 1] + 25.0 * x[idx] / 12.0) / dt;
      break;
    case 5:
      dxdt[idx] = (-x[idx - 5] / 5.0 + 1.2 * x[idx - 4] - 10.0*x[idx - 3] / 3.0 + 5.0 * x[idx - 2] - 5.0 * x[idx - 1] + 137.0 * x[idx] / 60.0) / dt;
      break;
    case 6:
      dxdt[idx] = (x[idx - 6] / 6.0 - 6.0 * x[idx - 5] / 5.0 + 15.0 * x[idx - 4] / 4.0 - 20.0 * x[idx - 3] / 3.0 + 7.5 * x[idx - 2] - 6.0 * x[idx - 1] + (49.0 / 20.0) * x[idx])/dt;
      break;
    default:
      break;
    }//end switch

    PetscFunctionReturn(0);
  }

  static PetscErrorCode ForwardFd(const PetscScalar *x, PetscScalar *dxdt, PetscInt idx, PetscReal dt, PetscInt order)
  {
    PetscFunctionBeginUser;
    if(order == 0 or order > 6){
      SETERRQ(PETSC_COMM_WORLD, 1, "Must supply forward finite difference order greater than 0 but less than 7.");
    }
    switch(order){
    case 1:
      dxdt[idx] = (x[idx + 1] - x[idx]) / dt;
      break;
    case 2:
      dxdt[idx] = (-0.5 * x[idx + 2] + 2.0 * x[idx + 1] - 1.5 * x[idx]) / dt;
      break;
    case 3:
      dxdt[idx] = (x[idx + 3] / 3.0 - 1.5 * x[idx + 2] + 3.0 * x[idx + 1] - 11.0 * x[idx] / 6.0) / dt;
      break;
    case 4:
      dxdt[idx] = (-x[idx + 4]/4.0 + 4 * x[idx + 3] / 3.0 - 3.0 * x[idx + 2] + 4.0 * x[idx + 1] - 25.0 * x[idx] / 12.0) / dt;
      break;
    case 5:
      dxdt[idx] = (x[idx + 5] / 5.0 - 1.2 * x[idx + 4] + 10.0*x[idx + 3] / 3.0 - 5.0 * x[idx + 2] + 5.0 * x[idx + 1] - 137.0 * x[idx] / 60.0) / dt;
      break;
    case 6:
      dxdt[idx] = (-x[idx + 6] / 6.0 + 6.0 * x[idx + 5] / 5.0 - 15.0 * x[idx + 4] / 4.0 + 20.0 * x[idx + 3] / 3.0 - 7.5 * x[idx + 2] + 6.0 * x[idx + 1] - (49.0 / 20.0) * x[idx])/dt;
      break;
    default:
      break;
    }//end switch

    PetscFunctionReturn(0);
  }

  static PetscErrorCode CentralFd(const PetscScalar *x, PetscScalar *dxdt, PetscInt idx, PetscReal dt, PetscInt order)
  {
    PetscFunctionBeginUser;
    if(order == 0 or order > 8 or order % 2 == 1){
      SETERRQ(PETSC_COMM_WORLD, 1, "Must supply an even central finite difference order greater than 0 but less than 9 (i.e. 2, 4, 6, or 8).");
    }
    switch(order){
    case 2:
      dxdt[idx] = (x[idx + 1] - x[idx - 1]) / (2 * dt);
      break;
    case 4:
      dxdt[idx] = (-x[idx + 2] / 12.0 + 2.0 * x[idx + 1] / 3.0 - 2.0 * x[idx - 1] / 3.0 + x[idx - 2] / 12.0) / dt;
      break;
    case 6:
      dxdt[idx] = (x[idx + 3] / 60.0 - 3.0 * x[idx + 2] / 20.0 + 0.75 * x[idx + 1] - 0.75 * x[idx - 1] + 3.0 * x[idx - 2] / 20.0 - x[idx - 3] / 60.0) / dt;
      break;
    case 8:
      dxdt[idx] = (-x[idx + 4] / 280.0 + 4.0 * x[idx + 3] / 105.0 - 0.2 * x[idx + 2] + 0.8 * x[idx + 1]
		   -0.8 * x[idx - 1] + 0.2 * x[idx - 2] - 4.0 * x[idx - 3] / 105.0 + x[idx - 4] / 280.0) / dt;
      break;
    default:
      break;
    }//end switch

    PetscFunctionReturn(0);
  }

  
  PetscErrorCode FiniteDifference(Vec X, Vec dXdt, PetscReal dt, PetscInt order,
				  StencilType direction)
  {
    PetscErrorCode ierr;
    PetscInt       i, nx, xs, xe, xps, xpe;
    const PetscScalar *x;
    PetscScalar       *dxdt;

    PetscFunctionBeginUser;
    ierr = VecGetOwnershipRange(X, &xs, &xe);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(dXdt, &xps, &xpe);CHKERRQ(ierr);
    if(xs != xps or xe != xpe){
      SETERRQ4(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "Must pass X and dXdt Vecs with the same local sizes! X has local ownwership range [%d, %d], but dXdt has range [%d, %d].\n", xs, xe, xps, xpe);
    }
    nx = xe - xs;
    if(nx < 2 * order){
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "Must pass array of (local) length at least twice the finite difference order!");
    }

    ierr = VecGetArrayRead(X, &x);CHKERRQ(ierr);
    ierr = VecGetArray(dXdt, &dxdt);CHKERRQ(ierr);

    for(i=xs; i < xe; ++i){
      switch(direction){
      case Forward:
	{
	  if(i < nx - order){
	    ierr = ForwardFd(x, dxdt, i - xs, dt, order);CHKERRQ(ierr);
	  } else {
	    /* near end of the array, use backward FD */
	    ierr = BackwardFd(x, dxdt, i - xs, dt, order);CHKERRQ(ierr);
	  }
	  break;
	}/*end case Forward*/
      case Backward:
	{
	  if(i - xs >= order){
	    ierr = BackwardFd(x, dxdt, i - xs, dt, order);CHKERRQ(ierr);
	  } else {
	    /* near start of array, use forward FD */
	    ierr = ForwardFd(x, dxdt, i - xs, dt, order);CHKERRQ(ierr);
	  }
	  break;
	}
      case Central:
	{
	  if(i - xs >= order/2 and i - xs < nx - order/2){
	    ierr = CentralFd(x, dxdt, i - xs, dt, order);CHKERRQ(ierr);
	  } else if(i < order / 2){
	    ierr = ForwardFd(x, dxdt, i - xs, dt, std::min(6, order));CHKERRQ(ierr);
	  } else {
	    ierr = BackwardFd(x, dxdt, i - xs, dt, std::min(6, order));CHKERRQ(ierr);
	  }
	  break;
	}
      }/* end switch */
    }
    ierr = VecRestoreArrayRead(X, &x);CHKERRQ(ierr);
    ierr = VecRestoreArray(dXdt, &dxdt);CHKERRQ(ierr);
      
    PetscFunctionReturn(ierr);
  }
        
	
}
	
      

    
