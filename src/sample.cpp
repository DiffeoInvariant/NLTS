#include <nlts/sample.h>

namespace nlts
{
  PetscErrorCode DecimatingDownsample(Vec X, Vec *Y, PetscInt M)
  {
    PetscInt nx, ny, i;
    PetscErrorCode ierr;
    PetscScalar    *y;
    const PetscScalar *x;
    PetscFunctionBeginUser;
    ierr = VecGetLocalSize(X, &nx);CHKERRQ(ierr);
    ny = nx / M;
    ierr = VecCreate(PETSC_COMM_WORLD, Y);CHKERRQ(ierr);
    ierr = VecSetFromOptions(*Y);CHKERRQ(ierr);
    ierr = VecSetSizes(*Y, ny, PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecGetArrayRead(X, &x);CHKERRQ(ierr);
    ierr = VecGetArray(*Y, &y);CHKERRQ(ierr);
    for(i=0; i<ny; ++i){
      y[i] = x[i * M];
    }
    ierr = VecRestoreArray(*Y, &y);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(X, &x);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(*Y);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*Y);CHKERRQ(ierr);
    PetscFunctionReturn(ierr);
  }
    
    


}
