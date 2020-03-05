#include <nlts/sys.h>

namespace nlts
{
  PetscErrorCode VecSetFromStd(Vec *v, const std::vector<PetscScalar>& data)
  {
    PetscErrorCode ierr;
    PetscScalar    *x;
    std::size_t    i;
    PetscFunctionBeginUser;
    ierr = VecCreate(PETSC_COMM_WORLD, v);CHKERRQ(ierr);
    ierr = VecSetFromOptions(*v);CHKERRQ(ierr);
    ierr = VecSetSizes(*v, data.size(), PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecGetArray(*v, &x);CHKERRQ(ierr);
    for(i=0; i<data.size(); ++i){
      x[i] = data[i];
    }
    ierr = VecRestoreArray(*v, &x);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(*v);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*v);CHKERRQ(ierr);
    PetscFunctionReturn(ierr);
  }


}
