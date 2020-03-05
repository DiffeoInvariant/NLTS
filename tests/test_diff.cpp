#include <nlts/diff.h>
#include <nlts/io.h>
#include <nlts/plot.h>
#include <petscmath.h>

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Vec            X, T, dXdT;
  PetscScalar    *x, xv;
  PetscReal      dt;
  PetscInt       order,i,j,nx;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = nlts::VecReadBinary("trajectory_data/data_t.dat", &T);CHKERRQ(ierr);
  ierr = nlts::VecReadBinary("trajectory_data/data_x.dat", &X);CHKERRQ(ierr);

  ierr = VecGetLocalSize(X, &nx);CHKERRQ(ierr);
  ierr = VecGetArray(X, &x);CHKERRQ(ierr);
  for(i=0; i<nx; ++i){
    xv = PetscFmodReal(x[i], 2 * M_PI);
    j=0;
    while(xv > 5.5){
      xv = PetscFmodReal(x[i] -  M_PI + j * 1.0e-8, 2 * M_PI);
      ++j;
      }
    x[i] = xv;
  }
  ierr = VecRestoreArray(X, &x);CHKERRQ(ierr);
  ierr = VecDuplicate(X, &dXdT);CHKERRQ(ierr);
  dt = 0.001;
  order = 6;
  ierr = nlts::FiniteDifference(X, dXdT, dt, order, nlts::Central);CHKERRQ(ierr);

  ierr = VecGetArray(dXdT, &x);CHKERRQ(ierr);
  for(i=0; i<nx; ++i){
    if(PetscAbsReal(x[i]) > 2.0){
      x[i] = 0.0;
    }
  }
  ierr = VecRestoreArray(dXdT, &x);CHKERRQ(ierr);
  
  ierr = nlts::PlotVecs(T, X, "x vs t", "t",
			"x", 2500, 2500, std::nullopt, "diff_plot.png",
			std::make_pair(0.0, 45.0), std::make_pair(2.9,3.2));CHKERRQ(ierr);

  ierr = nlts::PlotVecs(T, dXdT, "dx/dt vs t", "t",
			"dx/dt", 2500, 2500);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "X:\n");
  ierr = VecView(X, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "\n\n\n\n dXdt:\n");
  ierr = VecView(dXdT, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&T);CHKERRQ(ierr);
  ierr = VecDestroy(&dXdT);CHKERRQ(ierr);
  PetscFinalize();
  return 0;
}
