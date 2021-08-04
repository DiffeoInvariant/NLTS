#include <nlts/diff.h>
#include <nlts/io.h>
#include <nlts/plot.h>
#include <nlts/sample.h>
#include <petscmath.h>
#include <utility>
#include <tuple>
#include <iostream>
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  PetscInt       rank, size;
  Vec            X, dsX, T, dsT, dXdT;
  PetscScalar    *x, xv;
  PetscReal      dt;
  PetscInt       order,i,j,nx, downsample;
  bool           do_ds;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  ierr = nlts::VecReadBinary("trajectory_data/data_t.dat", &T);CHKERRQ(ierr);
  ierr = nlts::VecReadBinary("trajectory_data/data_x.dat", &X);CHKERRQ(ierr);

  std::tie(downsample, do_ds) = nlts::get_petsc_option<PetscInt>("-dr");
  if(!do_ds){
    std::tie(downsample, do_ds) = nlts::get_petsc_option<PetscInt>("--downsample-rate");
  }
  ierr = VecGetLocalSize(X, &nx);CHKERRQ(ierr);
  std::cout << "Process " << rank << " of " << size << " owns " << nx << " elements of X.\n";
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
  if(do_ds){
    PetscPrintf(PETSC_COMM_WORLD, "Using decimating downsampling rate %d.\n", downsample);
    ierr = nlts::DecimatingDownsample(X, &dsX, downsample);CHKERRQ(ierr);
    ierr = nlts::DecimatingDownsample(T, &dsT, downsample);CHKERRQ(ierr);
    ierr = VecDuplicate(dsX, &dXdT);CHKERRQ(ierr);
    dt = downsample * 0.002;
    order = 6;
    ierr = nlts::FiniteDifference(dsX, dXdT, dt, order, nlts::Central);CHKERRQ(ierr);
  } else {
    ierr = VecDuplicate(X, &dXdT);CHKERRQ(ierr);
    dt = 0.002;
    order = 6;
    ierr = nlts::FiniteDifference(X, dXdT, dt, order, nlts::Central);CHKERRQ(ierr);
  }

  /*
  ierr = VecGetArray(dXdT, &x);CHKERRQ(ierr);
  ierr = VecGetLocalSize(dXdT, &nx);CHKERRQ(ierr);
  for(i=0; i<nx; ++i){
    if(PetscAbsReal(x[i]) > 3.5 or PetscAbsReal(x[i]) < 1.0){
      x[i] = 0.0;
    }
  }
  ierr = VecRestoreArray(dXdT, &x);CHKERRQ(ierr);
  */
  if(do_ds){
    ierr = nlts::PlotVecs(dsT, dsX, "x vs t", "t",
			  "x", 2500, 2500, std::nullopt, "diff_plot.png",
			  std::make_pair(0.0, 45.0/(PetscReal)downsample),
			  std::make_pair(1.0, 3.5));CHKERRQ(ierr);
    ierr = nlts::PlotVecs(T, dXdT, "dx/dt vs t", "t",
			  "dx/dt", 2500, 2500);CHKERRQ(ierr);
    ierr = nlts::PlotVecs(dsX, dXdT, "theta versus omega", "theta",
			"omega", 2500, 2500);CHKERRQ(ierr);
  }
  else {
    ierr = nlts::PlotVecs(T, X, "x vs t", "t",
			  "x", 2500, 2500, std::nullopt, "diff_plot.png",
			  std::make_pair(0.0, 45.0), std::make_pair(2.9,3.3));CHKERRQ(ierr);
    ierr = nlts::PlotVecs(T, dXdT, "dx/dt vs t", "t",
			  "dx/dt", 2500, 2500);CHKERRQ(ierr);
  }

  

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
