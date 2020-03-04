#include <nlts/io.h>

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Vec            x, t, tt;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);

  auto filename = std::string{"data/ps8data/data2"};

  auto [xv, tv] = nlts::read_scalar_trajectory(filename);
  std::cout << "Data has " << std::scientific << double(xv.size()) << " rows.\n";
  ierr = nlts::VecReadScalarTrajectory(filename, &x, &t);CHKERRQ(ierr);
  ierr = nlts::VecWriteScalarTrajectory("data", x, t);CHKERRQ(ierr);
  ierr = VecCreateShared(PETSC_COMM_WORLD, PETSC_DECIDE, xv.size(), &tt);
  //ierr = nlts::VecReadBinary("trajectory_data/data_t.dat", tt);CHKERRQ(ierr);

  ierr = VecView(tt, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);

  return 0;
}
