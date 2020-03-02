#include <nlts/io.h>

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Vec            x, t;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);

  auto filename = std::string{"data/ps8data/data2.first250sec"};

  ierr = nlts::VecReadScalarTrajectory(filename, &x, &t);CHKERRQ(ierr);
  ierr = nlts::VecWriteScalarTrajectory("data", x, t);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);

  return 0;
}
