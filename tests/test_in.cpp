#include <nlts/io.h>

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Vec            t, x;

   ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);

   ierr = nlts::VecReadBinary("trajectory_data/data_t.dat", &t);
   ierr = VecView(t,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   ierr = VecDestroy(&t);CHKERRQ(ierr);

   return 0;
}
