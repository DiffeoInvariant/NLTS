#include <nlts/io.h>
#include <list>

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Vec            x, t;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  std::list<std::string> filenames = {std::string{"data/ps8data/data1"},
				      std::string{"data/ps8data/data2"},
				      std::string{"data/ps8data/data3"},
                                      std::string{"data/ps8data/data2.first250sec"},};

  for(const auto& fn : filenames){
    ierr = nlts::VecReadScalarTrajectory(fn, &x, &t);CHKERRQ(ierr);
    ierr = nlts::VecWriteScalarTrajectory("data", x, t);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&t);CHKERRQ(ierr);
  }

  ierr = PetscFinalize();

  return ierr;
}
  
						  
