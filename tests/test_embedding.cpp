#include <nlts/plot.h>
#include <nlts/io.h>
#include <nlts/embedding.h>

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Vec X, T;
  const PetscScalar *t;
  PetscScalar       *x, xv;
  PetscReal t0, tf, dt, tau;
  PetscInt m, nx, i, j, k;
  PetscBool flg;
  std::string xfile, tfile, data_dir, has_data_dir;

  ierr = nlts::Initialize(&argc, &argv);CHKERRQ(ierr);

  if(nlts::has_petsc_option("-data_dir")){
    std::tie(data_dir, has_data_dir) = nlts::get_petsc_option<std::string>("-data_dir");
  } else {
    data_dir = std::string{"trajectory_data"};
  }
  xfile = data_dir + std::string{"/data_x.dat"};
  tfile = data_dir + std::string{"/data_t.dat"};

  ierr = nlts::VecReadBinary(tfile, &T);CHKERRQ(ierr);
  ierr = nlts::VecReadBinary(xfile, &X);CHKERRQ(ierr);
  ierr = VecGetArrayRead(T, &t);CHKERRQ(ierr);
  t0 = t[0];
  ierr = VecRestoreArrayRead(T, &t);CHKERRQ(ierr);
  
  ierr = VecMax(T, NULL, &tf);CHKERRQ(ierr);
  
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
  ierr = PetscOptionsGetReal(NULL, NULL, "-tau", &tau, &flg);CHKERRQ(ierr);
  if(!flg){
    tau = 0.15;
  }
  ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &dt, &flg);CHKERRQ(ierr);
  if(!flg){
    dt = 0.001;
  }
  ierr = PetscOptionsGetInt(NULL, NULL, "-m", &m, &flg);CHKERRQ(ierr);
  if(!flg){
    m = 7;
  }
  ierr = PetscOptionsGetInt(NULL, NULL, "-j", &j, &flg);CHKERRQ(ierr);
  if(!flg){
    j = 0;
  }
  ierr = PetscOptionsGetInt(NULL, NULL, "-k", &k, &flg);CHKERRQ(ierr);
  if(!flg){
    k = 2;
  }

  {
  nlts::TakensEmbedding embedder(X, t0, tf, dt, tau, m);
  ierr = embedder.embed();CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "done embedding\n");
  std::string title{"Takens embedding with tau = "};
  title += std::to_string(tau);
  title += std::string{", m = "};
  title += std::to_string(m);
  title += std::string{", j = "};
  title += std::to_string(j);
  title += std::string{", k = "};
  title += std::to_string(k);
  title += std::string{"."};
  embedder.plotEmbedding(j, k, title, "j", "k", 2500, 2500);
  }
  PetscFinalize();
  return 0;
}
  
  
