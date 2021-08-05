static char help[] = "NLTS implementation of TISEAN's mutual program.\n\
Options:\n\
     [-f --file] : Path/name of binary file containing the data as a PETSc Vec.\n\
     [-x --ignore-first] : Ignore the first this-many lines of the data file (default 0).\n\
     [-b -nbox --num-box] : Number of boxes for the partition. (default 16).\n\
     [-D --max-tau --max-lag] : Maximum lag (in number of timesteps, defaults to 200).\n\
     [-o] : output file (default stdout)\n\
     [-p --plot] : plot mutual information versus lag? (default false)\n\
     [-n --no-entropy-vs-lag] : Do not print the entropy versus lag to stdout, even if the output file is not specified (only print the minimal entropy and corresponding lag).\n";

#include <nlts/mutual.h>
#include <nlts/io.h>
#include <nlts/plot.h>
#include <nlts/sys.h>
#include <algorithm>

int main(int argc, char **argv)
{
  PetscErrorCode    ierr;
  Vec               X, Y, Tau, H;
  PetscScalar       *vals;
  IS                is;
  PetscInt          taumax, nbox, ignore_first, nx, i;
  std::string       infile, outfile;
  bool              flg, outfile_flg, plot_flg, no_entropy;
  std::optional<
    std::vector<PetscReal>
    >               entropy;
  PetscMPIInt       size;
 

  ierr = nlts::Initialize(&argc, &argv, NULL, help);CHKERRQ(ierr);
  if(nlts::has_petsc_option("-h")) return 0;
  
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if(size > 1){
    SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "Error: mutual can currently only run with one process, not %d!\n", size);
  }

  std::tie(infile, flg) = nlts::get_petsc_option<std::string>("-f");
  if(!flg){
    std::tie(infile, flg) = nlts::get_petsc_option<std::string>("--file");
    if(!flg){
      PetscPrintf(PETSC_COMM_WORLD, "WARNING: input file not found in options (use flags -f or --file). Trying to read argv[1]=%s as an input file. This may crash.\n", argv[1]);
      infile = std::string{argv[1]};
    }
  }

  std::tie(ignore_first, flg) = nlts::get_petsc_option<PetscInt>("-x");
  if(!flg){
    std::tie(ignore_first, flg) = nlts::get_petsc_option<PetscInt>("--ignore-first");
    if(!flg){
      ignore_first = 0;
    }
  }

  std::tie(nbox, flg) = nlts::get_petsc_option<PetscInt>("-b");
  if(!flg){
    std::tie(nbox, flg) = nlts::get_petsc_option<PetscInt>("-nbox");
    if(!flg){
      std::tie(nbox, flg) = nlts::get_petsc_option<PetscInt>("--num-box");
      if(!flg){
	nbox = 16;
      }
    }
  }

  std::tie(taumax, flg) = nlts::get_petsc_option<PetscInt>("-D");
  if(!flg){
    std::tie(taumax, flg) = nlts::get_petsc_option<PetscInt>("--max-tau");
    if(!flg){
      std::tie(taumax, flg) = nlts::get_petsc_option<PetscInt>("--max-lag");
      if(!flg){
	taumax = 200;
      }
    }
  }

  ierr = VecCreateSeq(PETSC_COMM_SELF, taumax, &Tau);CHKERRQ(ierr);
  ierr = VecGetArray(Tau, &vals);CHKERRQ(ierr);
  for(i=0; i<taumax; ++i){
    vals[i] = i;
  }
  ierr = VecRestoreArray(Tau, &vals);CHKERRQ(ierr);
  
  std::tie(outfile, outfile_flg) = nlts::get_petsc_option<std::string>("-o");
  no_entropy = nlts::has_petsc_option("-n");
  if (!no_entropy) {
    no_entropy = nlts::has_petsc_option("--no-entropy-vs-lag");
  }

  plot_flg = nlts::has_petsc_option("-p");
  if(!plot_flg){
    plot_flg = nlts::has_petsc_option("--plot");
  }

  ierr = nlts::VecReadBinary(infile, &X);CHKERRQ(ierr);
  if(ignore_first > 0){
    ierr = VecGetSize(X, &nx);CHKERRQ(ierr);
    ierr = ISCreateStride(PETSC_COMM_WORLD, nx - ignore_first, ignore_first, 1, &is);CHKERRQ(ierr);
    ierr = VecGetSubVector(X, is, &Y);
    ierr = ISDestroy(&is);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Minimizing mutual information.\n");
    auto [hmtau, oerr] = nlts::MinimizeMutualInformation(Y, taumax,
							 2, nbox, 0.01);
    PetscPrintf(PETSC_COMM_WORLD, "Tau that minimizes mutual information: %d.\n", hmtau);
    std::tie(entropy, ierr) = nlts::MutualInformation(Y, std::max(taumax, hmtau),
						      nbox,
						      (outfile_flg ? std::optional(outfile)
						       : std::nullopt),
						      (!outfile_flg or plot_flg));CHKERRQ(ierr);
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Minimizing mutual information.\n");
    auto [hmtau, oerr] = nlts::MinimizeMutualInformation(X, taumax,
							 2, nbox, 0.01);
    PetscPrintf(PETSC_COMM_WORLD, "Tau that minimizes mutual information: %d.\n", hmtau);
    std::tie(entropy, ierr) = nlts::MutualInformation(X, taumax, nbox,
						      (outfile_flg ? std::optional(outfile)
						       : std::nullopt),
						      (!outfile_flg or plot_flg));CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "entropy has %D elements\n", entropy->size());
    PetscPrintf(PETSC_COMM_WORLD, "Minimal mutual information is %g nats.\n",*std::min_element(entropy->begin(), entropy->end()));
  }

  if(plot_flg and entropy){
    ierr = nlts::VecSetFromStd(&H, *entropy);CHKERRQ(ierr);
    ierr = nlts::PlotVecs(Tau, H, "Mutual Information versus Lag",
			  "Lag (tau)", "Mutual Information (Shannon Entropy)",
			  2000,2000, std::nullopt, "mutual_info.png");CHKERRQ(ierr);
    ierr = VecDestroy(&H);
  } else if (entropy and !outfile_flg and !no_entropy) {
    /* if there is an output file, no need to print to stdout */
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Entropy vs lag:\n");CHKERRQ(ierr);
    PetscInt i = 0;
    for (const auto& h : *entropy) {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Lag %D: entropy = %g nats\n", i, h);CHKERRQ(ierr);
      ++i;
    }
  }

  VecDestroy(&Tau);
  VecDestroy(&X);
  ierr = nlts::Finalize();
  return ierr;
}

  
  
  
    
  
