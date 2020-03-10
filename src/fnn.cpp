#include <nlts/fnn.h>
#include <tuple>
namespace nlts
{
  static PetscErrorCode ctor_compute_sizes(Vec X, PetscInt *nx, PetscInt *rank, PetscInt *size)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = VecGetLocalSize(X, nx);CHKERRQ(ierr);
    MPI_Comm_rank(PETSC_COMM_WORLD, rank);
    MPI_Comm_size(PETSC_COMM_WORLD, size);
    PetscFunctionReturn(ierr);
  }

  PetscErrorCode FNN::setFromOptions()
  {
    bool flg;
    PetscErrorCode ierr;
    std::string ifs, ofs;
    PetscFunctionBeginUser;
    std::tie(ifs, flg) = nlts::get_petsc_option<std::string>("--data-file");
    if(flg){
      infile = ifs;
      if(!X){
	ierr = VecReadBinary(infile, &X);CHKERRQ(ierr);
      }
    }
    std::tie(ofs, flg) = nlts::get_petsc_option<std::string>("-o");
    if(flg){
      outfile = ofs;
    }
    
    std::tie(verbose, flg) = nlts::get_petsc_option<PetscInt>("--verbose");
    if(!flg){
      std::tie(verbose, flg) = nlts::get_petsc_option<PetscInt>("-v");
      if(!flg){
	verbose=0;
      }
    }
    std::tie(ignore_first, flg) = nlts::get_petsc_option<PetscInt>("-x");
    if(!flg){
      ignore_first=0;
    }
    std::tie(maxdim, flg) = nlts::get_petsc_option<PetscInt>("--max-dim");
    if(!flg){
      std::tie(maxdim, flg) = nlts::get_petsc_option<PetscInt>("-dmax");
      if(!flg){
	maxdim=10;
      }
    }

    std::tie(mindim, flg) = nlts::get_petsc_option<PetscInt>("--min-dim");
    if(!flg){
      std::tie(mindim, flg) = nlts::get_petsc_option<PetscInt>("-dmin");
      if(!flg){
	mindim=1;
      }
    }

    std::tie(nbox, flg) = nlts::get_petsc_option<PetscInt>("-b");
    if(!flg){
      nbox = 1024;
    }

    std::tie(escape, flg) = nlts::get_petsc_option<PetscReal>("--escape");
    if(!flg){
      std::tie(escape, flg) = nlts::get_petsc_option<PetscReal>("-ef");
      if(!flg){
	escape = 10.0;
      }
    }

    std::tie(theiler, flg) = nlts::get_petsc_option<PetscInt>("--theiler-window");
    if(!flg){
      std::tie(theiler, flg) = nlts::get_petsc_option<PetscInt>("-tw");
      if(!flg){
	theiler = 0;
      }
    }

    PetscFunctionReturn(0);
  }
    
	     
	      

  FNN::FNN(Vec data, bool from_options)
  {
    X = data;
    ierr = ctor_compute_sizes(X, &nx, &rank, &size);
    if(from_options){
      ierr = setFromOptions();
    }
  }
      
  FNN::FNN(bool from_options)
  {
    ierr = setFromOptions();
    ierr = ctor_compute_sizes(X, &nx, &rank, &size);
  }
  
  FNN::FNN(std::string infilename)
  {
    infile = infilename;
    ierr = VecReadBinary(infile, &X);
    ierr = ctor_compute_sizes(X, &nx, &rank, &size);
  }

  FNN::FNN(std::string infilename, std::string outfilename)
  {
    infile = infilename;
    outfile = outfilename;
    ierr = VecReadBinary(infile, &X);
    ierr = ctor_compute_sizes(X, &nx, &rank, &size);
  }

  PetscErrorCode FNN::GatherProgress()
  {
    PetscFunctionBeginUser;
    MPI_Allreduce(&local_done_yet, &done_yet, 1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(&local_done, &done, 1, MPI_LONG, MPI_LAND, PETSC_COMM_WORLD);
    MPI_Allreduce(local_nearest.data(), nearest.data(), nx, MPI_LONG, MPI_LOR, PETSC_COMM_WORLD);
    PetscFunctionReturn(0);
  }

  /*PetscErrorCode FNN::nearestNeighbor(PetscInt idx, PetscInt dim, PetscReal eps, PetscInt *nearest_idx)
  {
    PetscInt    xv, yv, x1, x2, y, i, i1, ibox;
    long        elem, which=-1;
    PetscReal   dx, dx_min, dx_max = 1.0, fact;
    PetscScalar *x;
    PetscFunctionBeginUser;
  */
    /* NOTE: function largely copied from P2NTSA implementation,
       which is why the indexing rules are so weird */
  /*
    ibox = nbox - 1;
    ierr = VecGetArray(X, &x);
    xv = (PetscInt)(x[idx - (dim-1)*lag]/eps)&ibox;
    yv = (PetscInt)(x[idx]/eps)&ibox;
    */


}
