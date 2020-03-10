#ifndef NLTS_FNN_H
#define NLTS_FNN_H
#include <petscvec.h>
#include <nlts/sys.h>
#include <nlts/io.h>
#include <string>
#include <optional>
#include <cstddef>
#include <petscsys.h>
#include <cmath>

namespace nlts
{

  struct FNN
  {
    FNN(std::string infilename);

    FNN(std::string infilename, std::string outfilename);

    FNN(bool from_options=true);

    FNN(Vec data, bool from_options=true);

    PetscErrorCode setFromOptions();

    PetscErrorCode nearestNeighbor(PetscInt idx, PetscInt dim, PetscReal eps, PetscInt *nearest_idx);
  private:

    PetscErrorCode GatherProgress();
    Vec                                  X=NULL;
    PetscErrorCode                       ierr;
    std::string                          infile, outfile;
    long                                 done, local_done, partial_nbox, nbox=1024;
    
    std::vector<std::vector<long>>       box;
    
    std::vector<long>                    list, nearest, local_nearest;
    
    long                                 filelen, maxdim, mindim,
                                         lag, col, done_yet, local_done_yet,
                                         too_far, local_too_far,
	                                 verbose=0, theiler=0, ignore_first=0;
    
    PetscInt                             nx, rank, size, dim;
    
    PetscReal                            var, escape, aveps, local_aveps, vareps, local_vareps, *x;
  };



}
#endif
