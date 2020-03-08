#ifndef NLTS_MUTUAL_H
#define NLTS_MUTUAL_H
#include <petscvec.h>
#include <petscmath.h>
#include <optional>
#include <utility>
#include <tuple>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
namespace nlts
{
  extern std::pair<std::optional<std::vector<PetscReal>>,PetscErrorCode>
  MutualInformation(Vec X, std::optional<PetscInt> max_tau=std::nullopt,
		    std::optional<PetscInt> partition_boxes=std::nullopt,
		    std::optional<std::string> outfile=std::nullopt,
		    bool return_info=true);

  extern std::pair<PetscInt, PetscErrorCode> MinimizeMutualInformation(Vec X,
								       std::optional<PetscInt> max_tau=std::nullopt,
								       std::optional<PetscInt> tau_grad_stride=std::nullopt,
								       std::optional<PetscInt> partition_boxes=std::nullopt,
								       std::optional<PetscReal> tau_grad_min=std::nullopt);
				   
  


}
#endif
