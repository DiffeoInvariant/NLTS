#ifndef NLTS_PLOT_H
#define NLTS_PLOT_H
#include <petscdraw.h>
#include <petscvec.h>
#include <optional>
#include <string>
#include <utility>


namespace nlts
{

  PetscErrorCode PlotVecs(Vec X, Vec Y,
			  std::optional<std::string> title=std::nullopt,
			  std::optional<std::string> xlabel=std::nullopt,
			  std::optional<std::string> ylabel=std::nullopt,
			  std::optional<int> xpixels=std::nullopt,
			  std::optional<int> ypixels=std::nullopt,
			  std::optional<int> hold_setting=std::nullopt,
			  std::optional<std::pair<PetscReal, PetscReal>> xlims=std::nullopt,
			  std::optional<std::pair<PetscReal, PetscReal>> ylims=std::nullopt);

}
#endif
