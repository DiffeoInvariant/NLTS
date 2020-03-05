#include <nlts/plot.h>
#include <algorithm>
#include <numeric>

namespace nlts
{

  PetscErrorCode PlotVecs(Vec X, Vec Y,
			  std::optional<std::string> title,
			  std::optional<std::string> xlabel,
			  std::optional<std::string> ylabel,
			  std::optional<int> xpixels,
			  std::optional<int> ypixels,
			  std::optional<int> hold_setting,
			  std::optional<std::string> save_file,
			  std::optional<std::pair<PetscReal, PetscReal>> xlims,
			  std::optional<std::pair<PetscReal, PetscReal>> ylims)
  {
    PetscErrorCode ierr;
    PetscDraw      draw;
    PetscDrawSP    scatter;
    PetscDrawAxis  ax;
    const char     *xlab, *ylab, *ttl;
    int            width, height, hold;
    PetscReal      *x, *y;
    PetscInt       i, nx, ny, n;

    PetscFunctionBeginUser;

    if(title){
      ttl = title->c_str();
    } else {
      ttl = NULL;
    }
    if(xlabel){
      xlab = xlabel->c_str();
    } else {
      xlab = NULL;
    }
    if(ylabel){
      ylab = ylabel->c_str();
    } else {
      ylab = NULL;
    }
    if(xpixels){
      width = *xpixels;
    } else {
      width = 500;
    }
    if(ypixels){
      height = *ypixels;
    } else {
      height = 500;
    }
    if(hold_setting){
      hold = *hold_setting;
    } else {
      hold = -2;
    }
      
    ierr = VecGetLocalSize(X, &nx);CHKERRQ(ierr);
    ierr = VecGetLocalSize(Y, &ny);CHKERRQ(ierr);
    n = std::min(nx, ny);

    ierr = PetscDrawCreate(PETSC_COMM_SELF, NULL, ttl, 0, 0, width, height, &draw);CHKERRQ(ierr);
    if(save_file){
      ierr = PetscDrawSetSaveFinalImage(draw, save_file->c_str());CHKERRQ(ierr);
    }
    ierr = PetscDrawPointSetSize(draw, 0.0);CHKERRQ(ierr);
    ierr = PetscDrawSetFromOptions(draw);CHKERRQ(ierr);
    ierr = PetscDrawSPCreate(draw, 1, &scatter);CHKERRQ(ierr);
    if(xlab or ylab){
      ierr = PetscDrawSPGetAxis(scatter, &ax);CHKERRQ(ierr);
      ierr = PetscDrawAxisSetLabels(ax, NULL, xlab, ylab);CHKERRQ(ierr);
    }

    ierr = VecGetArray(X, &x);CHKERRQ(ierr);
    ierr = VecGetArray(Y, &y);CHKERRQ(ierr);
    for(i = 0; i < n; ++i){
      ierr = PetscDrawSPAddPoint(scatter, &x[i], &y[i]);CHKERRQ(ierr);
      /*ierr = PetscDrawSPDraw(scatter, PETSC_FALSE);CHKERRQ(ierr);*/
    }
    ierr = PetscDrawSPDraw(scatter, PETSC_FALSE);CHKERRQ(ierr);
    if(xlims and ylims){
      ierr = PetscDrawSPSetLimits(scatter, xlims->first, xlims->second,
				  ylims->first, ylims->second);CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(X, &x);CHKERRQ(ierr);
    ierr = VecRestoreArray(Y, &y);CHKERRQ(ierr);
    if(save_file){
      ierr = PetscDrawSave(draw);CHKERRQ(ierr);
    }
    ierr = PetscDrawSetPause(draw, hold);CHKERRQ(ierr);
    ierr = PetscDrawSPDestroy(&scatter);CHKERRQ(ierr);
    ierr = PetscDrawDestroy(&draw);CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
  }


}
