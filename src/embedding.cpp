#include <nlts/embedding.h>
#include <cmath>
#include <algorithm>
namespace nlts
{
  PetscInt TakensEmbedding::DataLength()
  {
    PetscInt N;
    ierr = VecGetSize(xdata, &N);CHKERRQ(ierr);
    return N;
  }
  
  
  PetscErrorCode TakensEmbedding::AllocateEmbeddingVecs()
  {
    int evlen, i = 0;
    PetscFunctionBeginUser;
    for(auto& evec : embedding){
      ierr = VecCreate(PETSC_COMM_WORLD, &evec);CHKERRQ(ierr);
      ierr = VecSetFromOptions(evec);CHKERRQ(ierr);
      evlen = n - i*std::ceil(tau/dt);
      ierr = VecSetSizes(evec, PETSC_DECIDE, evlen);CHKERRQ(ierr);
      ierr = VecGetLocalSize(evec, &evlen);
      evlens.push_back(evlen);
      ++i;
    }
    PetscFunctionReturn(ierr);
  }

  PetscErrorCode TakensEmbedding::FreeEmbeddingVecs()
  {
    PetscFunctionBeginUser;
    for(auto& evec : embedding){
      ierr = VecDestroy(&evec);CHKERRQ(ierr);
    }
    PetscFunctionReturn(ierr);
  }
    
  
  TakensEmbedding::TakensEmbedding(Vec state_data, PetscReal start_time,
		    PetscReal end_time, PetscReal sample_dt,
				   PetscReal Tau, PetscInt M,
				   std::optional<PetscReal> data_period)
		    : tau{Tau}, t0{start_time}, tf{end_time}, dt{sample_dt},
		      m{M}, xdata{state_data}
    {
      if(data_period){
	T = *data_period;
      } else {
	T = 0.0;
      }
      
      n = DataLength();
      embedding.resize(n);
      AllocateEmbeddingVecs();
    };

  TakensEmbedding::~TakensEmbedding()
  {
    FreeEmbeddingVecs();
  }

  PetscErrorCode TakensEmbedding::mapEmbedding(const std::function<PetscReal(PetscReal)>& f)
  {
    PetscScalar *x;
    PetscFunctionBeginUser;
    for(int i = 0; i < m; ++i){
      ierr = VecGetArray(embedding[i], &x);CHKERRQ(ierr);
      ierr = VecGetLocalSize(embedding[i], &evlens[i]);CHKERRQ(ierr);
      for(int j = 0; j < evlens[i]; ++j){
	x[i] = f(x[i]);
      }
      ierr = VecRestoreArray(embedding[i], &x);CHKERRQ(ierr);
    }
    PetscFunctionReturn(ierr);
  }

  std::pair<int, bool> TakensEmbedding::IndexFromTime(PetscReal t)
  {
    int id = std::round((t - t0) / dt);
    bool exact_time = (std::abs((PetscReal)(id * dt) + t0 - t) < 1.0e-10);
    return {id, exact_time};
  }

  
  PetscErrorCode TakensEmbedding::embed()
  {
    int i, j;
    PetscInt  nx, rank, xlow, xhigh;
    PetscReal t;
    const PetscScalar *x;
    PetscScalar *ex;
    PetscFunctionBeginUser;
    ierr = VecGetArrayRead(xdata, &x);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(xdata, &xlow, &xhigh);CHKERRQ(ierr);
    nx = xhigh - xlow;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for(i=0; i < m; ++i){
      ierr = VecGetLocalSize(embedding[i], &evlens[i]);CHKERRQ(ierr);
      if(evlens[i] > nx){
	SETERRQ2(PETSC_COMM_WORLD, 1, "Error, embedding Vec has larger local ownership range (%d elements) than data Vec (%d elements).", evlens[i], nx);
      }
      ierr = VecGetArray(embedding[i], &ex);CHKERRQ(ierr);
      t = t0 + (PetscReal)xlow * dt;
      for(j=0; j < evlens[i]; ++j){
	auto [id, is_exact_time] = IndexFromTime(t + (PetscReal)i * tau);
	if(is_exact_time and id < n){
	  ex[j] = x[id];
	} else if(id >= nx){
	  SETERRQ(PETSC_COMM_WORLD, 1, "Error, somehow index id in TakensEmbedding::embed is out of range.");
	} else {
	  PetscPrintf(PETSC_COMM_WORLD, "WARNING: interpolation currently not implemented. ignoring this index.");CHKERRQ(ierr);
	}
	t += dt;
      }
      ierr = VecRestoreArray(embedding[i], &ex);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(xdata, &x);CHKERRQ(ierr);
    PetscFunctionReturn(ierr);
  }

  PetscErrorCode TakensEmbedding::plotEmbedding(PetscInt x_eid, PetscInt y_eid,
						std::optional<std::string> title,
						std::optional<std::string> xlabel,
						std::optional<std::string> ylabel,
						std::optional<int> xpixels,
						std::optional<int> ypixels,
						std::optional<int> hold_setting)
  {
    VecScatter           ctx;
    PetscInt             i, nx, ny, rank;
    const PetscScalar    *x;
    PetscScalar          *xcpy;
    PetscFunctionBeginUser;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* get relevant data onto rank 0 */
    ierr = VecScatterCreateToZero(embedding[x_eid], &ctx, &plotx);
    ierr = VecScatterBegin(ctx, embedding[x_eid], plotx, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, embedding[x_eid], plotx, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

    ierr = VecScatterCreateToZero(embedding[y_eid], &ctx, &ploty);
    ierr = VecScatterBegin(ctx, embedding[y_eid], ploty, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, embedding[y_eid], ploty, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
    
    ierr = VecGetSize(plotx, &nx);CHKERRQ(ierr);
    ierr = VecGetSize(ploty, &ny);CHKERRQ(ierr);
    if(nx > ny and !rank){
      /* shorten plotx */
      ierr = VecCreateSeq(PETSC_COMM_WORLD, ny, &plot_buff);CHKERRQ(ierr);
      ierr = VecGetArrayRead(plotx, &x);
      ierr = VecGetArray(plot_buff, &xcpy);
      for(i=0; i < ny; ++i){
	xcpy[i] = x[i];
      }
      ierr = VecRestoreArray(plot_buff, &xcpy);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(plotx, &x);CHKERRQ(ierr);

      ierr = VecDestroy(&plotx);CHKERRQ(ierr);
      ierr = VecDuplicate(plot_buff, &plotx);CHKERRQ(ierr);
      ierr = VecCopy(plot_buff, plotx);CHKERRQ(ierr);
      ierr = VecDestroy(&plot_buff);CHKERRQ(ierr);
    } else if(ny > nx and !rank){
      /* shorten ploty */
      ierr = VecCreateSeq(PETSC_COMM_WORLD, nx, &plot_buff);CHKERRQ(ierr);
      ierr = VecGetArrayRead(ploty, &x);
      ierr = VecGetArray(plot_buff, &xcpy);
      for(i=0; i < nx; ++i){
	xcpy[i] = x[i];
      }
      ierr = VecRestoreArray(plot_buff, &xcpy);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(ploty, &x);CHKERRQ(ierr);

      ierr = VecDestroy(&ploty);CHKERRQ(ierr);
      ierr = VecDuplicate(plot_buff, &ploty);CHKERRQ(ierr);
      ierr = VecCopy(plot_buff, ploty);CHKERRQ(ierr);
      ierr = VecDestroy(&plot_buff);CHKERRQ(ierr);
    }

    ierr = VecGetSize(plotx, &nx);CHKERRQ(ierr);
    ierr = VecGetSize(ploty, &ny);CHKERRQ(ierr);

    if(nx != ny){
      SETERRQ2(PETSC_COMM_WORLD, 1, "Error, after resizing, plotx has %d elements but ploty has %d elements.\n", nx, ny);
    }

    ierr = PlotVecs(plotx, ploty, title, xlabel, ylabel,
		    xpixels, ypixels, hold_setting);CHKERRQ(ierr);

    PetscFunctionReturn(ierr);
  }

  /*
  std::vector<Vec> TakensEmbedding::embeddedTrajectory(bool truncate)
  {
    if(!truncate){
      return embedding;
    } else {
      std::vector<Vec> trunc_embed(embedding.size());
      int max_len = *std::min_element(evlens.begin(), evlens.end());
      for(auto i = 0; i < embedding.size(); ++i){
	ierr = VecCreate(PETSC_COMM_WORLD, &trunc_embed[i]);
	ierr = VecSetFromOptions(trunc_embed[i]);
	ierr = VecSetSizes(trunc_embed[i], PETSC_DECIDE, max_len);
	PetscScalar *tx;
	const PetscScalar *x;
	ierr = VecGetArrayRead(embedding[i], &x);
	ierr = VecGetArray(trunc_embed[i], &tx);
	for(auto j = 0; j < max_len; ++j){
	  tx[j] = x[j];
	}
	ierr = VecRestoreArray(trunc_embed[i], &tx);
	ierr = VecRestoreArrayRead(embedding[i], &x);
      }
      return trunc_embed;
    }
  }
  */

  PetscErrorCode TakensEmbedding::status() const noexcept
  {
    return ierr;
  }

  PetscReal TakensEmbedding::embeddingInterval() const noexcept
  {
    return tau;
  }

  PetscReal TakensEmbedding::timeWindowLength() const noexcept
  {
    return tf;
  }

  PetscReal TakensEmbedding::initialTime() const noexcept
  {
    return t0;
  }

  PetscReal TakensEmbedding::samplingInterval() const noexcept
  {
    return dt;
  }
  
	

}
