#ifndef NLTS_EMBEDDING_H
#define NLTS_EMBEDDING_H
#include <petscvec.h>
#include <vector>
#include <optional>
#include <string>
#include <functional>
#include <utility>
#include <nlts/plot.h>

namespace nlts
{

  struct TakensEmbedding
  {
    TakensEmbedding(Vec state_data, PetscReal start_time,
		    PetscReal end_time, PetscReal sample_dt,
		    PetscReal Tau, PetscInt M, std::optional<PetscReal> data_period=std::nullopt);

    PetscErrorCode setEmbeddingInterval(PetscReal);
    
    PetscErrorCode setSamplingInterval(PetscReal);
    
    PetscErrorCode setTimeWindow(PetscReal, PetscReal);
    
    PetscErrorCode setTimeWindow(PetscReal);
    
    PetscErrorCode setTrajectoryData(Vec);

    PetscErrorCode readTrajectoryData(std::string);

    PetscErrorCode writeEmbedding(std::optional<std::string>) const;

    PetscErrorCode mapEmbedding(const std::function<PetscReal(PetscReal)>&);
    
    PetscErrorCode embed();

    PetscErrorCode plotEmbedding(PetscInt x_eid, PetscInt y_eid,
				 std::optional<std::string> title=std::nullopt,
				 std::optional<std::string> xlabel=std::nullopt,
				 std::optional<std::string> ylabel=std::nullopt,
				 std::optional<int> xpixels=std::nullopt,
				 std::optional<int> ypixels=std::nullopt,
				 std::optional<int> hold_setting=std::nullopt);

    /* if you set truncate=true, you WILL have to destroy the Vecs returned from this function. */
    /*
      std::vector<Vec> embeddedTrajectory(bool truncate=false, bool gather_to_root=false);*/
    
    ~TakensEmbedding();

    PetscErrorCode status() const noexcept;

    PetscReal      embeddingInterval() const noexcept;

    PetscReal      timeWindowLength() const noexcept;

    PetscReal      initialTime() const noexcept;

    PetscReal      samplingInterval() const noexcept;

  private:
    PetscReal        tau, t0, tf, dt, T;
    PetscInt         m, n;
    Vec              xdata;
    std::vector<int> evlens;
    std::vector<Vec> embedding;
    Vec              plotx, ploty, plot_buff;
    PetscErrorCode   ierr;
    PetscInt         DataLength();
    /* gets the index of the nearest time and a bool indicating if the 
       index matches the requested time within 10^{-10} or not */
    std::pair<int, bool>  IndexFromTime(PetscReal t);
    
    PetscErrorCode   AllocateEmbeddingVecs();
    PetscErrorCode   FreeEmbeddingVecs();

  };


}
#endif
