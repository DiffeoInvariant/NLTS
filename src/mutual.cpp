#include <nlts/mutual.h>

namespace nlts
{

  static PetscErrorCode RescaleData(Vec X, PetscInt n, double& min, double& interval)
  {
    PetscInt i;
    PetscScalar *x;
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = VecGetArray(X, &x);CHKERRQ(ierr);
    min = x[0];
    interval = x[0];
    for(i=0; i<n; ++i){
      if(x[i] < min){
	min = x[i];
      }
      if(x[i] > interval){
	interval=x[i];
      }
    }
    interval -= min;
    if(interval == 0.0){
      SETERRQ(PETSC_COMM_WORLD, 1, "Error, the data is constant! Cannot continue with Shannon entropy calculation.\n");
    }
    for(i=0; i<n; ++i){
      x[i] = (x[i] - min)/interval;
    }
    ierr = VecRestoreArray(X, &x);CHKERRQ(ierr);
    PetscFunctionReturn(ierr);
  }
      
  
  static PetscErrorCode PartitionArray(Vec X, PetscInt num_part, PetscInt& nx, std::vector<long>& arr)
  {
    PetscErrorCode ierr;
    PetscInt       n;
    long           i;
    const PetscScalar *x;
    PetscFunctionBeginUser;
    ierr = VecGetLocalSize(X, &n);CHKERRQ(ierr);
    nx = n;
    arr.resize(n);
    ierr = VecGetArrayRead(X, &x);CHKERRQ(ierr);
    for(i=0; i < n; ++i){
      if(x[i] < 1.0){
	arr[i] = (long)(x[i] * (PetscReal)num_part);
      } else {
	arr[i] = num_part - 1;
      }
    }
    ierr = VecRestoreArrayRead(X, &x);CHKERRQ(ierr);
    PetscFunctionReturn(ierr);
  }

  static std::tuple<std::vector<PetscInt>,
		    std::vector<PetscInt>,
		    std::vector<std::vector<PetscInt>> >
  AllocateEntropyVecs(PetscInt num_part)
  {
    std::vector<PetscInt> hx, hy;
    std::vector<std::vector<PetscInt>> hxy;
    hx.resize(num_part);
    hy.resize(num_part);
    hxy.resize(num_part);
    for(auto& h : hxy){
      h.resize(num_part);
    }
    return std::make_tuple(hx, hy, hxy);
  }

  static PetscErrorCode ShannonEntropy(PetscInt idx, long num_part,
				       const std::vector<long>& arr,
				       std::vector<PetscInt>& hx,
				       std::vector<PetscInt>& hy,
				       std::vector<std::vector<PetscInt>>& hxy,
				       PetscReal *H)
  {
    long i, j, count=0;
    PetscReal shannon_norm, px, py, pxy=0.0;
    PetscFunctionBeginUser;
    for(i=0; i < num_part; ++i){
      hx[i] = 0;
      hy[i] = 0;
      for(j=0; j < num_part; ++j){
	hxy[i][j] = 0;
      }
    }
    for(i=idx; i < static_cast<PetscInt>(arr.size()); ++i){
      auto hfront = arr[i];
      auto hlag = arr[i - idx];
      hx[hlag]++;
      hy[hfront]++;
      hxy[hlag][hfront]++;
      count++;
    }
    shannon_norm = 1.0/(PetscReal)count;
    *H = 0.0;
    for(i=0; i < num_part; ++i){
      px = (PetscReal)(hx[i]) * shannon_norm;
      if(px > 0.0){
	for(j=0; j < num_part; ++j){
	  py = (PetscReal)(hy[j]) * shannon_norm;
	  if(py > 0.0){
	    pxy = (PetscReal)(hxy[i][j]) * shannon_norm;
	    if(pxy > 0.0){
	      *H += pxy * PetscLogReal(pxy / (px * py));
	    }
	  }
	}
      }
    }

    PetscFunctionReturn(0);
  }
    
      
    
  
  std::pair<std::optional<std::vector<PetscReal>>,PetscErrorCode>
  MutualInformation(Vec X, std::optional<PetscInt> max_tau,
		    std::optional<PetscInt> partition_boxes,
		    std::optional<std::string> outfile,
		    bool return_info)
  {
    PetscErrorCode ierr;
    std::vector<long> arr;
    PetscInt num_part, nx, tau, taumax;
    PetscReal min, interval;
    std::vector<PetscReal> entropy;
    
    if(partition_boxes){
      num_part = *partition_boxes;
    } else {
      num_part = 16;
    }
    if(max_tau){
      taumax = *max_tau;
    } else {
      taumax = 20;
    }

    ierr = VecGetLocalSize(X, &nx);
    ierr = RescaleData(X, nx, min, interval);
    ierr = PartitionArray(X, num_part, nx, arr);
    
    if(taumax >= nx){
      taumax = nx - 1;
    }
    entropy.resize(taumax);
    auto [hx, hy, hxy] = AllocateEntropyVecs(num_part);
    for(tau=0; tau<taumax; ++tau){
      ierr = ShannonEntropy(tau, num_part, arr, hx, hy, hxy, &entropy[tau]);
    }

    if(outfile){
      std::ofstream of(*outfile);
      auto precision = std::numeric_limits<double>::digits10;
      
      of << "Lag " << "Shannon_Entropy\n";
      of << std::setprecision(precision);
      for(auto i=0; i<entropy.size(); ++i){
	of << i << " " << entropy[i] << '\n';
      }
    }
      

    if(return_info){
      auto pair = std::make_pair(entropy, ierr);
      PetscFunctionReturn(pair);
    } else {
      auto pair = std::make_pair(std::nullopt, ierr);
      PetscFunctionReturn(pair);
    }
  }
}
