#include <nlts/mutual.h>
#include <cmath>
#include <algorithm>

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

  
  

  static PetscErrorCode ShannonEntropy(PetscInt tau, long num_part,
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
    for(i=tau; i < static_cast<PetscInt>(arr.size()); ++i){
      auto hfront = arr[i];
      auto hlag = arr[i - tau];
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
    
      
  static PetscErrorCode ShannonEntropyAndGradTau(PetscInt tau, PetscInt tau_stride, long num_part,
						 const std::vector<long>& arr,
						 std::vector<PetscInt>& hx,
						 std::vector<PetscInt>& hy,
						 std::vector<std::vector<PetscInt>>& hxy,
						 PetscReal *Htau, PetscReal *Htaunext,
						 PetscReal *Hgradtau)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = ShannonEntropy(tau, num_part, arr, hx, hy, hxy, Htau);CHKERRQ(ierr);
    ierr = ShannonEntropy(tau + tau_stride, num_part, arr, hx, hy, hxy, Htaunext);CHKERRQ(ierr);
    *Hgradtau = (*Htaunext - *Htau)/(PetscReal)tau_stride;
    PetscFunctionReturn(ierr);
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
    taumax = max_tau.value_or(20);

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

  static PetscErrorCode ComputeTauUpdateFromGrad(PetscInt tau, PetscReal htau, PetscReal grad_tau, PetscInt *domain_min, PetscInt *domain_max,
						 PetscInt *next_tau)
  {
    PetscFunctionBeginUser;
    if(grad_tau < 0.0){
      *domain_min = std::max(*domain_min, tau);
    } else if(grad_tau > 0.0){
      *domain_max = std::min(*domain_max, tau);
    } else {
      *next_tau = tau;
      PetscFunctionReturn(0);
    }
    *next_tau = tau - (PetscInt)htau/grad_tau;
    if(*next_tau > *domain_max){
      *next_tau = *domain_max;
    } else if(*next_tau < *domain_min){
      *next_tau = *domain_min;
    }

    PetscFunctionReturn(0);
  }

  static PetscErrorCode ExploreAdjacentTaus(PetscInt tau, PetscReal *htau,
					    long num_part,
					    const std::vector<long>& arr,
					    std::vector<PetscInt>& hx,
					    std::vector<PetscInt>& hy,
					    std::vector<std::vector<PetscInt>>& hxy,
					    PetscReal *hlow, PetscReal *hhigh, bool *found_local_min)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    if(tau == 0){
      SETERRQ(PETSC_COMM_WORLD, 1, "Error, ExploreAdjacentTaus called at tau=0, this makes no sense! You should probably not be algorithmically minimizing the mutual entropy of your data. Try plotting it first.\n");
    }
    ierr = ShannonEntropy(tau, num_part, arr, hx, hy, hxy, htau);CHKERRQ(ierr);
    ierr = ShannonEntropy(tau-1, num_part, arr, hx, hy, hxy, hlow);CHKERRQ(ierr);
    ierr = ShannonEntropy(tau+1, num_part, arr, hx, hy, hxy, hhigh);CHKERRQ(ierr);
    
    *found_local_min = (*htau < *hlow and *htau < *hhigh);

    PetscFunctionReturn(0);
  }
   
    
    
					 

  std::pair<PetscInt, PetscErrorCode> MinimizeMutualInformation(Vec X,
								std::optional<PetscInt> max_tau,
								std::optional<PetscInt> tau_grad_stride,
								std::optional<PetscInt> partition_boxes,
								std::optional<PetscReal> tau_grad_min)
  {
    PetscErrorCode ierr;
    std::vector<long> arr;
    PetscInt num_part, nx, tau, tau_next, hmin_tau, domain_min, taumax, tau_step;
    PetscReal tau_grad, htau, htau_next, hmin, grad_min, htleft, htright, min, interval;
    bool      at_local_min=false;

    PetscFunctionBeginUser;
    ierr = VecGetLocalSize(X, &nx);
    num_part = partition_boxes.value_or(16);
    taumax = std::min(nx-1, max_tau.value_or(400));
    tau_step = tau_grad_stride.value_or(2);
    grad_min = tau_grad_min.value_or(0.05);

    ierr = RescaleData(X, nx, min, interval);
    ierr = PartitionArray(X, num_part, nx, arr);
    auto [hx, hy, hxy] = AllocateEntropyVecs(num_part);
    tau_grad = 10 * grad_min;/* initialize with any value larger than grad_min */
    tau=1;
    domain_min=1;
    while(PetscAbsReal(tau_grad) > grad_min){
      /* take first-order Newton steps until derivative is small */
      /*PetscPrintf(PETSC_COMM_WORLD, "Doing Newton steps, tau=%d, tau_grad=%f.\n", tau, tau_grad);*/
      ierr = ShannonEntropyAndGradTau(tau, tau_step, num_part, arr, hx, hy, hxy, &htau, &htau_next, &tau_grad);
      ierr = ComputeTauUpdateFromGrad(tau, htau, tau_grad, &domain_min, &taumax, &tau_next);
      /*PetscPrintf(PETSC_COMM_WORLD, "htau=%f, tau_next=%d.\n", htau, tau_next);*/
      if(tau == tau_next){
	break;
      } else {
	tau = tau_next;
      }
    }
    /*PetscPrintf(PETSC_COMM_WORLD, "Breaking Newton loop.\n");*/
    while(!at_local_min){
      /* explore tau space one step at a time */
      ierr = ExploreAdjacentTaus(tau, &htau, num_part, arr, hx, hy, hxy, &htleft, &htright, &at_local_min);
      /*PetscPrintf(PETSC_COMM_WORLD, "htau=%f, tau=%d.\n", htau, tau);*/
      if(std::abs(htright - htau) < 1.0e-12 and (std::abs(htleft - htau) < 1.0e-12)){
	at_local_min = true;
      }
      if(at_local_min){
	break;
      }
      /*PetscPrintf(PETSC_COMM_WORLD, "htleft=%f, htright=%f, htau=%f\n", htleft, htright, htau);*/
	
      if(htleft < htright){
	  tau -= 1;
      } else if(htright < htleft){
	  tau += 1;
      } else {
	break;
      }
	
    }
  

    std::pair<PetscInt, PetscErrorCode> retval{tau, ierr};
    PetscFunctionReturn(retval);
  }
      
    

    
}

