#include <nlts/io.h>
#include <algorithm>
#include <numeric>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace nlts
{
  inline namespace io
  {
    std::pair<std::vector<double>, std::vector<double>> read_scalar_trajectory(std::string filename)
    {
      double x, t;
      std::pair<std::vector<double>, std::vector<double>> data;
      std::ifstream infile(filename);
      while(infile >> x >> t){
        data.first.push_back(x);
	data.second.push_back(t);
      }
      return data;
    }


    PetscErrorCode VecReadScalarTrajectory(std::string filename,
					   Vec *X, Vec *T)
    {
      PetscErrorCode ierr;
      PetscMPIInt    rank, size;
      std::vector<PetscInt> idx;
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      MPI_Comm_size(PETSC_COMM_WORLD, &size);

      if(size > 1 and rank){
	PetscPrintf(PETSC_COMM_WORLD, "WARNING: creating vectors on multiple processes or on a non-root process! This is probably not gonna work.\n");
      }
      
      auto [x, t] = read_scalar_trajectory(filename);
      idx.resize(x.size());
      std::iota(idx.begin(), idx.end(), 0);
      
      ierr = VecCreateShared(PETSC_COMM_WORLD, PETSC_DECIDE, x.size(), X);CHKERRQ(ierr);
      ierr = VecCreateShared(PETSC_COMM_WORLD, PETSC_DECIDE, t.size(), T);CHKERRQ(ierr);
      ierr = VecSetValues(*X, x.size(), idx.data(), x.data(), INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(*T, t.size(), idx.data(), t.data(), INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(*X);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(*X);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(*T);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(*T);CHKERRQ(ierr);
    

      return 0;
    }

    PetscErrorCode VecWriteScalarTrajectory(std::string filename_pre, Vec X, Vec T, bool make_new_directory)
    {
      PetscErrorCode ierr;
      PetscViewer    viewer;
      std::string    xflname, tflname, ndname;
      namespace bfs = boost::filesystem;

      xflname = filename_pre + std::string{"_x.dat"};
      tflname = filename_pre + std::string{"_t.dat"};
      if(make_new_directory){
	bfs::path new_dir_candidate = bfs::current_path();
	auto dir_base = new_dir_candidate;
        ndname = std::string{"trajectory_data"};
	new_dir_candidate /= ndname;
	bool okay_to_mkdir = false;
	int  n = 0;
	while(!okay_to_mkdir){
	  if(bfs::exists(bfs::status(new_dir_candidate))){
	    if(n == 0){
	      ndname += "_0";
	    }
	    ++n;
	    ndname.pop_back();
	    if(n > 10){
	      ndname.pop_back();
	      if(n > 100){
		ndname.pop_back();
	      }
	    }
	    ndname += std::to_string(n);
	    new_dir_candidate = dir_base;
	    new_dir_candidate /= ndname;
	  } else {
	    okay_to_mkdir = true;
	  }
	}

	auto dirname =  new_dir_candidate.string();
	auto syscall = std::string{"mkdir "} + dirname;
	std::system(syscall.c_str());
        xflname.insert(0, std::string(dirname + std::string{"/"}));
	tflname.insert(0, std::string(dirname + std::string{"/"}));
      }
	    
	  
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, xflname.c_str(), FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
      ierr = VecView(X, viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, tflname.c_str(), FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
      ierr = VecView(T, viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

      PetscFunctionReturn(ierr);
    }


    PetscErrorCode VecReadBinary(std::string filename, Vec *x)
    {
      PetscErrorCode ierr;
      PetscViewer    viewer;
      ierr = VecCreate(PETSC_COMM_WORLD, x);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer);CHKERRQ(ierr);
      ierr = VecLoad(*x, viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
      return 0;
    }

  }
}
