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
      
      ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, x.size(), X);CHKERRQ(ierr);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, t.size(), T);CHKERRQ(ierr);
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

    PetscErrorCode Initialize(int *argc, char ***argv,
			      const char file[],
			      const char help[])
    {
      PetscFunctionBeginUser;
      auto ierr = PetscInitialize(argc, argv, file, help);CHKERRQ(ierr);
      PetscFunctionReturn(ierr);
    }

    PetscErrorCode Finalize()
    {
      PetscFunctionBeginUser;
      auto ierr = PetscFinalize();
      PetscFunctionReturn(ierr);
    }


    bool has_petsc_option(std::string name,
			  std::optional<std::string> prepend,
			  PetscOptions opts_db)
    {
      PetscBool has_opt;
      if(prepend){
	auto ierr = PetscOptionsHasName(opts_db, prepend->c_str(),
					name.c_str(), &has_opt);CHKERRQ(ierr);
      } else {
	auto ierr = PetscOptionsHasName(opts_db, NULL,
					name.c_str(), &has_opt);CHKERRQ(ierr);
      }
      return has_opt;
    }

  
    template<>
    std::pair<PetscReal, bool> get_petsc_option<PetscReal>(std::string name,
							   std::optional<std::string> prepend,
							   PetscOptions opts_db)
    {
      PetscBool has_opt;
      PetscReal val=0;
      if(prepend){
	auto ierr = PetscOptionsGetReal(opts_db, prepend->c_str(),
					name.c_str(), &val, &has_opt); 
      } else {
	auto ierr = PetscOptionsGetReal(opts_db, NULL,
					name.c_str(), &val, &has_opt); 
      }
      
      return {val, has_opt};
    }

    template<>
    std::pair<PetscInt, bool> get_petsc_option<PetscInt>(std::string name,
							 std::optional<std::string> prepend,
							 PetscOptions opts_db)
    {
      PetscBool has_opt;
      PetscInt val=0;
      if(prepend){
	auto ierr = PetscOptionsGetInt(opts_db, prepend->c_str(),
				       name.c_str(), &val, &has_opt); 
      } else {
	auto ierr = PetscOptionsGetInt(opts_db, NULL,
				       name.c_str(), &val, &has_opt); 
      }
      
      return {val, has_opt};
    }

    template<>
    std::pair<std::string, bool> get_petsc_option<std::string>(std::string name,
							       std::optional<std::string> prepend,
							       PetscOptions opts_db)
    {
      PetscBool has_opt;
      char val[100];

      if(prepend){
	auto ierr = PetscOptionsGetString(opts_db, prepend->c_str(),
					  name.c_str(), val, 100, &has_opt); 
      } else {
	auto ierr = PetscOptionsGetString(opts_db, NULL,
					  name.c_str(), val, 100, &has_opt); 
      }

      return {std::string{val}, has_opt};

    }


#define PETRBF_MAX_OPT_ARR_LEN 10000
    template<>
    std::pair<std::vector<PetscInt>, bool>
    get_petsc_option<std::vector<PetscInt>>(std::string name,
					    std::optional<std::string> prepend,
					    PetscOptions opts_db)
    {
      PetscBool has_opt;
      PetscInt size = PETRBF_MAX_OPT_ARR_LEN;
      PetscInt *vals;
      if(prepend){
	auto ierr = PetscOptionsGetIntArray(opts_db, prepend->c_str(), name.c_str(),
					    vals, &size, &has_opt); 
      } else {
	auto ierr = PetscOptionsGetIntArray(opts_db, NULL, name.c_str(),
					    vals, &size, &has_opt); 
      }
      
      if(has_opt){
	return std::make_pair(std::vector(vals, vals+size), true);
      } else {
	return std::make_pair(std::vector<PetscInt>(), false);
      }
    }

    template<>
    std::pair<std::vector<PetscReal>, bool>
    get_petsc_option<std::vector<PetscReal>>(std::string name,
					     std::optional<std::string> prepend,
					     PetscOptions opts_db)
    {
      PetscBool has_opt;
      PetscInt size = PETRBF_MAX_OPT_ARR_LEN;
      PetscReal *vals;
      if(prepend){
	auto ierr = PetscOptionsGetRealArray(opts_db, prepend->c_str(), name.c_str(),
					     vals, &size, &has_opt); 
      } else {
	auto ierr = PetscOptionsGetRealArray(opts_db, NULL, name.c_str(),
					     vals, &size, &has_opt); 
      }
      if(has_opt){
	auto retval = std::make_pair(std::vector(vals, vals+size), true);
	if(vals){
	  delete[] vals;
	  vals = NULL;
	}
	return retval;
      } else {
	return std::make_pair(std::vector<PetscReal>(), false);
      }
  }

    

  }
}
