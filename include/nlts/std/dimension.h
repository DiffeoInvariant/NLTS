#ifndef NLTS_STD_DIMENSION_H
#define NLTS_STD_DIMENSION_H
#include <unordered_map>
#include <array>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>

namespace nlts
{

  template<typename T, std::size_t N>
  struct array_hash
  {
    auto operator() (const std::array<T, N>& arr) const
    {
      std::hash<T> hasher;
      std::size_t res=N;
      for(auto i=0; i<N; ++i){
	res = (res << 1) ^ hasher(arr[i]);
      }
      return res;
    }
  };

  /* type param T: the type of the elements of points in your trajectory.
     a trajectory is a vector of pairs of (point, time), where `point` is
     an std::array<T, N> and `time` is a T. T should be arithmetic,
     but this requirement is not currently expressed with std::enable_if<>
     because the relative simplicity of the class makes diagnosing type errors
     easier without it.

     type param N: how many T's make u a point? N of them!

     The BoxCounter<T, N> class takes a trajectory (as represented by
     a BoxCounter<T, N>::Trajectory_t) and a box size of epsilon, partitions
     the minimal N-box enclosing the entire trajectory into N-cubes of side
     length epsilon, and counts how many boxes have at least one point in 
     them.
     
    */
  template<typename T, std::size_t N>
  class BoxCounter
  {
  public:
    using Trajectory_t = std::vector<std::pair<std::array<T, N>, T>>;

    BoxCounter(const Trajectory_t& traj, T epsilon) : trajectory{traj}, eps{epsilon}, box_counts{traj.size()}
    {
      compute_bounds();
      add_all_points();
      count_nz_box();
    };

    /* how many boxes have at least one point in them? */
    long nonzeroBoxes() const noexcept
    {
      return nnzbox;
    }

    /* the bounds of the minimal N-box enclosing the points */
    std::array<std::pair<T, T>,N> bounds() const noexcept
    {
      return min_max;
    }

    /* the box-counting dimension of the trajectory */
    T boxCountingDimension() const
    {
      return -std::log(static_cast<double>(nnzbox)) / std::log(eps);
    }
    

  private:
    using Key_t = std::array<long, N>;
    
    Trajectory_t                   trajectory;
    T                              eps;
    std::unordered_map<Key_t, int, array_hash<long,N>> box_counts;
    long                           nnzbox;
    std::array<std::pair<T, T>,N>  min_max;

    void compute_bounds()
    {
      for(auto& mm : min_max){
	mm.first = std::numeric_limits<T>::max();
	mm.second = -std::numeric_limits<T>::max();
      }

      for(const auto& xt : trajectory){
	const auto& [x, t] = xt;
	for(auto i=0; i<N; ++i){
	  if(x[i] < min_max[i].first){
	    min_max[i].first = x[i];
	  }
	  if(x[i] > min_max[i].second){
	    min_max[i].second = x[i];
	  }
	}
      }
    }

    long index1d(T point, std::size_t dim)
    {
      return static_cast<long>(std::floor((point - min_max[dim].first)/eps));
    }

    Key_t index(const std::array<T, N>& point)
    {
      Key_t id;
      for(auto i=0; i<N; ++i){
	id[i] = index1d(point[i], i);
      }
      return id;
    }

    void add_point(const std::array<T, N>& point)
    {
      box_counts.insert_or_assign(index(point), 1);
    }

    void add_all_points()
    {
      for(const auto& xt : trajectory){
	const auto& [x, t] = xt;
	add_point(x);
      }
    }

    void count_nz_box()
    {
      nnzbox = 0;
      for(const auto& it : box_counts){
	nnzbox += it.second;
      }
    }
  };
	    

}
#endif
