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

    int nonzeroBoxes() const noexcept
    {
      return nnzbox;
    }

    auto bounds() const noexcept
    {
      return min_max;
    }
    

  private:
    using Key_t = std::array<long, N>;
    
    Trajectory_t                   trajectory;
    T                              eps;
    std::unordered_map<Key_t, int, array_hash<long,N>> box_counts;
    int                            nnzbox;
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
