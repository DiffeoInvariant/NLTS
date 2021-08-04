#include <nlts/std/dimension.h>
#include <odeint/std/rkfehlberg.hpp>
#include <odeint/std/rk4.hpp>
#include <cmath>
#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
using State = std::array<double, 3>;

State Lorenz63RHS(double t, const State& x, double a, double r, double b)
{
  State xprime = x;
  xprime[0] = a * (x[1] - x[0]);
  xprime[1] = r * x[0] - x[1] - x[0] * x[2];
  xprime[2] = x[0] * x[1] - b * x[2];
  return xprime;
}

double bcdim(int N, double eps)
{
  return std::log(static_cast<double>(N))/(std::log(1.0/eps));
}

int main()
{
  namespace ph = std::placeholders;
  double a=16.0,r=45.0,b=4.0, dt=0.001, tf=100.0;

  State x0{10.0, -5.0, 2.0};

  auto rhs = std::bind(Lorenz63RHS, ph::_1, ph::_2, a, r, b);
  auto trajectory = odeint::rk4_integrate(rhs, x0, dt, 0.0, tf);

  std::vector<double> eps{0.5, 1.0e-1, 7.5e-2, 5.0e-2, 2.5e-2, 1.0e-2, 8.5e-3, 7.5e-3,
			  6.5e-3,  5.0e-3, 3.5e-3, 2.5e-3,1.5e-3, 1.0e-3,2.0e-4,
			  5.0e-4, 7.5e-4, 1.0e-4, 5.0e-5, 1.0e-5, 1.0e-6,
			  1.0e-7, 1.0e-8};
  namespace mp = boost::multiprecision;
  std::vector<int> N;
  std::vector<mp::uint1024_t> memusage;
  std::vector<double> bcd;

  for(const auto& e : eps){
    nlts::BoxCounter<double, 3> bctr(trajectory, e);
    /* memory that would be used for the naive calculation */
    N.push_back(bctr.nonzeroBoxes());
    auto bounds = bctr.bounds();
    mp::uint1024_t mem = 1;
    for(const auto& bd : bounds){
      mem *= static_cast<mp::uint1024_t>(sizeof(int)*std::floor((bd.second - bd.first)/e));
    }
    memusage.push_back(mem);
    bcd.push_back(bcdim(bctr.nonzeroBoxes(), e));
  }

  for(auto i=0; i<eps.size(); ++i){
    std::cout << "With epsilon=" << eps[i] << ", N(epsilon) = "
	      << N[i] << ", and log(N(epsilon))/log(1/epsilon) = "
	      << bcd[i] << ".\n"
	      << "If we had used the naive calculation, we would have used "
	      << std::setprecision(5) << std::scientific
	      << mp::cpp_dec_float_50(memusage[i]) / 1'000'000 << " MB of memory.\n";
  }

  return 0;
}
