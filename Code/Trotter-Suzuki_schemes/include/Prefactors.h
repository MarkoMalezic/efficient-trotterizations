#ifndef _PREFACTORS_H_
#define _PREFACTORS_H_

#include <vector>
#include <type_traits>
#include <boost/multiprecision/cpp_bin_float.hpp>

using quad = boost::multiprecision::cpp_bin_float_quad;
using namespace std;

// Class which stores the prefactors needed to compute Coefficients
// using template for the Tensor to take in any class (double, ...)
template <typename RealT>
class Prefactors
{
public:
  // Order n=2
  static const RealT alpha;
  static const RealT beta;
  // Order n=4
  static const vector<RealT> gammas;
  // Order n=6
  static const vector<RealT> deltas;
  // Order n=8
  static const vector<RealT> epsilons;

  // Prefactors default constructor
  Prefactors();

  // Prefactors deconstructor
  ~Prefactors();
};

#endif // _PREFACTORS_H_
