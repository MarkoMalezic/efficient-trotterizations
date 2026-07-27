#ifndef _NMINIMIZATION_H_
#define _NMINIMIZATION_H_

#include "Minim_helpers.h"
#include "NScheme.h"

// Class which finds the optimal NScheme
template <typename Vec>
class NMinimization
{
public:
  // Deduce the Precision type
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;
  // Deduce scalar type from Vec
  using Scalar = typename Vec::value_type;
  // Deduce matrix type from Vec
  using Mat = MatrixForVector_t<Vec>;

  // Vectors of the current ai and bi values
  Vec a_vec;
  Vec b_vec;

  // Scheme parameters
  int n;
  int q;

  // Dimensions
  int m;             // Number of polynomials
  array<int, 2> nps; // Number of parameters

  // Weights
  DiagonalMatrix<RealT, Dynamic> W;

  // Derivative parameters
  RealT step;

  // Minimization parameters
  int n_iter;          // Number of iterations
  array<RealT, 2> Ls;  // L_up, L_down
  array<RealT, 4> eps; // Convergence parameters

  // NMinimization constructor
  NMinimization(const Vec &a_init, const Vec &b_init, const RealT step,
                const int &n, const int &q, const vector<RealT> &W_vec,
                const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls = {9.0, 11.0});

  // Overloaded NMinimization constructor
  NMinimization(const RealT step, const int &n, const int &q, const vector<RealT> &W_vec,
                const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls = {9.0, 11.0});

  // NMinimization default constructor
  NMinimization();

  // NMinimization deconstructor
  ~NMinimization();

  // Method to construct the vector y
  Vec y_eval(const Vec &a_eval, const Vec &b_eval);

  // Method to construct the vector y with additional term (distance from the origin)
  Vec y_eval_origin(const Vec &a_eval, const Vec &b_eval, const RealT &ratio);

  // Method to construct the Jacobian
  Mat jacobian();

  // Method to construct the Jacobian with the additional term (distance from the origin)
  Mat jacobian_origin(const RealT &ratio);

  // Method to minimize the polynomial manifold
  MinResult<Vec> minimize(RealT lambda, const HessMatrixForVector_t<Vec> *hess = nullptr, const bool &verbose = false);

  // Method to minimize the polynomial manifold with the additional term (distance from the origin)
  MinResult<Vec> minimize_origin(RealT lambda, const RealT ratio, const HessMatrixForVector_t<Vec> *hess = nullptr, const bool &verbose = false);

  // Method to minimize the polynomial manifold in two steps
  MinResult<Vec> min_twostep(RealT lambda, const array<RealT, 4> *eps2 = nullptr, const bool &verbose = false, const bool &freeze = false);

  // Method to minimize the polynomial manifold in two steps with the additional term (distance from the origin)
  MinResult<Vec> min_twostep_origin(RealT lambda, const RealT ratio, const array<RealT, 4> *eps2 = nullptr, const bool &verbose = false, const bool &freeze = false);

  // Method to find as many minima of the polynomial manifold
  pair<array<vector<Vec>, 4>, vector<int>> find(const int &N, RealT &lambda, const int &steps,
                                                const Scalar &mu, const Scalar &sigma,
                                                const array<RealT, 4> *eps2 = nullptr, const bool &bsort = true, const bool &verbose = false,
                                                const RealT &tol = RealT(1e-15), const bool &freeze = false);
};

#endif // _NMINIMIZATION_H_
