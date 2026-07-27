#ifndef _MINIMIZATION_H_
#define _MINIMIZATION_H_

#include "Minim_helpers.h"
#include "Scheme.h"

// Class which minimizes the provided polynomial manifold
template <typename Vec>
class Minimization
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

  // The loaded or computed scheme
  const Scheme<RealT> scheme;
  // Scheme parameters
  int n;
  int q;

  // Matrix of the Tensor derivatives
  vector<vector<Tensor<RealT>>> *derivs;

  // Dimensions
  int m;             // Number of polynomials
  array<int, 2> nps; // Number of parameters

  // Weights
  DiagonalMatrix<RealT, Dynamic> W;

  // Minimization parameters
  int n_iter;           // Number of iterations
  array<RealT, 2> Ls;  // L_up, L_down
  array<RealT, 4> eps; // Convergence parameters

  // Minimization constructor
  Minimization(const Vec &a_init, const Vec &b_init,
               const Scheme<RealT> &scheme, const vector<RealT> &W_vec,
               const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls = {9.0, 11.0}, const bool &bderivs = true);

  // Overloaded Minimization constructor
  Minimization(const Scheme<RealT> &scheme, const vector<RealT> &W_vec,
               const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls = {9.0, 11.0}, const bool &bderivs = true);

  // Minimization deconstructor
  ~Minimization();

  // Method to construct and evaluate the vector y
  Vec y_eval(const Vec &a_eval, const Vec &b_eval);

  // Method to construct and evaluate the Jacobian
  Mat jacobian();

  // Method to minimize the polynomial manifold
  MinResult<Vec> minimize(RealT lambda, const HessMatrixForVector_t<Vec> *hess = nullptr, const bool &verbose = false);

  // Method to minimize the polynomial manifold in two steps
  MinResult<Vec> min_twostep(RealT lambda, const array<RealT, 4> *eps2 = nullptr, const bool &verbose = false, const bool &freeze = false);

  // Method to find as many minima of the polynomial manifold
  pair<array<vector<Vec>, 4>, vector<int>> find(const int &N, RealT &lambda, const int &steps,
                                                const Scalar &mu, const Scalar &sigma,
                                                const array<RealT, 4> *eps2 = nullptr, const bool &bsort = true, const bool &verbose = false,
                                                const RealT &tol = RealT(1e-15), const bool &freeze = false);
};

#endif // _MINIMIZATION_H_
