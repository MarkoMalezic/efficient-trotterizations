#ifndef _MIN_HELPERS_H_
#define _MIN_HELPERS_H_

#include <iostream>
#include <array>
#include <random>
#include <iomanip>
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_bin_float.hpp>

using quad = boost::multiprecision::cpp_bin_float_quad;
using namespace std;

using Eigen::DiagonalMatrix;
using Eigen::Dynamic;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using MatrixXld = Eigen::Matrix<long double, Dynamic, Dynamic>;
using MatrixXcld = Eigen::Matrix<complex<long double>, Dynamic, Dynamic>;
using MatrixXQ = Eigen::Matrix<quad, Dynamic, Dynamic>;
using MatrixXcQ = Eigen::Matrix<complex<quad>, Dynamic, Dynamic>;

using Eigen::VectorXcd;
using Eigen::VectorXd;
using VectorXld = Eigen::Matrix<long double, Dynamic, 1>;
using VectorXcld = Eigen::Matrix<complex<long double>, Dynamic, 1>;
using VectorXQ = Eigen::Matrix<quad, Dynamic, 1>;
using VectorXcQ = Eigen::Matrix<complex<quad>, Dynamic, 1>;

// Type trait to deduce the matrix type for a given vector type
template <typename Vec>
struct MatrixForVector;

template <>
struct MatrixForVector<VectorXd>
{
  using type = MatrixXd;
};

template <>
struct MatrixForVector<VectorXcd>
{
  using type = MatrixXcd;
};

template <>
struct MatrixForVector<VectorXld>
{
  using type = MatrixXld;
};

template <>
struct MatrixForVector<VectorXcld>
{
  using type = MatrixXcld;
};

template <>
struct MatrixForVector<VectorXQ>
{
  using type = MatrixXQ;
};

template <>
struct MatrixForVector<VectorXcQ>
{
  using type = MatrixXcQ;
};

template <typename Vec>
using MatrixForVector_t = typename MatrixForVector<Vec>::type;

// Type trait to deduce the hessian matrix type for a given vector type
template <typename Vec>
struct HessMatrixForVector;

template <>
struct HessMatrixForVector<VectorXd>
{
  using type = MatrixXd;
};

template <>
struct HessMatrixForVector<VectorXcd>
{
  using type = MatrixXd;
};

template <>
struct HessMatrixForVector<VectorXld>
{
  using type = MatrixXld;
};

template <>
struct HessMatrixForVector<VectorXcld>
{
  using type = MatrixXld;
};

template <>
struct HessMatrixForVector<VectorXQ>
{
  using type = MatrixXQ;
};

template <>
struct HessMatrixForVector<VectorXcQ>
{
  using type = MatrixXQ;
};

template <typename Vec>
using HessMatrixForVector_t = typename HessMatrixForVector<Vec>::type;


// Helper function to compute the Jt * W
template <typename Mat, typename RealT>
Mat compute_JtW(const Mat &J, const DiagonalMatrix<RealT, Dynamic> &W)
{
  if constexpr (
    std::is_same_v<Mat, MatrixXd> ||
    std::is_same_v<Mat, MatrixXld> ||
    std::is_same_v<Mat, MatrixXQ>
  )
  {
    // Compute the transposed Jacobian * Weights
    return J.transpose() * W;
  }
  else if constexpr (
    std::is_same_v<Mat, MatrixXcd> ||
    std::is_same_v<Mat, MatrixXcld> ||
    std::is_same_v<Mat, MatrixXcQ>
  )
  {
    // Compute the adjoint Jacobian * Weights
    return J.adjoint() * W;
  }
}

// Helper function to compute the approximate Hessian matrix from the Jacobian
template <typename Mat, typename Vec>
HessMatrixForVector_t<Vec> compute_hess(const Mat &J, const Mat &JtW)
{
  if constexpr (
    std::is_same_v<Mat, MatrixXd> ||
    std::is_same_v<Mat, MatrixXld> ||
    std::is_same_v<Mat, MatrixXQ>
  )
  {
    // Compute the regular product for real-valued matrices
    return JtW * J;
  }
  else if constexpr (
    std::is_same_v<Mat, MatrixXcd> ||
    std::is_same_v<Mat, MatrixXcld> ||
    std::is_same_v<Mat, MatrixXcQ>
  )
  {
    // Compute the Hermitian product and take the real part
    return (JtW * J).real();
  }
}

// Helper function to compute the Chi2 depending on the Vector type
template <typename Vec, typename RealT>
auto compute_chi2(const Vec &y, const DiagonalMatrix<RealT, Dynamic> &W)
{
  using ReturnType = conditional_t<is_same_v<Vec, VectorXd> || is_same_v<Vec, VectorXcd>, double, long double>;

  if constexpr (
    std::is_same_v<Vec, VectorXd> ||
    std::is_same_v<Vec, VectorXld> ||
    std::is_same_v<Vec, VectorXQ>
  )
  {
    auto chi2{y.transpose() * W * y}; // Regular transpose for real vectors
    return static_cast<ReturnType>(chi2(0, 0));
  }
  else if constexpr (
    std::is_same_v<Vec, VectorXcd> ||
    std::is_same_v<Vec, VectorXcld> ||
    std::is_same_v<Vec, VectorXcQ>
  )
  {
    auto chi2{y.adjoint() * W * y}; // Conjugate transpose for complex vectors
    return static_cast<ReturnType>(chi2(0, 0).real());
  }
}

// Helper function to compute the scalar product, considering complex conjugation
template <typename Vec>
auto abs_dot_prod(const Vec &v1, const Vec &v2)
{
  using ReturnType = conditional_t<is_same_v<Vec, VectorXd> || is_same_v<Vec, VectorXcd>, double, long double>;

  if constexpr (
    std::is_same_v<Vec, VectorXd> ||
    std::is_same_v<Vec, VectorXld> ||
    std::is_same_v<Vec, VectorXQ>
  )
  {
    // Real vector: Use regular transpose
    auto dot_prod{v1.transpose() * v2};
    return static_cast<ReturnType>(abs(dot_prod(0, 0)));
  }
  else if constexpr (
    std::is_same_v<Vec, VectorXcd> ||
    std::is_same_v<Vec, VectorXcld> ||
    std::is_same_v<Vec, VectorXcQ>
  )
  {
    // Complex vector: Use conjugate transpose
    auto dot_prod{v1.adjoint() * v2};
    return static_cast<ReturnType>(abs(dot_prod(0, 0)));
  }
}

// Helper function to randomly sample initial conditions
template <typename Vec, typename Scalar>
array<Vec, 2> sample(const Scalar &mu, const Scalar &sigma, const array<int, 2> &nps)
{
  // Deduce the precision from the vector
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;
  // Set up the random number generator
  random_device rd;  // Seed generator
  mt19937 gen(rd()); // Mersenne Twister engine

  // Declare the normal distributions
  normal_distribution<long double> dist_real;
  normal_distribution<long double> dist_imag;
  // If the Scalar is double or long double initialize only the real part,
  if constexpr (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
  {
    long double mu_ld(static_cast<long double>(mu));
    long double sigma_ld(static_cast<long double>(sigma));
    dist_real = normal_distribution<long double>(mu_ld, sigma_ld);
  }
  else if constexpr (is_same_v<Scalar, quad>)
  {
    long double mu_ld(mu.template convert_to<long double>());
    long double sigma_ld(sigma.template convert_to<long double>());
    dist_real = normal_distribution<long double>(mu_ld, sigma_ld);
  }
  else if constexpr (is_same_v<Scalar, complex<double>> || is_same_v<Scalar, complex<long double>>)
  {
    complex<long double> mu_ld(static_cast<long double>(mu.real()), static_cast<long double>(mu.imag()));
    complex<long double> sigma_ld(static_cast<long double>(sigma.real()), static_cast<long double>(sigma.imag()));
    dist_real = normal_distribution<long double>(mu_ld.real(), sigma_ld.real());
    dist_imag = normal_distribution<long double>(mu_ld.imag(), sigma_ld.imag());
  }
  else if constexpr (is_same_v<Scalar, complex<quad>>)
  {
    complex<long double> mu_ld(mu.real().template convert_to<long double>(), mu.imag().template convert_to<long double>());
    complex<long double> sigma_ld(sigma.real().template convert_to<long double>(), sigma.imag().template convert_to<long double>());
    dist_real = normal_distribution<long double>(mu_ld.real(), sigma_ld.real());
    dist_imag = normal_distribution<long double>(mu_ld.imag(), sigma_ld.imag());
  }
  // Declare the initial conditions
  Vec a_vec(nps[1]);
  Vec b_vec(nps[0]);

  // Sample and fill the initial conditions for a_vec
  for (int i{0}; i < nps[1] - 1; ++i)
  {
    if constexpr (
    is_same_v<Scalar, double> ||
    is_same_v<Scalar, long double> ||
    is_same_v<Scalar, quad>
    )
    {
      long double ai = dist_real(gen);
      a_vec(i) = RealT(ai);
    }
    if constexpr (
    is_same_v<Scalar, complex<double>> ||
    is_same_v<Scalar, complex<long double>> ||
    is_same_v<Scalar, complex<quad>>
    )
    {
      long double real_part = dist_real(gen);
      long double imag_part = dist_imag(gen);
      a_vec(i) = complex<RealT>(real_part, imag_part);
    }
  }
  // Ensure that parameters ai add up to 1
  a_vec[nps[1] - 1] = static_cast<Scalar>(1.0) - a_vec.head(nps[1] - 1).sum();

  // Sample and fill the initial conditions for b_vec
  for (int i{0}; i < nps[0] - 1; ++i)
  {
    if constexpr (
    is_same_v<Scalar, double> ||
    is_same_v<Scalar, long double> ||
    is_same_v<Scalar, quad>
    )
    {
      long double bi = dist_real(gen);
      b_vec(i) = RealT(bi);
    }
    if constexpr (
    is_same_v<Scalar, complex<double>> ||
    is_same_v<Scalar, complex<long double>> ||
    is_same_v<Scalar, complex<quad>>
    )
    {
      long double real_part = dist_real(gen);
      long double imag_part = dist_imag(gen);
      b_vec(i) = complex<RealT>(real_part, imag_part);
    }
  }
  // Ensure that parameters bi add up to 1
  b_vec[nps[0] - 1] = static_cast<Scalar>(1.0) - b_vec.head(nps[0] - 1).sum();

  // Return the array of initial conditions
  return array<Vec, 2>{{a_vec, b_vec}};
}

// Helper function to compare two vectors
template <typename Vec>
bool compare_min(const Vec &v1, const Vec &v2, const typename Eigen::NumTraits<typename Vec::Scalar>::Real &tol);

// Helper to transfrom from the symmetric basis to the standard basis
template <typename Vec>
pair<Vec, Vec> to_standard(const int &q, const Vec &a_vec, const Vec &b_vec);

// Helper to transfrom from the symmetric basis to the ramp basis
template <typename Vec>
pair<Vec, Vec> to_ramp(const int &q, const Vec &a_vec, const Vec &b_vec);


// Class which stores the minimization results
template <typename Vec>
class MinResult
{
public:
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;
  int conv;       // Whether the algorithm converged and which condition was met (0, 1, 2, 3)
  int iter;       // Number of iterations performed
  HessMatrixForVector_t<Vec> hess;  // Final Hessian matrix
  RealT err2;     // Error of the desired order
  Vec a_vec;      // Final a_vec
  Vec b_vec;      // Final b_vec
  Vec da_vec;     // Final step for a_vec
  Vec db_vec;     // Final step for b_vec

  // MinResult constructor
  MinResult(const array<int, 2> &nps);

  // MinResult default constructor
  MinResult();

  // MinResult deconstructor
  ~MinResult();

  // Method to display the minimization result
  void display();
};

#endif // _MIN_HELPERS_H_
