#ifndef _NCOEFFICIENTS_H_
#define _NCOEFFICIENTS_H_

#include <iostream>
#include <stdexcept>
#include "Prefactors.h"
#include <Eigen/Dense>

using Eigen::VectorXcd;
using Eigen::VectorXd;
using VectorXld = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using VectorXcld = Eigen::Matrix<complex<long double>, Eigen::Dynamic, 1>;
using VectorXQ = Eigen::Matrix<quad, Eigen::Dynamic, 1>;
using VectorXcQ = Eigen::Matrix<complex<quad>, Eigen::Dynamic, 1>;

// Class which stores the numerical values of coefficients
// using template to take in any Eigen vector class (VectorXd, VectorXcd, VectorXld, VectorXcld)
template <typename Vec>
class NCoefficients
{
public:
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;
  using Scalar = typename Vec::value_type;

  // Prefactors of the symmetric BCH formula
  static const Prefactors<RealT> prefacs;

  // Evaluation vectors a_eval and b_eval
  Vec a_eval;
  Vec b_eval;

  // Order n=0
  Scalar nu{};
  Scalar sigma{};

  // Order n=2
  Scalar alpha{};
  Scalar beta{};

  // Order n=4
  vector<Scalar> gammas;

  // Order n=6
  vector<Scalar> deltas;

  // Order n=8
  vector<Scalar> epsilons;

  // NCoefficients constructor
  NCoefficients(const Vec &a_eval, const Vec &b_eval);

  // NCoefficients default constructor
  NCoefficients();

  // NCoefficients deconstructor
  ~NCoefficients();

  // Methods to update the alpha coefficient
  void alpha_stepA(int &ind_A);
  void alpha_stepB(int &ind_B);

  // Methods to update the beta coefficient
  void beta_stepA(int &ind_A);
  void beta_stepB(int &ind_B);

  // Method to update the gamma coefficients
  void *gammas_step(const int &gamma_ind, int &ind);

  // Method to update the delta coefficients
  void *deltas_step(const int &delta_ind, int &ind);

  // Method to update the epsilon coefficients
  void *epsilons_step(const int &epsilon_ind, int &ind);
};

#endif // _NCOEFFICIENTS_H_
