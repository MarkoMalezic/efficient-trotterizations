#ifndef _NSCHEME_H_
#define _NSCHEME_H_

#include <iomanip>
#include "NCoefficients.h"

// Class which calculates the numerical values of the Trotter-Suzuki scheme by recursively iterating from the inside out
// using template to take in any Eigen vector class (VectorXd, VectorXcd, VectorXld, VectorXcld, VectorXQ, VectorXcQ)
template <typename Vec>
class NScheme
{
public:
  // Determine precision type
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;

  // indices tracking the steps taken
  int ind_A{0};
  int ind_B{0};

  // the switch bool checks which step to take (false = step_A, true = step_B)
  const int n;              // Order O(h^(n+1)) of the scheme
  const int q;              // Number of cycles
  NCoefficients<Vec> coefs; // Coefficients of all orders

  // Evaluation vectors a_eval and b_eval
  Vec a_eval;
  Vec b_eval;

  // Check what to print: - 0 = Nothing
  //                      - 1 = Cycles
  //                      - 2 = Coefficient values
  char verbose;

  // Scheme constructor
  NScheme(const int &order, const int &no_cycles,
          const Vec &a_eval, const Vec &b_eval, const char &verbose = 0);

  // Scheme destructor
  ~NScheme();

  // Methods to update scheme
  void step_A();
  void step_B();

  // Method to iterate through the cycles
  void iterate();

  // Method, which resets the scheme, introduces new evaluation vectors a
  // and b and iterates the scheme with them
  void reiterate(const Vec &a_vec, const Vec &b_vec);

  // Method, which calculates the error for the desired order
  RealT err2(const int &order);

  // Method, which calculates the efficiency for the desired order
  RealT eff(const int &order);

  // Method, which displays squared errors of different order (if known)
  void display_err2();

  // Method, which displays efficiencies of different order (if known)
  void display_eff();
};

#endif // _NSCHEME_H_
