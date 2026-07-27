#ifndef _COEFFICIENTS_H_
#define _COEFFICIENTS_H_

#include <filesystem>
#include "Multinomial.h"
#include "Prefactors.h"
#include "Tensor.h"

namespace fs = std::filesystem;

// Helper function to collect the terms from the sums with their prefactors
vector<pair<int, vector<int>>> get_vecpair(const int &ind, const int &pow, const int &sum_pow, const int &np);


// Class which stores the coefficients of different order Tensor classes
template <typename RealT>
class Coefficients
{
public:
  // Prefactors of the symmetric BCH formula
  static const Prefactors<RealT> prefacs;

  // indices tracking the steps taken
  int &ind_A;
  int &ind_B;

  // Boolean tracking which step to take. false = stepA, true = stepB
  bool &swtch;

  // Order n=2
  Tensor<RealT> alpha;
  Tensor<RealT> beta;

  // Order n=4
  vector<Tensor<RealT>> gammas;

  // Order n=6
  vector<Tensor<RealT>> deltas;

  // Order n=8
  vector<Tensor<RealT>> epsilons;

  // Coefficients constructor
  Coefficients(const int &n, const int &q,
               int &ind_A, int &ind_B, bool &swtch);

  // Coefficients load constructor
  Coefficients(const string &directory, const int &n, const int &q,
               int &ind_A, int &ind_B, bool &swtch);

  // Copy constructor
  Coefficients(const Coefficients &other);

  // Copy assignment operator
  Coefficients &operator=(const Coefficients &other);

  // Move assignment operator
  Coefficients &operator=(Coefficients &&other) noexcept;

  // Coefficients deconstructor
  ~Coefficients();

  // Method to add all combinations to a tensor
  void add_combs(Tensor<RealT> &tensor, const RealT &Q,
                 array<int, 2> powsA, array<int, 2> powsB);

  // Method to add all combinations to a tensor using a vector of transformation tensors of lower order
  void add_combs_transf(Tensor<RealT> &tensor, vector<Tensor<RealT>> tensor_transf, const RealT &Q,
                        array<int, 2> powsA, array<int, 2> powsB);

  // Methods to update the alpha coefficient
  void alpha_stepA();
  void alpha_stepB();

  // Methods to update the beta coefficient
  void beta_stepA();
  void beta_stepB();

  // Method to update the gamma coefficients
  void *gammas_step(const int &gamma_ind);

  // Method to update the delta coefficients
  void *deltas_step(const int &delta_ind);

  // Method to update the epsilon coefficients
  void *epsilons_step(const int &epsilon_ind);
};

#endif // _COEFFICIENTS_H_
