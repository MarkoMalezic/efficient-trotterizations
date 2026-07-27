#include "Coefficients.h"

#include <limits>

// Define the static member for double specialization
template <>
const Prefactors<double> Coefficients<double>::prefacs = Prefactors<double>();

// Define the static member for long double specialization
template <>
const Prefactors<long double> Coefficients<long double>::prefacs = Prefactors<long double>();

// Helper function to collect the terms from the sums with their prefactors
// vector<prefactor, combination>
// e.g. a_ind^pow * eta_ind^sum_pow
vector<pair<int, vector<int>>> get_vecpair(const int &ind, const int &pow, const int &sum_pow, const int &np)
{
  vector<pair<int, vector<int>>> vecpair{};
  // If the index on the sum is 0 then there is no possible combination
  if (ind == 0 && sum_pow > 0)
  {
    return vecpair;
  }
  // For null powers the prefactor is simply 1
  if (pow == 0 && sum_pow == 0)
  {
    pair<int, vector<int>> pair_term(1, vector<int>(np, 0));
    vecpair.push_back(pair_term);
    // If only the powers of the sum are zero then we add the powers of the single parameter and the prefactor is still 1
  }
  else if (pow > 0 && sum_pow == 0)
  {
    vector<int> comb(np, 0);
    comb[ind] = pow;
    pair<int, vector<int>> pair_term(1, comb); // The prefactor is simply 1
    vecpair.push_back(pair_term);
    // If the powers of the sum aren't zero we need to find all the terms in the multinomial.
    // If the powers of the parameter are zero, then we don't add them to the combination
  }
  else if (pow == 0 && sum_pow > 0)
  {
    Multinomial multinom(sum_pow, np, ind); // Calculate all the combinations and their prefactors
    for (int pairi{0}; pairi < multinom.prefacs.size(); ++pairi)
    {
      pair<int, vector<int>> pair_term(multinom.prefacs[pairi], multinom.combs[pairi]);
      vecpair.push_back(pair_term);
    }
    // Otherwise we also add the powers of the parameter to every combination
  }
  else
  {
    Multinomial multinom(sum_pow, np, ind); // Calculate all the combinations and their prefactors
    for (int pairi{0}; pairi < multinom.prefacs.size(); ++pairi)
    {
      vector<int> combsi{multinom.combs[pairi]};
      combsi[ind] += pow;
      pair<int, vector<int>> pair_term(multinom.prefacs[pairi], combsi);
      vecpair.push_back(pair_term);
    }
  }
  return vecpair;
}

// Constructor for Coefficients
template <typename RealT>
Coefficients<RealT>::Coefficients(const int &n, const int &q,
                                  int &ind_A, int &ind_B, bool &swtch)
    : ind_A(ind_A), ind_B(ind_B), swtch(swtch),
      alpha(q, {3, 1}), beta(q, {3, 2}), gammas(6), deltas(18), epsilons(56)
{
  // Initialization of higher order only if desired
  if (n >= 4)
  {
    for (int i = 0; i < gammas.size(); ++i)
    {
      array<int, 2> coef5i = {5, static_cast<int>(i + 1)};
      gammas[i] = Tensor<RealT>(q, coef5i);
    }
  }
  if (n >= 6)
  {
    for (int i = 0; i < deltas.size(); ++i)
    {
      array<int, 2> coef7i = {7, static_cast<int>(i + 1)};
      deltas[i] = Tensor<RealT>(q, coef7i);
    }
  }
  if (n >= 8)
  {
    cout << "Initializing epsilons..." << endl;
    for (int i = 0; i < epsilons.size(); ++i)
    {
      cout << "Epsilon " << i + 1 << endl;
      array<int, 2> coef9i = {9, static_cast<int>(i + 1)};
      epsilons[i] = Tensor<RealT>(q, coef9i);
    }
  }
}

// Load constructor for Coefficients
template <typename RealT>
Coefficients<RealT>::Coefficients(const string &directory, const int &n, const int &q,
                                  int &ind_A, int &ind_B, bool &swtch)
    : ind_A(ind_A), ind_B(ind_B), swtch(swtch),
      alpha(q, {3, 1}), beta(q, {3, 2}), gammas(6), deltas(18), epsilons(56)
{
  int order{};
  if (q < 3)
  {
    order = 2;
  }
  else if (q < 7)
  {
    order = 4;
  }
  else if (q < 17)
  {
    order = 6;
  }
  else
  {
    order = 8;
  }

  string dir{directory + "n" + to_string(order) + "_q" + to_string(q) + "/"};

  if (!fs::exists(dir))
  {
    throw runtime_error("No loadable data found.");
  }

  alpha.load(dir, "alpha");
  beta.load(dir, "beta");

  // Initialization of higher order only if desired
  if (n >= 4)
  {
    for (int i = 0; i < gammas.size(); ++i)
    {
      array<int, 2> coefi = {5, static_cast<int>(i + 1)};
      gammas[i] = Tensor<RealT>(q, coefi);
      gammas[i].load(dir, "gamma" + to_string(i + 1));
    }
  }
  if (n >= 6)
  {
    for (int i = 0; i < deltas.size(); ++i)
    {
      array<int, 2> coefi = {7, static_cast<int>(i + 1)};
      deltas[i] = Tensor<RealT>(q, coefi);
      deltas[i].load(dir, "delta" + to_string(i + 1));
    }
  }
  if (n >= 8)
  {
    for (int i = 0; i < epsilons.size(); ++i)
    {
      array<int, 2> coefi = {9, static_cast<int>(i + 1)};
      epsilons[i] = Tensor<RealT>(q, coefi);
      epsilons[i].load(dir, "epsilon" + to_string(i + 1));
    }
  }
}

// Copy constructor for Coefficients
template <typename RealT>
Coefficients<RealT>::Coefficients(const Coefficients<RealT> &other)
    : ind_A(other.ind_A), ind_B(other.ind_B), swtch(other.swtch),
      alpha(other.alpha), beta(other.beta), gammas(other.gammas), deltas(other.deltas), epsilons(other.epsilons)
{
}

// Copy assignment operator
template <typename RealT>
Coefficients<RealT> &Coefficients<RealT>::operator=(const Coefficients<RealT> &other)
{
  if (this != &other)
  {
    // Copy the members from the other object
    // Assuming Coefficients has members like alpha, beta, gammas, deltas, epsilons, etc.
    alpha = other.alpha;
    beta = other.beta;
    gammas = other.gammas;
    deltas = other.deltas;
    epsilons = other.epsilons;
    // Copy other members as needed
  }
  return *this;
}

// Move assignment operator
template <typename RealT>
Coefficients<RealT> &Coefficients<RealT>::operator=(Coefficients<RealT> &&other) noexcept
{
  if (this != &other)
  {
    // Move the members from the other object
    alpha = std::move(other.alpha);
    beta = std::move(other.beta);
    gammas = std::move(other.gammas);
    deltas = std::move(other.deltas);
    epsilons = std::move(other.epsilons);
    // Move other members as needed
  }
  return *this;
}

// Destructor for Coefficients
template <typename RealT>
Coefficients<RealT>::~Coefficients()
{
}

// Method to add all combinations to a desired Tensor
template <typename RealT>
void Coefficients<RealT>::add_combs(Tensor<RealT> &tensor, const RealT &Q,
                                    array<int, 2> pows_A, array<int, 2> pows_B)
{
  // Declare the vector of pairs
  vector<pair<int, vector<int>>> vecpair_A{};
  vector<pair<int, vector<int>>> vecpair_B{};
  if (swtch)
  {
    vecpair_A = get_vecpair(ind_A, pows_B[0], pows_B[1], tensor.metadata.nps[1]);
    vecpair_B = get_vecpair(ind_B, pows_A[0], pows_A[1], tensor.metadata.nps[0]);
  }
  else
  {
    vecpair_A = get_vecpair(ind_A, pows_A[0], pows_A[1], tensor.metadata.nps[1]);
    vecpair_B = get_vecpair(ind_B, pows_B[0], pows_B[1], tensor.metadata.nps[0]);
  }
  // Iterate through all terms and add them to the tensor
  for (int termA{0}; termA < vecpair_A.size(); ++termA)
  {
    pair<int, vector<int>> pairA = vecpair_A[termA];
    for (int termB{0}; termB < vecpair_B.size(); ++termB)
    {
      pair<int, vector<int>> pairB = vecpair_B[termB];

      // Generate the Powers class
      vector<int> a_pow(pairA.second);
      vector<int> b_pow(pairB.second);
      a_pow.resize(tensor.metadata.nps[1]);
      b_pow.resize(tensor.metadata.nps[0]);
      Powers powers(a_pow, b_pow);

      // Generate the prefactor
      RealT prefac = Q * pairA.first * pairB.first;

      // Add the new term to the existing coefficient
      array<const int, 2> lex_inds{powers.get_lex_inds(tensor.metadata.dims, tensor.metadata.pows)};
      tensor.T[lex_inds[0]][lex_inds[1]] += prefac;
    }
  }
}

// Method to add all combinations to a tensor using a vector of transformation tensors of lower order
// The switch between vecpair_A and vecpair_B for step B is done in the method call
template <typename RealT>
void Coefficients<RealT>::add_combs_transf(Tensor<RealT> &tensor, vector<Tensor<RealT>> transf_tensors, const RealT &Q,
                                           array<int, 2> pows_A, array<int, 2> pows_B)
{
  // Declare the vector of pairs
  vector<pair<int, vector<int>>> vecpair_A{};
  vector<pair<int, vector<int>>> vecpair_B{};
  if (swtch)
  {
    vecpair_A = get_vecpair(ind_A, pows_B[0], pows_B[1], tensor.metadata.nps[1]);
    vecpair_B = get_vecpair(ind_B, pows_A[0], pows_A[1], tensor.metadata.nps[0]);
  }
  else
  {
    vecpair_A = get_vecpair(ind_A, pows_A[0], pows_A[1], tensor.metadata.nps[1]);
    vecpair_B = get_vecpair(ind_B, pows_B[0], pows_B[1], tensor.metadata.nps[0]);
  }

  // Initialize the Tensor to be transformed
  Tensor<RealT> multi_tens{transf_tensors[0]};

  // Multiply the Tensors before adding the combinations to the final solution (If there are more than 1 Tensors to be transformed)
  for (size_t ti{1}; ti < transf_tensors.size(); ++ti)
  {
    multi_tens = multi_tens.multiply(transf_tensors[ti]);
  }

  // Add all the combinations
  for (size_t termA{0}; termA < vecpair_A.size(); ++termA)
  {
    pair<int, vector<int>> pairA = vecpair_A[termA];
    for (size_t termB{0}; termB < vecpair_B.size(); ++termB)
    {
      pair<int, vector<int>> pairB = vecpair_B[termB];

      // Generate the Powers class
      vector<int> a_pow(pairA.second);
      vector<int> b_pow(pairB.second);
      a_pow.resize(tensor.metadata.nps[1]);
      b_pow.resize(tensor.metadata.nps[0]);
      Powers powers(a_pow, b_pow);

      // Generate the prefactor
      RealT prefac = Q * pairA.first * pairB.first;

      // Add the transformed tensor to the desired tensor
      tensor += multi_tens.transform(tensor.metadata, powers, prefac);
    }
  }
}

// Method to update the alpha coefficient with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename RealT>
void Coefficients<RealT>::alpha_stepA()
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;

  // First term
  array<int, 2> pows_A1{{2, 0}};
  array<int, 2> pows_B1{{0, 1}};
  add_combs(alpha, alpha_Q, pows_A1, pows_B1);

  // Second term
  array<int, 2> pows_A2{{1, 1}};
  array<int, 2> pows_B2{{0, 1}};
  add_combs(alpha, -beta_Q, pows_A2, pows_B2);
}

// Method to update the alpha coefficient with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename RealT>
void Coefficients<RealT>::alpha_stepB()
{
  // Prefactor initialization
  RealT beta_Q = prefacs.beta;

  // First term
  array<int, 2> pows_A1{{1, 0}};
  array<int, 2> pows_B1{{0, 2}};
  add_combs(alpha, beta_Q, pows_A1, pows_B1);
}

// Method to update the beta coefficient with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename RealT>
void Coefficients<RealT>::beta_stepA()
{
  // Prefactor initialization
  RealT beta_Q = prefacs.beta;

  // First term
  array<int, 2> pows_A1{{1, 0}};
  array<int, 2> pows_B1{{0, 2}};
  add_combs(beta, beta_Q, pows_A1, pows_B1);
}

// Method to update the beta coefficient with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename RealT>
void Coefficients<RealT>::beta_stepB()
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;

  // First term
  array<int, 2> pows_A1{{2, 0}};
  array<int, 2> pows_B1{{0, 1}};
  add_combs(beta, alpha_Q, pows_A1, pows_B1);

  // Second term
  array<int, 2> pows_A2{{1, 1}};
  array<int, 2> pows_B2{{0, 1}};
  add_combs(beta, -beta_Q, pows_A2, pows_B2);
}

// Method to update the gamma coefficients by switching between operators A, B
template <typename RealT>
void *Coefficients<RealT>::gammas_step(const int &gamma_ind)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;
  vector<RealT> gammas_Q = prefacs.gammas;
  switch (gamma_ind)
  {
  case 1:
  {
    // First term
    array<int, 2> pows_A1{{4, 0}};
    array<int, 2> pows_B1{{0, 1}};
    add_combs(gammas[gamma_ind-1], gammas_Q[0], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{3, 1}};
    array<int, 2> pows_B2{{0, 1}};
    add_combs(gammas[gamma_ind-1], gammas_Q[1] + gammas_Q[2], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{2, 2}};
    array<int, 2> pows_B3{{0, 1}};
    add_combs(gammas[gamma_ind-1], -gammas_Q[3] - gammas_Q[4], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 3}};
    array<int, 2> pows_B4{{0, 1}};
    add_combs(gammas[gamma_ind-1], -gammas_Q[5], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{2, 0}};
    array<int, 2> pows_B5{{0, 0}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {beta}, alpha_Q, pows_A5, pows_B5) :
            add_combs_transf(gammas[gamma_ind-1], {alpha}, alpha_Q, pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{1, 1}};
    array<int, 2> pows_B6{{0, 0}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {beta}, -beta_Q, pows_A6, pows_B6) :
            add_combs_transf(gammas[gamma_ind-1], {alpha}, -beta_Q, pows_A6, pows_B6);
    break;
  }
  case 2:
  {
    // First term
    array<int, 2> pows_A1{{3, 0}};
    array<int, 2> pows_B1{{0, 2}};
    add_combs(gammas[gamma_ind-1], gammas_Q[1], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{2, 1}};
    array<int, 2> pows_B2{{0, 2}};
    add_combs(gammas[gamma_ind-1], -2*gammas_Q[3] - gammas_Q[4], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 2}};
    array<int, 2> pows_B3{{0, 2}};
    add_combs(gammas[gamma_ind-1], -2*gammas_Q[5], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 0}};
    array<int, 2> pows_B4{{0, 0}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {alpha}, -alpha_Q, pows_A4, pows_B4) :
            add_combs_transf(gammas[gamma_ind-1], {beta}, -alpha_Q, pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 0}};
    array<int, 2> pows_B5{{0, 1}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {beta}, beta_Q, pows_A5, pows_B5) :
            add_combs_transf(gammas[gamma_ind-1], {alpha}, beta_Q, pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{1, 1}};
    array<int, 2> pows_B6{{0, 0}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {alpha}, beta_Q, pows_A6, pows_B6) :
            add_combs_transf(gammas[gamma_ind-1], {beta}, beta_Q, pows_A6, pows_B6);
    break;
  }
  case 3:
  {
    // First term
    array<int, 2> pows_A1{{3, 0}};
    array<int, 2> pows_B1{{0, 2}};
    add_combs(gammas[gamma_ind-1], gammas_Q[2], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{2, 1}};
    array<int, 2> pows_B2{{0, 2}};
    add_combs(gammas[gamma_ind-1], -gammas_Q[4], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 2}};
    array<int, 2> pows_B3{{0, 2}};
    add_combs(gammas[gamma_ind-1], -gammas_Q[5], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 0}};
    array<int, 2> pows_B4{{0, 1}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {beta}, -2*beta_Q, pows_A4, pows_B4) :
            add_combs_transf(gammas[gamma_ind-1], {alpha}, -2*beta_Q, pows_A4, pows_B4);
    break;
  }
  case 4:
  {
    array<int, 2> pows_A1{{2, 0}};
    array<int, 2> pows_B1{{0, 3}};
    // First term
    add_combs(gammas[gamma_ind-1], gammas_Q[3], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{1, 1}};
    array<int, 2> pows_B2{{0, 3}};
    add_combs(gammas[gamma_ind-1], gammas_Q[5], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 0}};
    array<int, 2> pows_B3{{0, 1}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {alpha}, beta_Q, pows_A3, pows_B3) :
            add_combs_transf(gammas[gamma_ind-1], {beta}, beta_Q, pows_A3, pows_B3);
    break;
  }
  case 5:
  {
    // First term
    array<int, 2> pows_A1{{2, 0}};
    array<int, 2> pows_B1{{0, 3}};
    add_combs(gammas[gamma_ind-1], gammas_Q[4], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{1, 1}};
    array<int, 2> pows_B2{{0, 3}};
    add_combs(gammas[gamma_ind-1], 2*gammas_Q[5], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 0}};
    array<int, 2> pows_B3{{0, 1}};
    swtch ? add_combs_transf(gammas[gamma_ind-1], {alpha}, -2*beta_Q, pows_A3, pows_B3) :
            add_combs_transf(gammas[gamma_ind-1], {beta}, -2*beta_Q, pows_A3, pows_B3);
    break;
  }
  case 6:
  {
    // First term
    array<int, 2> pows_A1{{1, 0}};
    array<int, 2> pows_B1{{0, 4}};
    add_combs(gammas[gamma_ind-1], gammas_Q[5], pows_A1, pows_B1);
    break;
  }
  default:
    throw invalid_argument("Invalid gamma index");
  }
  return nullptr;
}

// Method to update the delta coefficients by switching between operators A, B
template <typename RealT>
void *Coefficients<RealT>::deltas_step(const int &delta_ind)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;
  vector<RealT> gammas_Q = prefacs.gammas;
  vector<RealT> deltas_Q = prefacs.deltas;
  switch (delta_ind)
  {
  case 1:
  {
    // First term
    array<int, 2> pows_A1{{6, 0}};
    array<int, 2> pows_B1{{0, 1}};
    add_combs(deltas[delta_ind-1], deltas_Q[0], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{5, 1}};
    array<int, 2> pows_B2{{0, 1}};
    add_combs(deltas[delta_ind-1], deltas_Q[1] + deltas_Q[2] + deltas_Q[3], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{4, 2}};
    array<int, 2> pows_B3{{0, 1}};
    add_combs(deltas[delta_ind-1], deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{3, 3}};
    array<int, 2> pows_B4{{0, 1}};
    add_combs(deltas[delta_ind-1], -deltas_Q[9] - deltas_Q[10] - deltas_Q[11] - deltas_Q[12] - deltas_Q[13], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{2, 4}};
    array<int, 2> pows_B5{{0, 1}};
    add_combs(deltas[delta_ind-1], -deltas_Q[14] - deltas_Q[15] - deltas_Q[16], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{1, 5}};
    array<int, 2> pows_B6{{0, 1}};
    add_combs(deltas[delta_ind-1], -deltas_Q[17], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{4, 0}};
    array<int, 2> pows_B7{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[0], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[0], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{3, 1}};
    array<int, 2> pows_B8{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[1] + gammas_Q[2], pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[1] + gammas_Q[2], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{2, 2}};
    array<int, 2> pows_B9{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[3] - gammas_Q[4], pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[3] - gammas_Q[4], pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{1, 3}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[5], pows_A10, pows_B10) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[5], pows_A10, pows_B10);

    // Eleventh term
    array<int, 2> pows_A11{{2, 0}};
    array<int, 2> pows_B11{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[5]}, alpha_Q, pows_A11, pows_B11) :
            add_combs_transf(deltas[delta_ind-1], {gammas[0]}, alpha_Q, pows_A11, pows_B11);

    // Twelfth term
    array<int, 2> pows_A12{{1, 1}};
    array<int, 2> pows_B12{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[5]}, -beta_Q, pows_A12, pows_B12) :
            add_combs_transf(deltas[delta_ind-1], {gammas[0]}, -beta_Q, pows_A12, pows_B12);
    break;
  }
  case 2:
  {
    // First term
    array<int, 2> pows_A1{{5, 0}};
    array<int, 2> pows_B1{{0, 2}};
    add_combs(deltas[delta_ind-1], deltas_Q[1], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{4, 1}};
    array<int, 2> pows_B2{{0, 2}};
    add_combs(deltas[delta_ind-1], 3*deltas_Q[4] + 2*deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + 2*deltas_Q[8], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{3, 2}};
    array<int, 2> pows_B3{{0, 2}};
    add_combs(deltas[delta_ind-1], -3*deltas_Q[9] - 4*deltas_Q[10] - 4*deltas_Q[11] - 3*deltas_Q[12] - 2*deltas_Q[13], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 3}};
    array<int, 2> pows_B4{{0, 2}};
    add_combs(deltas[delta_ind-1], -5*deltas_Q[14] - 5*deltas_Q[15] - 4*deltas_Q[16], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 4}};
    array<int, 2> pows_B5{{0, 2}};
    add_combs(deltas[delta_ind-1], -5*deltas_Q[17], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{4, 0}};
    array<int, 2> pows_B6{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[0], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[0], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{3, 0}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[2], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[2], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{3, 1}};
    array<int, 2> pows_B8{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[1] - 2*gammas_Q[2], pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[1] - 2*gammas_Q[2], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{2, 1}};
    array<int, 2> pows_B9{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[4], pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[4], pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{2, 2}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 2*gammas_Q[3] + 2*gammas_Q[4], pows_A10, pows_B10) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 2*gammas_Q[3] + 2*gammas_Q[4], pows_A10, pows_B10);

    // Eleventh term
    array<int, 2> pows_A11{{1, 2}};
    array<int, 2> pows_B11{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[5], pows_A11, pows_B11) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[5], pows_A11, pows_B11);

    // Twelfth term
    array<int, 2> pows_A12{{1, 3}};
    array<int, 2> pows_B12{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 2*gammas_Q[5], pows_A12, pows_B12) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 2*gammas_Q[5], pows_A12, pows_B12);

    // Thirteenth term
    array<int, 2> pows_A13{{2, 0}};
    array<int, 2> pows_B13{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[4]}, 2*alpha_Q, pows_A13, pows_B13) :
            add_combs_transf(deltas[delta_ind-1], {gammas[1]}, 2*alpha_Q, pows_A13, pows_B13);

    // Fourtheenth term
    array<int, 2> pows_A14{{2, 0}};
    array<int, 2> pows_B14{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[3]}, alpha_Q, pows_A14, pows_B14) :
            add_combs_transf(deltas[delta_ind-1], {gammas[2]}, alpha_Q, pows_A14, pows_B14);

    // Fiftheenth term
    array<int, 2> pows_A15{{1, 1}};
    array<int, 2> pows_B15{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[4]}, -2*beta_Q, pows_A15, pows_B15) :
            add_combs_transf(deltas[delta_ind-1], {gammas[1]}, -2*beta_Q, pows_A15, pows_B15);

    // Sixteenth term
    array<int, 2> pows_A16{{1, 1}};
    array<int, 2> pows_B16{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[3]}, -beta_Q, pows_A16, pows_B16) :
            add_combs_transf(deltas[delta_ind-1], {gammas[2]}, -beta_Q, pows_A16, pows_B16);

    // Seventeenth term
    array<int, 2> pows_A17{{1, 0}};
    array<int, 2> pows_B17{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta, beta}, -beta_Q, pows_A17, pows_B17) :
            add_combs_transf(deltas[delta_ind-1], {alpha, alpha}, -beta_Q, pows_A17, pows_B17);
    break;
  }
  case 3:
  {
    // First term
    array<int, 2> pows_A1{{5, 0}};
    array<int, 2> pows_B1{{0, 2}};
    add_combs(deltas[delta_ind-1], deltas_Q[2], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{4, 1}};
    array<int, 2> pows_B2{{0, 2}};
    add_combs(deltas[delta_ind-1], -deltas_Q[4] + deltas_Q[6] - deltas_Q[8], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{3, 2}};
    array<int, 2> pows_B3{{0, 2}};
    add_combs(deltas[delta_ind-1], deltas_Q[10] + 2*deltas_Q[11] + deltas_Q[12], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 3}};
    array<int, 2> pows_B4{{0, 2}};
    add_combs(deltas[delta_ind-1], deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 4}};
    array<int, 2> pows_B5{{0, 2}};
    add_combs(deltas[delta_ind-1], deltas_Q[17], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{4, 0}};
    array<int, 2> pows_B6{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[0], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[0], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{3, 0}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[1] - 2*gammas_Q[2], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[1] - 2*gammas_Q[2], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{3, 1}};
    array<int, 2> pows_B8{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[1] + gammas_Q[2], pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[1] + gammas_Q[2], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{2, 1}};
    array<int, 2> pows_B9{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[3] + gammas_Q[4], pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[3] + gammas_Q[4], pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{2, 2}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[3] - gammas_Q[4], pows_A10, pows_B10) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[3] - gammas_Q[4], pows_A10, pows_B10);

    // Eleventh term
    array<int, 2> pows_A11{{1, 3}};
    array<int, 2> pows_B11{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[5], pows_A11, pows_B11) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[5], pows_A11, pows_B11);

    // Twelfth term
    array<int, 2> pows_A12{{2, 0}};
    array<int, 2> pows_B12{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[4]}, -alpha_Q, pows_A12, pows_B12) :
            add_combs_transf(deltas[delta_ind-1], {gammas[1]}, -alpha_Q, pows_A12, pows_B12);

    // Thirteenth term
    array<int, 2> pows_A13{{1, 0}};
    array<int, 2> pows_B13{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[5]}, beta_Q, pows_A13, pows_B13) :
            add_combs_transf(deltas[delta_ind-1], {gammas[0]}, beta_Q, pows_A13, pows_B13);

    // Fourtheenth term
    array<int, 2> pows_A14{{1, 1}};
    array<int, 2> pows_B14{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[4]}, beta_Q, pows_A14, pows_B14) :
            add_combs_transf(deltas[delta_ind-1], {gammas[1]}, beta_Q, pows_A14, pows_B14);

    // Fiftheenth term
    array<int, 2> pows_A15{{1, 0}};
    array<int, 2> pows_B15{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta, beta}, 2*beta_Q, pows_A15, pows_B15) :
            add_combs_transf(deltas[delta_ind-1], {alpha, alpha}, 2*beta_Q, pows_A15, pows_B15);
    break;
  }
  case 4:
  {
    // First term
    array<int, 2> pows_A1{{5, 0}};
    array<int, 2> pows_B1{{0, 2}};
    add_combs(deltas[delta_ind-1], deltas_Q[3], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{4, 1}};
    array<int, 2> pows_B2{{0, 2}};
    add_combs(deltas[delta_ind-1], deltas_Q[7] + deltas_Q[8], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{3, 2}};
    array<int, 2> pows_B3{{0, 2}};
    add_combs(deltas[delta_ind-1], -deltas_Q[11] - deltas_Q[12] - deltas_Q[13], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 3}};
    array<int, 2> pows_B4{{0, 2}};
    add_combs(deltas[delta_ind-1], -deltas_Q[15] - deltas_Q[16], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 4}};
    array<int, 2> pows_B5{{0, 2}};
    add_combs(deltas[delta_ind-1], -deltas_Q[17], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{3, 0}};
    array<int, 2> pows_B6{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, 2*gammas_Q[2], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, 2*gammas_Q[2], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{2, 1}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[4], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[4], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{1, 2}};
    array<int, 2> pows_B8{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[5], pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[5], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{1, 0}};
    array<int, 2> pows_B9{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[5]}, -2*beta_Q, pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {gammas[0]}, -2*beta_Q, pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{1, 0}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta, beta}, -beta_Q, pows_A10, pows_B10) :
            add_combs_transf(deltas[delta_ind-1], {alpha, alpha}, -beta_Q, pows_A10, pows_B10);
    break;
  }
  case 5:
  {
    // First term
    array<int, 2> pows_A1{{4, 0}};
    array<int, 2> pows_B1{{0, 3}};
    add_combs(deltas[delta_ind-1], deltas_Q[4], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{3, 1}};
    array<int, 2> pows_B2{{0, 3}};
    add_combs(deltas[delta_ind-1], -deltas_Q[9] - 3*deltas_Q[10] - 3*deltas_Q[11] - deltas_Q[12], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{2, 2}};
    array<int, 2> pows_B3{{0, 3}};
    add_combs(deltas[delta_ind-1], -5*deltas_Q[14] - 5*deltas_Q[15] - 3*deltas_Q[16], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 3}};
    array<int, 2> pows_B4{{0, 3}};
    add_combs(deltas[delta_ind-1], -5*deltas_Q[17], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{3, 0}};
    array<int, 2> pows_B5{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[1], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[1], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{2, 0}};
    array<int, 2> pows_B6{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[3] - gammas_Q[4], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[3] - gammas_Q[4], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{2, 1}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[3] - gammas_Q[4], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[3] - gammas_Q[4], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{1, 1}};
    array<int, 2> pows_B8{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -3*gammas_Q[5], pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -3*gammas_Q[5], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{1, 2}};
    array<int, 2> pows_B9{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[5], pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[5], pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{2, 0}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, -3*alpha_Q, pows_A10, pows_B10) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, -3*alpha_Q, pows_A10, pows_B10);

    // Eleventh term
    array<int, 2> pows_A11{{2, 0}};
    array<int, 2> pows_B11{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[1]}, -alpha_Q, pows_A11, pows_B11) :
            add_combs_transf(deltas[delta_ind-1], {gammas[4]}, -alpha_Q, pows_A11, pows_B11);

    // Twelfth term
    array<int, 2> pows_A12{{1, 1}};
    array<int, 2> pows_B12{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, 3*beta_Q, pows_A12, pows_B12) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, 3*beta_Q, pows_A12, pows_B12);

    // Thirteenth term
    array<int, 2> pows_A13{{1, 1}};
    array<int, 2> pows_B13{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[1]}, beta_Q, pows_A13, pows_B13) :
            add_combs_transf(deltas[delta_ind-1], {gammas[4]}, beta_Q, pows_A13, pows_B13);

    // Fourtheenth term
    array<int, 2> pows_A14{{1, 0}};
    array<int, 2> pows_B14{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta, alpha}, beta_Q, pows_A14, pows_B14) :
            add_combs_transf(deltas[delta_ind-1], {alpha, beta}, beta_Q, pows_A14, pows_B14);
    break;
  }
  case 6:
  {
    // First term
    array<int, 2> pows_A1{{4, 0}};
    array<int, 2> pows_B1{{0, 3}};
    add_combs(deltas[delta_ind-1], deltas_Q[5], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{3, 1}};
    array<int, 2> pows_B2{{0, 3}};
    add_combs(deltas[delta_ind-1], -deltas_Q[9] + deltas_Q[10] + 3*deltas_Q[11] - deltas_Q[13], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{2, 2}};
    array<int, 2> pows_B3{{0, 3}};
    add_combs(deltas[delta_ind-1], deltas_Q[14] + 3*deltas_Q[15] + deltas_Q[16], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 3}};
    array<int, 2> pows_B4{{0, 3}};
    add_combs(deltas[delta_ind-1], deltas_Q[17], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{3, 0}};
    array<int, 2> pows_B5{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -3*gammas_Q[1], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -3*gammas_Q[1], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{2, 0}};
    array<int, 2> pows_B6{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, 3*gammas_Q[3] + 2*gammas_Q[4], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, 3*gammas_Q[3] + 2*gammas_Q[4], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{2, 1}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 6*gammas_Q[3] + 3*gammas_Q[4], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 6*gammas_Q[3] + 3*gammas_Q[4], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{1, 1}};
    array<int, 2> pows_B8{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, 7*gammas_Q[5], pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, 7*gammas_Q[5], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{1, 2}};
    array<int, 2> pows_B9{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 6*gammas_Q[5], pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 6*gammas_Q[5], pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{2, 0}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, 3*alpha_Q, pows_A10, pows_B10) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, 3*alpha_Q, pows_A10, pows_B10);

    // Eleventh term
    array<int, 2> pows_A11{{1, 0}};
    array<int, 2> pows_B11{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[4]}, beta_Q, pows_A11, pows_B11) :
            add_combs_transf(deltas[delta_ind-1], {gammas[1]}, beta_Q, pows_A11, pows_B11);

    // Twelfth term
    array<int, 2> pows_A12{{1, 1}};
    array<int, 2> pows_B12{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, -3*beta_Q, pows_A12, pows_B12) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, -3*beta_Q, pows_A12, pows_B12);

    // Thirteenth term
    array<int, 2> pows_A13{{1, 0}};
    array<int, 2> pows_B13{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta, alpha}, -2*beta_Q, pows_A13, pows_B13) :
            add_combs_transf(deltas[delta_ind-1], {alpha, beta}, -2*beta_Q, pows_A13, pows_B13);
    break;
  }
  case 7:
  {
    // First term
    array<int, 2> pows_A1{{4, 0}};
    array<int, 2> pows_B1{{0, 3}};
    add_combs(deltas[delta_ind-1], deltas_Q[6], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{3, 1}};
    array<int, 2> pows_B2{{0, 3}};
    add_combs(deltas[delta_ind-1], -deltas_Q[9] - deltas_Q[10] - deltas_Q[11], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{2, 2}};
    array<int, 2> pows_B3{{0, 3}};
    add_combs(deltas[delta_ind-1], -2*deltas_Q[14] - deltas_Q[15] - deltas_Q[16], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 3}};
    array<int, 2> pows_B4{{0, 3}};
    add_combs(deltas[delta_ind-1], -2*deltas_Q[17], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{3, 0}};
    array<int, 2> pows_B5{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[1] + gammas_Q[2], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[1] + gammas_Q[2], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{2, 0}};
    array<int, 2> pows_B6{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -3*gammas_Q[3], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -3*gammas_Q[3], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{2, 1}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[3] - 2*gammas_Q[4], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[3] - 2*gammas_Q[4], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{1, 1}};
    array<int, 2> pows_B8{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -3*gammas_Q[5], pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -3*gammas_Q[5], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{1, 2}};
    array<int, 2> pows_B9{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -3*gammas_Q[5], pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -3*gammas_Q[5], pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{2, 0}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, -alpha_Q, pows_A10, pows_B10) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, -alpha_Q, pows_A10, pows_B10);

    // Eleventh term
    array<int, 2> pows_A11{{1, 0}};
    array<int, 2> pows_B11{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[3]}, beta_Q, pows_A11, pows_B11) :
            add_combs_transf(deltas[delta_ind-1], {gammas[2]}, beta_Q, pows_A11, pows_B11);

    // Twelfth term
    array<int, 2> pows_A12{{1, 1}};
    array<int, 2> pows_B12{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, beta_Q, pows_A12, pows_B12) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, beta_Q, pows_A12, pows_B12);

    // Thirteenth term
    array<int, 2> pows_A13{{1, 0}};
    array<int, 2> pows_B13{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta, alpha}, -beta_Q, pows_A13, pows_B13) :
            add_combs_transf(deltas[delta_ind-1], {alpha, beta}, -beta_Q, pows_A13, pows_B13);
    break;
  }
  case 8:
  {
    // First term
    array<int, 2> pows_A1{{4, 0}};
    array<int, 2> pows_B1{{0, 3}};
    add_combs(deltas[delta_ind-1], deltas_Q[7], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{3, 1}};
    array<int, 2> pows_B2{{0, 3}};
    add_combs(deltas[delta_ind-1], -deltas_Q[12] - 2*deltas_Q[13], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{2, 2}};
    array<int, 2> pows_B3{{0, 3}};
    add_combs(deltas[delta_ind-1], -deltas_Q[15] - 2*deltas_Q[16], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 3}};
    array<int, 2> pows_B4{{0, 3}};
    add_combs(deltas[delta_ind-1], -3*deltas_Q[17], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{2, 0}};
    array<int, 2> pows_B5{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[4], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[4], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{1, 1}};
    array<int, 2> pows_B6{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -4*gammas_Q[5], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -4*gammas_Q[5], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{1, 0}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[3]}, -2*beta_Q, pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {gammas[2]}, -2*beta_Q, pows_A7, pows_B7);
    break;
  }
  case 9:
  {
    // First term
    array<int, 2> pows_A1{{4, 0}};
    array<int, 2> pows_B1{{0, 3}};
    add_combs(deltas[delta_ind-1], deltas_Q[8], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{3, 1}};
    array<int, 2> pows_B2{{0, 3}};
    add_combs(deltas[delta_ind-1], -2*deltas_Q[11] - deltas_Q[12], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{2, 2}};
    array<int, 2> pows_B3{{0, 3}};
    add_combs(deltas[delta_ind-1], -2*deltas_Q[15] - deltas_Q[16], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 3}};
    array<int, 2> pows_B4{{0, 3}};
    add_combs(deltas[delta_ind-1], -deltas_Q[17], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{3, 0}};
    array<int, 2> pows_B5{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[2], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[2], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{2, 1}};
    array<int, 2> pows_B6{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 2*gammas_Q[4], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 2*gammas_Q[4], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{1, 2}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 2*gammas_Q[5], pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 2*gammas_Q[5], pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{1, 0}};
    array<int, 2> pows_B8{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[4]}, -2*beta_Q, pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {gammas[1]}, -2*beta_Q, pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{1, 0}};
    array<int, 2> pows_B9{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta, alpha}, 2*beta_Q, pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {alpha, beta}, 2*beta_Q, pows_A9, pows_B9);
    break;
  }
  case 10:
  {
    // First term
    array<int, 2> pows_A1{{3, 0}};
    array<int, 2> pows_B1{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[9], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{2, 1}};
    array<int, 2> pows_B2{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[14] - deltas_Q[15], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 2}};
    array<int, 2> pows_B3{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[17], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 0}};
    array<int, 2> pows_B4{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[3] + gammas_Q[4], pows_A4, pows_B4) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[3] + gammas_Q[4], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 0}};
    array<int, 2> pows_B5{{0, 3}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -gammas_Q[5], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -gammas_Q[5], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{2, 0}};
    array<int, 2> pows_B6{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[0]}, -alpha_Q, pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {gammas[5]}, -alpha_Q, pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{1, 0}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[1]}, beta_Q, pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {gammas[4]}, beta_Q, pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{1, 1}};
    array<int, 2> pows_B8{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[0]}, beta_Q, pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {gammas[5]}, beta_Q, pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{1, 0}};
    array<int, 2> pows_B9{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha, alpha}, -beta_Q, pows_A9, pows_B9) :
            add_combs_transf(deltas[delta_ind-1], {beta, beta}, -beta_Q, pows_A9, pows_B9);
    break;
  }
  case 11:
  {
    // First term
    array<int, 2> pows_A1{{3, 0}};
    array<int, 2> pows_B1{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[10], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{2, 1}};
    array<int, 2> pows_B2{{0, 4}};
    add_combs(deltas[delta_ind-1], 3*deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 2}};
    array<int, 2> pows_B3{{0, 4}};
    add_combs(deltas[delta_ind-1], 3*deltas_Q[17], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 0}};
    array<int, 2> pows_B4{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[3], pows_A4, pows_B4) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[3], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 1}};
    array<int, 2> pows_B5{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[5], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[5], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{2, 0}};
    array<int, 2> pows_B6{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[0]}, 2*alpha_Q, pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {gammas[5]}, 2*alpha_Q, pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{1, 0}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, beta_Q, pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, beta_Q, pows_A7, pows_B7);

    // Eigth term
    array<int, 2> pows_A8{{1, 1}};
    array<int, 2> pows_B8{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[0]}, -2*beta_Q, pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {gammas[5]}, -2*beta_Q, pows_A8, pows_B8);
    break;
  }
  case 12:
  {
    // First term
    array<int, 2> pows_A1{{3, 0}};
    array<int, 2> pows_B1{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[11], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{2, 1}};
    array<int, 2> pows_B2{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[15] + deltas_Q[16], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 2}};
    array<int, 2> pows_B3{{0, 4}};
    add_combs(deltas[delta_ind-1], 2*deltas_Q[17], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 0}};
    array<int, 2> pows_B4{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[4], pows_A4, pows_B4) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[4], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 1}};
    array<int, 2> pows_B5{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 2*gammas_Q[5], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 2*gammas_Q[5], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{1, 0}};
    array<int, 2> pows_B6{{0, 3}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, 4*gammas_Q[5], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, 4*gammas_Q[5], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{1, 0}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[2]}, -2*beta_Q, pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {gammas[3]}, -2*beta_Q, pows_A7, pows_B7);
    break;
  }
  case 13:
  {
    // First term
    array<int, 2> pows_A1{{3, 0}};
    array<int, 2> pows_B1{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[12], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{2, 1}};
    array<int, 2> pows_B2{{0, 4}};
    add_combs(deltas[delta_ind-1], 2*deltas_Q[15], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 2}};
    array<int, 2> pows_B3{{0, 4}};
    add_combs(deltas[delta_ind-1], -deltas_Q[17], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 0}};
    array<int, 2> pows_B4{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -4*gammas_Q[4], pows_A4, pows_B4) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -4*gammas_Q[4], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 1}};
    array<int, 2> pows_B5{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -8*gammas_Q[5], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -8*gammas_Q[5], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{1, 0}};
    array<int, 2> pows_B6{{0, 3}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, -8*gammas_Q[5], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, -8*gammas_Q[5], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{1, 0}};
    array<int, 2> pows_B7{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[1]}, -2*beta_Q, pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {gammas[4]}, -2*beta_Q, pows_A7, pows_B7);

    // Seventh term
    array<int, 2> pows_A8{{1, 0}};
    array<int, 2> pows_B8{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha, alpha}, 2*beta_Q, pows_A8, pows_B8) :
            add_combs_transf(deltas[delta_ind-1], {beta, beta}, 2*beta_Q, pows_A8, pows_B8);
    break;
  }
  case 14:
  {
    // First term
    array<int, 2> pows_A1{{3, 0}};
    array<int, 2> pows_B1{{0, 4}};
    add_combs(deltas[delta_ind-1], deltas_Q[13], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{2, 1}};
    array<int, 2> pows_B2{{0, 4}};
    add_combs(deltas[delta_ind-1], 2*deltas_Q[16], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 2}};
    array<int, 2> pows_B3{{0, 4}};
    add_combs(deltas[delta_ind-1], 5*deltas_Q[17], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{2, 0}};
    array<int, 2> pows_B4{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[4], pows_A4, pows_B4) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[4], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{1, 1}};
    array<int, 2> pows_B5{{0, 2}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, 2*gammas_Q[5], pows_A5, pows_B5) :
            add_combs_transf(deltas[delta_ind-1], {beta}, 2*gammas_Q[5], pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{1, 0}};
    array<int, 2> pows_B6{{0, 3}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {beta}, 6*gammas_Q[5], pows_A6, pows_B6) :
            add_combs_transf(deltas[delta_ind-1], {alpha}, 6*gammas_Q[5], pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{1, 0}};
    array<int, 2> pows_B7{{0, 0}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha, alpha}, -beta_Q, pows_A7, pows_B7) :
            add_combs_transf(deltas[delta_ind-1], {beta, beta}, -beta_Q, pows_A7, pows_B7);
    break;
  }
  case 15:
  {
    // First term
    array<int, 2> pows_A1{{2, 0}};
    array<int, 2> pows_B1{{0, 5}};
    add_combs(deltas[delta_ind-1], deltas_Q[14], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{1, 1}};
    array<int, 2> pows_B2{{0, 5}};
    add_combs(deltas[delta_ind-1], deltas_Q[17], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 0}};
    array<int, 2> pows_B3{{0, 3}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, gammas_Q[5], pows_A3, pows_B3) :
            add_combs_transf(deltas[delta_ind-1], {beta}, gammas_Q[5], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{1, 0}};
    array<int, 2> pows_B4{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[0]}, beta_Q, pows_A4, pows_B4) :
            add_combs_transf(deltas[delta_ind-1], {gammas[5]}, beta_Q, pows_A4, pows_B4);
    break;
  }
  case 16:
  {
    // First term
    array<int, 2> pows_A1{{2, 0}};
    array<int, 2> pows_B1{{0, 5}};
    add_combs(deltas[delta_ind-1], deltas_Q[15], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{1, 1}};
    array<int, 2> pows_B2{{0, 5}};
    add_combs(deltas[delta_ind-1], -deltas_Q[17], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 0}};
    array<int, 2> pows_B3{{0, 1}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {gammas[0]}, -2*beta_Q, pows_A3, pows_B3) :
            add_combs_transf(deltas[delta_ind-1], {gammas[5]}, -2*beta_Q, pows_A3, pows_B3);
    break;
  }
  case 17:
  {
    // First term
    array<int, 2> pows_A1{{2, 0}};
    array<int, 2> pows_B1{{0, 5}};
    add_combs(deltas[delta_ind-1], deltas_Q[16], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{1, 1}};
    array<int, 2> pows_B2{{0, 5}};
    add_combs(deltas[delta_ind-1], 5*deltas_Q[17], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{1, 0}};
    array<int, 2> pows_B3{{0, 3}};
    swtch ? add_combs_transf(deltas[delta_ind-1], {alpha}, -2*gammas_Q[5], pows_A3, pows_B3) :
            add_combs_transf(deltas[delta_ind-1], {beta}, -2*gammas_Q[5], pows_A3, pows_B3);
    break;
  }
  case 18:
  {
    // First term
    array<int, 2> pows_A1{{1, 0}};
    array<int, 2> pows_B1{{0, 6}};
    add_combs(deltas[delta_ind-1], deltas_Q[17], pows_A1, pows_B1);
    break;
  }
  default:
    throw invalid_argument("Invalid delta index");
  }
  return nullptr;
}

// Method to update the epsilon coefficients by switching between operators A, B
template <typename RealT>
void *Coefficients<RealT>::epsilons_step(const int &epsilon_ind)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;
  vector<RealT> gammas_Q = prefacs.gammas;
  vector<RealT> deltas_Q = prefacs.deltas;
  vector<RealT> epsilons_Q = prefacs.epsilons;
  // !!! Very memory intensive !!!
  // ** Not implemented fully **
  switch (epsilon_ind)
  {
  case 1:
  {
    // First term
    array<int, 2> pows_A1{{8, 0}};
    array<int, 2> pows_B1{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], epsilons_Q[0], pows_A1, pows_B1);

    // Second term
    array<int, 2> pows_A2{{7, 1}};
    array<int, 2> pows_B2{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], epsilons_Q[1] + epsilons_Q[2] + epsilons_Q[3] + epsilons_Q[4], pows_A2, pows_B2);

    // Third term
    array<int, 2> pows_A3{{6, 2}};
    array<int, 2> pows_B3{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], epsilons_Q[5] + epsilons_Q[6] + epsilons_Q[7] + epsilons_Q[8] + epsilons_Q[9] + epsilons_Q[10] + epsilons_Q[11] + epsilons_Q[12] + epsilons_Q[13], pows_A3, pows_B3);

    // Fourth term
    array<int, 2> pows_A4{{5, 3}};
    array<int, 2> pows_B4{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], epsilons_Q[14] + epsilons_Q[15] + epsilons_Q[16] + epsilons_Q[17] + epsilons_Q[18] + epsilons_Q[19] + epsilons_Q[20] + epsilons_Q[21] + epsilons_Q[22] + epsilons_Q[23] + epsilons_Q[24] + epsilons_Q[25] + epsilons_Q[26] + epsilons_Q[27], pows_A4, pows_B4);

    // Fifth term
    array<int, 2> pows_A5{{4, 4}};
    array<int, 2> pows_B5{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], -(epsilons_Q[28] + epsilons_Q[29] + epsilons_Q[30] + epsilons_Q[31] + epsilons_Q[32] + epsilons_Q[33] + epsilons_Q[34] + epsilons_Q[35] + epsilons_Q[36] + epsilons_Q[37] + epsilons_Q[38] + epsilons_Q[39] + epsilons_Q[40] + epsilons_Q[41]), pows_A5, pows_B5);

    // Sixth term
    array<int, 2> pows_A6{{3, 5}};
    array<int, 2> pows_B6{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], -(epsilons_Q[42] + epsilons_Q[43] + epsilons_Q[44] + epsilons_Q[45] + epsilons_Q[46] + epsilons_Q[47] + epsilons_Q[48] + epsilons_Q[49] + epsilons_Q[50]), pows_A6, pows_B6);

    // Seventh term
    array<int, 2> pows_A7{{2, 6}};
    array<int, 2> pows_B7{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], -(epsilons_Q[51] + epsilons_Q[52] + epsilons_Q[53] + epsilons_Q[54]), pows_A7, pows_B7);

    // Eighth term
    array<int, 2> pows_A8{{1, 7}};
    array<int, 2> pows_B8{{0, 1}};
    add_combs(epsilons[epsilon_ind-1], -epsilons_Q[55], pows_A8, pows_B8);

    // Ninth term
    array<int, 2> pows_A9{{6, 0}};
    array<int, 2> pows_B9{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {beta}, deltas_Q[0], pows_A9, pows_B9) :
            add_combs_transf(epsilons[epsilon_ind-1], {alpha}, deltas_Q[0], pows_A9, pows_B9);

    // Tenth term
    array<int, 2> pows_A10{{5, 1}};
    array<int, 2> pows_B10{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {beta}, deltas_Q[1] + deltas_Q[2] + deltas_Q[3], pows_A10, pows_B10) :
            add_combs_transf(epsilons[epsilon_ind-1], {alpha}, deltas_Q[1] + deltas_Q[2] + deltas_Q[3], pows_A10, pows_B10);

    // Eleventh term
    array<int, 2> pows_A11{{4, 2}};
    array<int, 2> pows_B11{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {beta}, deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8], pows_A11, pows_B11) :
            add_combs_transf(epsilons[epsilon_ind-1], {alpha}, deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8], pows_A11, pows_B11);

    // Twelfth term
    array<int, 2> pows_A12{{3, 3}};
    array<int, 2> pows_B12{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {beta}, -(deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]), pows_A12, pows_B12) :
            add_combs_transf(epsilons[epsilon_ind-1], {alpha}, -(deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]), pows_A12, pows_B12);

    // Thirteenth term
    array<int, 2> pows_A13{{2, 4}};
    array<int, 2> pows_B13{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {beta}, -(deltas_Q[14] + deltas_Q[15] + deltas_Q[16]), pows_A13, pows_B13) :
            add_combs_transf(epsilons[epsilon_ind-1], {alpha}, -(deltas_Q[14] + deltas_Q[15] + deltas_Q[16]), pows_A13, pows_B13);

    // Fourteenth term
    array<int, 2> pows_A14{{1, 5}};
    array<int, 2> pows_B14{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {beta}, -deltas_Q[17], pows_A14, pows_B14) :
            add_combs_transf(epsilons[epsilon_ind-1], {alpha}, -deltas_Q[17], pows_A14, pows_B14);

    // Fifteenth term
    array<int, 2> pows_A15{{4, 0}};
    array<int, 2> pows_B15{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {gammas[5]}, gammas_Q[0], pows_A15, pows_B15) :
            add_combs_transf(epsilons[epsilon_ind-1], {gammas[0]}, gammas_Q[0], pows_A15, pows_B15);

    // Sixteenth term
    array<int, 2> pows_A16{{3, 1}};
    array<int, 2> pows_B16{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {gammas[5]}, gammas_Q[1] + gammas_Q[2], pows_A16, pows_B16) :
            add_combs_transf(epsilons[epsilon_ind-1], {gammas[0]}, gammas_Q[1] + gammas_Q[2], pows_A16, pows_B16);

    // Seventeenth term
    array<int, 2> pows_A17{{2, 2}};
    array<int, 2> pows_B17{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {gammas[5]}, -(gammas_Q[3] + gammas_Q[4]), pows_A17, pows_B17) :
            add_combs_transf(epsilons[epsilon_ind-1], {gammas[0]}, -(gammas_Q[3] + gammas_Q[4]), pows_A17, pows_B17);

    // Eighteenth term
    array<int, 2> pows_A18{{1, 3}};
    array<int, 2> pows_B18{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {gammas[5]}, -gammas_Q[5], pows_A18, pows_B18) :
            add_combs_transf(epsilons[epsilon_ind-1], {gammas[0]}, -gammas_Q[5], pows_A18, pows_B18);

    // Nineteenth term
    array<int, 2> pows_A19{{2, 0}};
    array<int, 2> pows_B19{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {deltas[17]}, alpha_Q, pows_A19, pows_B19) :
            add_combs_transf(epsilons[epsilon_ind-1], {deltas[0]}, alpha_Q, pows_A19, pows_B19);

    // Twentieth term
    array<int, 2> pows_A20{{1, 1}};
    array<int, 2> pows_B20{{0, 0}};
    swtch ? add_combs_transf(epsilons[epsilon_ind-1], {deltas[17]}, -beta_Q, pows_A20, pows_B20) :
            add_combs_transf(epsilons[epsilon_ind-1], {deltas[0]}, -beta_Q, pows_A20, pows_B20);
    break;
  }
  case 2:
  {
    break;
  }
  case 3:
  {
    break;
  }
  case 4:
  {
    break;
  }
  case 5:
  {
    break;
  }
  case 6:
  {
    break;
  }
  case 7:
  {
    break;
  }
  case 8:
  {
    break;
  }
  case 9:
  {
    break;
  }
  case 10:
  {
    break;
  }
  case 11:
  {
    break;
  }
  case 12:
  {
    break;
  }
  case 13:
  {
    break;
  }
  case 14:
  {
    break;
  }
  case 15:
  {
    break;
  }
  case 16:
  {
    break;
  }
  case 17:
  {
    break;
  }
  case 18:
  {
    break;
  }
  case 56:
  {
    // First term
    array<int, 2> pows_A1{{1, 0}};
    array<int, 2> pows_B1{{0, 8}};
    add_combs(epsilons[epsilon_ind-1], epsilons_Q[55], pows_A1, pows_B1);
    break;
  }
  default:
    throw invalid_argument("Invalid epsilon index");
  }
  return nullptr;
}


// Explicit instantiation for the template class
template class Coefficients<double>;
template class Coefficients<long double>;
