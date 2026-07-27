#include "NCoefficients.h"

// Define the static member for double specialization
template <>
const Prefactors<double> NCoefficients<VectorXd>::prefacs = Prefactors<double>();
template <>
const Prefactors<double> NCoefficients<VectorXcd>::prefacs = Prefactors<double>();

// Define the static member for long double specialization
template <>
const Prefactors<long double> NCoefficients<VectorXld>::prefacs = Prefactors<long double>();
template <>
const Prefactors<long double> NCoefficients<VectorXcld>::prefacs = Prefactors<long double>();

// Define the static member for quad specialization
template <>
const Prefactors<quad> NCoefficients<VectorXQ>::prefacs = Prefactors<quad>();
template <>
const Prefactors<quad> NCoefficients<VectorXcQ>::prefacs = Prefactors<quad>();

// Constructor for NCoefficients
template <typename Vec>
NCoefficients<Vec>::NCoefficients(const Vec &a_eval, const Vec &b_eval)
    : a_eval(a_eval), b_eval(b_eval), gammas(6), deltas(18), epsilons(56)
{
}

// Default constructor for NCoefficients
template <typename Vec>
NCoefficients<Vec>::NCoefficients()
{
}

// Destructor for NCoefficients
template <typename Vec>
NCoefficients<Vec>::~NCoefficients()
{
}

// Method to update the alpha coefficient with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename Vec>
void NCoefficients<Vec>::alpha_stepA(int &ind_A)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;

  // First term
  alpha += alpha_Q * pow(a_eval[ind_A], 2) * sigma;
  // Second term
  alpha -= beta_Q * a_eval[ind_A] * nu * sigma;
}

// Method to update the alpha coefficient with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename Vec>
void NCoefficients<Vec>::alpha_stepB(int &ind_B)
{
  // Prefactor initialization
  RealT beta_Q = prefacs.beta;

  // First term
  alpha += beta_Q * b_eval[ind_B] * pow(nu, 2);
}

// Method to update the beta coefficient with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename Vec>
void NCoefficients<Vec>::beta_stepA(int &ind_A)
{
  // Prefactor initialization
  RealT beta_Q = prefacs.beta;

  // First term
  beta += beta_Q * a_eval[ind_A] * pow(sigma, 2);
}

// Method to update the beta coefficient with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename Vec>
void NCoefficients<Vec>::beta_stepB(int &ind_B)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;

  // First term
  beta += alpha_Q * pow(b_eval[ind_B], 2) * nu;
  // Second term
  beta -= beta_Q * b_eval[ind_B] * nu * sigma;
}

// Method to update the gamma coefficients by switching between operators A, B
template <typename Vec>
void *NCoefficients<Vec>::gammas_step(const int &gamma_ind, int &ind)
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
    gammas[0] += gammas_Q[0] * pow(a_eval[ind], 4) * sigma;

    // Second term
    gammas[0] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * sigma;

    // Third term
    gammas[0] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * sigma;

    // Fourth term
    gammas[0] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * sigma;

    // Fifth term
    gammas[0] += alpha_Q * pow(a_eval[ind], 2) * alpha;

    // Sixth term
    gammas[0] -= beta_Q * a_eval[ind] * nu * alpha;
    break;
  }
  case 2:
  {
    // First term
    gammas[1] += gammas_Q[1] * pow(a_eval[ind], 3) * pow(sigma, 2);

    // Second term
    gammas[1] -= (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * pow(sigma, 2);

    // Third term
    gammas[1] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    gammas[1] -= alpha_Q * pow(a_eval[ind], 2) * beta;

    // Fifth term
    gammas[1] += beta_Q * a_eval[ind] * nu * beta;

    // Sixth term
    gammas[1] += beta_Q * a_eval[ind] * sigma * alpha;
    break;
  }
  case 3:
  {
    // First term
    gammas[2] += gammas_Q[2] * pow(a_eval[ind], 3) * pow(sigma, 2);

    // Second term
    gammas[2] -= gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(sigma, 2);

    // Third term
    gammas[2] -= gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    gammas[2] -= 2*beta_Q * a_eval[ind] * sigma * alpha;
    break;
  }
  case 4:
  {
    // First term
    gammas[3] += gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 3);

    // Second term
    gammas[3] += gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 3);

    // Third term
    gammas[3] += beta_Q * a_eval[ind] * sigma * beta;
    break;
  }
  case 5:
  {
    // First term
    gammas[4] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 3);

    // Second term
    gammas[4] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 3);

    // Third term
    gammas[4] -= 2*beta_Q * a_eval[ind] * sigma * beta;
    break;
  }
  case 6:
  {
    // First term
    gammas[5] += gammas_Q[5] * a_eval[ind] * pow(sigma, 4);
    break;
  }
  default:
    throw invalid_argument("Invalid gamma index");
  }
  return nullptr;
}

// Method to update the gamma coefficients by switching between operators A, B
template <typename Vec>
void *NCoefficients<Vec>::deltas_step(const int &delta_ind, int &ind)
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
    deltas[0] += deltas_Q[0] * pow(a_eval[ind], 6) * sigma;

    // Second term
    deltas[0] += (deltas_Q[1] + deltas_Q[2] + deltas_Q[3]) * pow(a_eval[ind], 5) * nu * sigma;

    // Third term
    deltas[0] += (deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(nu, 2) * sigma;

    // Fourth term
    deltas[0] -= (deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 3) * sigma;

    // Fifth term
    deltas[0] -= (deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 4) * sigma;

    // Sixth term
    deltas[0] -= deltas_Q[17] * a_eval[ind] * pow(nu, 5) * sigma;

    // Seventh term
    deltas[0] += gammas_Q[0] * pow(a_eval[ind], 4) * alpha;

    // Eighth term
    deltas[0] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * alpha;

    // Ninth term
    deltas[0] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * alpha;

    // Tenth term
    deltas[0] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * alpha;

    // Eleventh term
    deltas[0] += alpha_Q * pow(a_eval[ind], 2) * gammas[0];

    // Twelfth term
    deltas[0] -= beta_Q * a_eval[ind] * nu * gammas[0];
    break;
  }
  case 2:
  {
    // First term
    deltas[1] += deltas_Q[1] * pow(a_eval[ind], 5) * pow(sigma, 2);

    // Second term
    deltas[1] += (3*deltas_Q[4] + 2*deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + 2*deltas_Q[8]) * pow(a_eval[ind], 4) * nu * pow(sigma, 2);

    // Third term
    deltas[1] -= (3*deltas_Q[9] + 4*deltas_Q[10] + 4*deltas_Q[11] + 3*deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    deltas[1] -= (5*deltas_Q[14] + 5*deltas_Q[15] + 4*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    deltas[1] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    deltas[1] -= 2*gammas_Q[0] * pow(a_eval[ind], 4) * beta;

    // Seventh term
    deltas[1] += gammas_Q[2] * pow(a_eval[ind], 3) * sigma * alpha;

    // Eighth term
    deltas[1] -= 2*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * beta;

    // Ninth term
    deltas[1] -= gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * alpha;

    // Tenth term
    deltas[1] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * beta;

    // Eleventh term
    deltas[1] -= gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * alpha;

    // Twelfth term
    deltas[1] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * beta;

    // Thirteenth term
    deltas[1] += 2*alpha_Q * pow(a_eval[ind], 2) * gammas[1];

    // Fourtheenth term
    deltas[1] += alpha_Q * pow(a_eval[ind], 2) * gammas[2];

    // Fiftheenth term
    deltas[1] -= 2*beta_Q * a_eval[ind] * nu * gammas[1];

    // Sixteenth term
    deltas[1] -= beta_Q * a_eval[ind] * nu * gammas[2];

    // Seventeenth term
    deltas[1] -= beta_Q * a_eval[ind] * pow(alpha, 2);
    break;
  }
  case 3:
  {
    // First term
    deltas[2] += deltas_Q[2] * pow(a_eval[ind], 5) * pow(sigma, 2);

    // Second term
    deltas[2] -= (deltas_Q[4] - deltas_Q[6] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * pow(sigma, 2);

    // Third term
    deltas[2] += (deltas_Q[10] + 2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    deltas[2] += (deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    deltas[2] += deltas_Q[17] * a_eval[ind] * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    deltas[2] += gammas_Q[0] * pow(a_eval[ind], 4) * beta;

    // Seventh term
    deltas[2] += (gammas_Q[1] - 2*gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * alpha;

    // Eighth term
    deltas[2] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * beta;

    // Ninth term
    deltas[2] -= (2*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * alpha;

    // Tenth term
    deltas[2] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * beta;

    // Eleventh term
    deltas[2] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * beta;

    // Twelfth term
    deltas[2] -= alpha_Q * pow(a_eval[ind], 2) * gammas[1];

    // Thirteenth term
    deltas[2] += beta_Q * a_eval[ind] * sigma * gammas[0];

    // Fourtheenth term
    deltas[2] += beta_Q * a_eval[ind] * nu * gammas[1];

    // Fiftheenth term
    deltas[2] += 2*beta_Q * a_eval[ind] * pow(alpha, 2);
    break;
  }
  case 4:
  {
    // First term
    deltas[3] += deltas_Q[3] * pow(a_eval[ind], 5) * pow(sigma, 2);

    // Second term
    deltas[3] += (deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * pow(sigma, 2);

    // Third term
    deltas[3] -= (deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    deltas[3] -= (deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    deltas[3] -= deltas_Q[17] * a_eval[ind] * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    deltas[3] += 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * alpha;

    // Seventh term
    deltas[3] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * alpha;

    // Eighth term
    deltas[3] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * alpha;

    // Ninth term
    deltas[3] -= 2*beta_Q * a_eval[ind] * sigma * gammas[0];

    // Tenth term
    deltas[3] -= beta_Q * a_eval[ind] * pow(alpha, 2);
    break;
  }
  case 5:
  {
    // First term
    deltas[4] += deltas_Q[4] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[4] -= (deltas_Q[9] + 3*deltas_Q[10] + 3*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[4] -= (5*deltas_Q[14] + 5*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[4] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[4] += gammas_Q[1] * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[4] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[4] -= (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Eighth term
    deltas[4] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Ninth term
    deltas[4] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Tenth term
    deltas[4] -= 3*alpha_Q * pow(a_eval[ind], 2) * gammas[3];

    // Eleventh term
    deltas[4] -= alpha_Q * pow(a_eval[ind], 2) * gammas[4];

    // Twelfth term
    deltas[4] += 3*beta_Q * a_eval[ind] * nu * gammas[3];

    // Thirteenth term
    deltas[4] += beta_Q * a_eval[ind] * nu * gammas[4];

    // Fourtheenth term
    deltas[4] += beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 6:
  {
    // First term
    deltas[5] += deltas_Q[5] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[5] -= (deltas_Q[9] - deltas_Q[10] - 3*deltas_Q[11] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[5] += (deltas_Q[14] + 3*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[5] += deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[5] -= 3*gammas_Q[1] * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[5] += (3*gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[5] += (6*gammas_Q[3] + 3*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Eighth term
    deltas[5] += 7*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Ninth term
    deltas[5] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Tenth term
    deltas[5] += 3*alpha_Q * pow(a_eval[ind], 2) * gammas[3];

    // Eleventh term
    deltas[5] += beta_Q * a_eval[ind] * sigma * gammas[1];

    // Twelfth term
    deltas[5] -= 3*beta_Q * a_eval[ind] * nu * gammas[3];

    // Thirteenth term
    deltas[5] -= 2*beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 7:
  {
    // First term
    deltas[6] += deltas_Q[6] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[6] -= (deltas_Q[9] + deltas_Q[10] + deltas_Q[11]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[6] -= (2*deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[6] -= 2*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[6] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[6] -= 3*gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[6] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Eighth term
    deltas[6] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Ninth term
    deltas[6] -= 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Tenth term
    deltas[6] -= alpha_Q * pow(a_eval[ind], 2) * gammas[3];

    // Eleventh term
    deltas[6] += beta_Q * a_eval[ind] * sigma * gammas[2];

    // Twelfth term
    deltas[6] += beta_Q * a_eval[ind] * nu * gammas[3];

    // Thirteenth term
    deltas[6] -= beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 8:
  {
    // First term
    deltas[7] += deltas_Q[7] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[7] -= (deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[7] -= (deltas_Q[15] + 2*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[7] -= 3*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[7] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Sixth term
    deltas[7] -= 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[7] -= 2*beta_Q * a_eval[ind] * sigma * gammas[2];
    break;
  }
  case 9:
  {
    // First term
    deltas[8] += deltas_Q[8] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[8] -= (2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[8] -= (2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[8] -= deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[8] -= 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[8] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Seventh term
    deltas[8] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Eighth term
    deltas[8] -= 2*beta_Q * a_eval[ind] * sigma * gammas[1];

    // Ninth term
    deltas[8] += 2*beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 10:
  {
    // First term
    deltas[9] += deltas_Q[9] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[9] += (deltas_Q[14] - deltas_Q[15]) * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[9] += deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[9] -= (2*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[9] -= gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[9] -= alpha_Q * pow(a_eval[ind], 2) * gammas[5];

    // Seventh term
    deltas[9] += beta_Q * a_eval[ind] * sigma * gammas[4];

    // Eighth term
    deltas[9] += beta_Q * a_eval[ind] * nu * gammas[5];

    // Ninth term
    deltas[9] -= beta_Q * a_eval[ind] * pow(beta, 2);
    break;
  }
  case 11:
  {
    // First term
    deltas[10] += deltas_Q[10] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[10] += (3*deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[10] += 3*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[10] += gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[10] += gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Sixth term
    deltas[10] += 2*alpha_Q * pow(a_eval[ind], 2) * gammas[5];

    // Seventh term
    deltas[10] += beta_Q * a_eval[ind] * sigma * gammas[3];

    // Eighth term
    deltas[10] -= 2*beta_Q * a_eval[ind] * nu * gammas[5];
    break;
  }
  case 12:
  {
    // First term
    deltas[11] += deltas_Q[11] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[11] += (deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[11] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[11] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[11] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[11] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Seventh term
    deltas[11] -= 2*beta_Q * a_eval[ind] * sigma * gammas[3];
    break;
  }
  case 13:
  {
    // First term
    deltas[12] += deltas_Q[12] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[12] += 2*deltas_Q[15] * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[12] -= deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[12] -= 4*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[12] -= 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[12] -= 8*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Seventh term
    deltas[12] -= 2*beta_Q * a_eval[ind] * sigma * gammas[4];

    // Ninth term
    deltas[12] += 2*beta_Q * a_eval[ind] * pow(beta, 2);
    break;
  }
  case 14:
  {
    // First term
    deltas[13] += deltas_Q[13] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[13] += 2*deltas_Q[16] * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[13] += 5*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[13] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[13] += 6*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[13] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Seventh term
    deltas[13] -= beta_Q * a_eval[ind] * pow(beta, 2);
    break;
  }
  case 15:
  {
    // First term
    deltas[14] += deltas_Q[14] * pow(a_eval[ind], 2) * pow(sigma, 5);

    // Second term
    deltas[14] += deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 5);

    // Third term
    deltas[14] += gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * beta;

    // Fourth term
    deltas[14] += beta_Q * a_eval[ind] * sigma * gammas[5];
    break;
  }
  case 16:
  {
    // First term
    deltas[15] += deltas_Q[15] * pow(a_eval[ind], 2) * pow(sigma, 5);

    // Second term
    deltas[15] -= deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 5);

    // Third term
    deltas[15] -= 2*beta_Q * a_eval[ind] * sigma * gammas[5];
    break;
  }
  case 17:
  {
    // First term
    deltas[16] += deltas_Q[16] * pow(a_eval[ind], 2) * pow(sigma, 5);

    // Second term
    deltas[16] += 5*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 5);

    // Third term
    deltas[16] -= 2*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * beta;
    break;
  }
  case 18:
  {
    // First term
    deltas[17] += deltas_Q[17] * a_eval[ind] * pow(sigma, 6);
    break;
  }
  default:
    throw invalid_argument("Invalid delta index");
  }
  return nullptr;
}

// Method to update the epsilon coefficients by switching between operators A, B
template <typename Vec>
void *NCoefficients<Vec>::epsilons_step(const int &epsilon_ind, int &ind)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;
  vector<RealT> gammas_Q = prefacs.gammas;
  vector<RealT> deltas_Q = prefacs.deltas;
  vector<RealT> epsilons_Q = prefacs.epsilons;
  
  switch (epsilon_ind)
  {
  case 1:
  {
    // First term
    epsilons[0] += epsilons_Q[0] * pow(a_eval[ind], 8) * sigma;

    // Second term
    epsilons[0] += (epsilons_Q[1] + epsilons_Q[2] + epsilons_Q[3] + epsilons_Q[4]) * pow(a_eval[ind], 7) * nu * sigma;

    // Third term
    epsilons[0] += (epsilons_Q[5] + epsilons_Q[6] + epsilons_Q[7] + epsilons_Q[8] + epsilons_Q[9] + epsilons_Q[10] + epsilons_Q[11] + epsilons_Q[12] + epsilons_Q[13]) * pow(a_eval[ind], 6) * pow(nu, 2) * sigma;

    // Fourth term
    epsilons[0] += (epsilons_Q[14] + epsilons_Q[15] + epsilons_Q[16] + epsilons_Q[17] + epsilons_Q[18] + epsilons_Q[19] + epsilons_Q[20] + epsilons_Q[21] + epsilons_Q[22] + epsilons_Q[23] + epsilons_Q[24] + epsilons_Q[25] + epsilons_Q[26] + epsilons_Q[27]) * pow(a_eval[ind], 5) * pow(nu, 3) * sigma;

    // Fifth term
    epsilons[0] -= (epsilons_Q[28] + epsilons_Q[29] + epsilons_Q[30] + epsilons_Q[31] + epsilons_Q[32] + epsilons_Q[33] + epsilons_Q[34] + epsilons_Q[35] + epsilons_Q[36] + epsilons_Q[37] + epsilons_Q[38] + epsilons_Q[39] + epsilons_Q[40] + epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 4) * sigma;

    // Sixth term
    epsilons[0] -= (epsilons_Q[42] + epsilons_Q[43] + epsilons_Q[44] + epsilons_Q[45] + epsilons_Q[46] + epsilons_Q[47] + epsilons_Q[48] + epsilons_Q[49] + epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 5) * sigma;

    // Seventh term
    epsilons[0] -= (epsilons_Q[51] + epsilons_Q[52] + epsilons_Q[53] + epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 6) * sigma;

    // Eighth term
    epsilons[0] -= epsilons_Q[55] * a_eval[ind] * pow(nu, 7) * sigma;

    // Ninth term
    epsilons[0] += deltas_Q[0] * pow(a_eval[ind], 6) * alpha;

    // Tenth term
    epsilons[0] += (deltas_Q[1] + deltas_Q[2] + deltas_Q[3]) * pow(a_eval[ind], 5) * nu * alpha;

    // Eleventh term
    epsilons[0] += (deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(nu, 2) * alpha;

    // Twelfth term
    epsilons[0] -= (deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 3) * alpha;

    // Thirteenth term
    epsilons[0] -= (deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 4) * alpha;

    // Fourteenth term
    epsilons[0] -= deltas_Q[17] * a_eval[ind] * pow(nu, 5) * alpha;

    // Fifteenth term
    epsilons[0] += gammas_Q[0] * pow(a_eval[ind], 4) * gammas[0];

    // Sixteenth term
    epsilons[0] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[0];

    // Seventeenth term
    epsilons[0] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[0];

    // Eighteenth term
    epsilons[0] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[0];

    // Nineteenth term
    epsilons[0] += alpha_Q * pow(a_eval[ind], 2) * deltas[0];

    // Twentieth term
    epsilons[0] -= beta_Q * a_eval[ind] * nu * deltas[0];
    break;
  }
  case 2:
  {
    // First term
    epsilons[1] += epsilons_Q[1] * pow(a_eval[ind], 7) * pow(sigma, 2);

    // Second term
    epsilons[1] += (4*epsilons_Q[5] + 3*epsilons_Q[6] + 5*epsilons_Q[7] + 3*epsilons_Q[8] + epsilons_Q[9] + 3*epsilons_Q[11] + epsilons_Q[12]) * pow(a_eval[ind], 6) * nu * pow(sigma, 2);

    // Third term
    epsilons[1] += (10*epsilons_Q[14] + 8*epsilons_Q[15] + 6*epsilons_Q[16] + 8*epsilons_Q[17] + 6*epsilons_Q[18] + 4*epsilons_Q[19] + 3*epsilons_Q[20] + 5*epsilons_Q[21] + 8*epsilons_Q[22] + 6*epsilons_Q[23] + 4*epsilons_Q[24] + 3*epsilons_Q[25] + 5*epsilons_Q[26] + epsilons_Q[27]) * pow(a_eval[ind], 5) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    epsilons[1] -= (13*epsilons_Q[28] + 9*epsilons_Q[29] + 11*epsilons_Q[30] + 10*epsilons_Q[31] + 8*epsilons_Q[32] + 6*epsilons_Q[33] + 9*epsilons_Q[34] + 11*epsilons_Q[35] + 10*epsilons_Q[36] + 8*epsilons_Q[37] + 6*epsilons_Q[38] + 8*epsilons_Q[39] + 6*epsilons_Q[40] + 4*epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    epsilons[1] -= (14*epsilons_Q[42] + 13*epsilons_Q[43] + 11*epsilons_Q[44] + 14*epsilons_Q[45] + 13*epsilons_Q[46] + 11*epsilons_Q[47] + 9*epsilons_Q[48] + 11*epsilons_Q[49] + 10*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    epsilons[1] -= (14*epsilons_Q[51] + 14*epsilons_Q[52] + 14*epsilons_Q[53] + 13*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 5) * pow(sigma, 2);

    // Seventh term
    epsilons[1] -= 14*epsilons_Q[55] * a_eval[ind] * pow(nu, 6) * pow(sigma, 2);

    // Eighth term
    epsilons[1] -= 5*deltas_Q[0] * pow(a_eval[ind], 6) * beta;

    // Ninth term
    epsilons[1] += (deltas_Q[1] + deltas_Q[2]) * pow(a_eval[ind], 5) * sigma * alpha;

    // Tenth term
    epsilons[1] -= 5*(deltas_Q[1] + deltas_Q[2] + deltas_Q[3]) * pow(a_eval[ind], 5) * nu * beta;

    // Eleventh term
    epsilons[1] += (2*deltas_Q[4] + 2*deltas_Q[5] + 2*deltas_Q[6] + deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * alpha;

    // Twelfth term
    epsilons[1] -= 5*(deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(nu, 2) * beta;

    // Thirteenth term
    epsilons[1] -= (3*deltas_Q[9] + 3*deltas_Q[10] + 2*deltas_Q[11] + 2*deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * alpha;

    // Fourteenth term
    epsilons[1] += 5*(deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 3) * beta;

    // Fifteenth term
    epsilons[1] -= (4*deltas_Q[14] + 3*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * alpha;

    // Sixteenth term
    epsilons[1] += 5*(deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 4) * beta;

    // Seventeenth term
    epsilons[1] -= 4*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * alpha;

    // Eighteenth term
    epsilons[1] += 5*deltas_Q[17] * a_eval[ind] * pow(nu, 5) * beta;

    // Nineteenth term
    epsilons[1] += 5*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[1];

    // Twentieth term
    epsilons[1] += 3*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[2];

    // Twenty-first term
    epsilons[1] -= (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[0];

    // Twenty-second term
    epsilons[1] += 5*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[1];

    // Twenty-third term
    epsilons[1] += 3*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[2];

    // Twenty-fourth term
    epsilons[1] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[0];

    // Twenty-fifth term
    epsilons[1] -= 5*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[1];

    // Twenty-sixth term
    epsilons[1] -= 3*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[2];

    // Twenty-seventh term
    epsilons[1] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[0];

    // Twenty-eighth term
    epsilons[1] -= 5*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[1];

    // Twenty-ninth term
    epsilons[1] -= 3*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[2];

    // Thirtieth term
    epsilons[1] += gammas_Q[1] * pow(a_eval[ind], 3) * pow(alpha, 2);

    // Thirty-first term
    epsilons[1] -= (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * pow(alpha, 2);

    // Thirty-second term
    epsilons[1] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(alpha, 2);

    // Thirty-third term
    epsilons[1] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[1];

    // Thirty-fourth term
    epsilons[1] += alpha_Q * pow(a_eval[ind], 2) * deltas[2];

    // Thirty-fifth term
    epsilons[1] -= 3*beta_Q * a_eval[ind] * nu * deltas[1];

    // Thirty-sixth term
    epsilons[1] -= beta_Q * a_eval[ind] * nu * deltas[2];

    // Thirty-seventh term
    epsilons[1] += beta_Q * a_eval[ind] * alpha * gammas[0];
    break;
  }
  case 3:
  {
    // First term
    epsilons[2] += epsilons_Q[2] * pow(a_eval[ind], 7) * pow(sigma, 2);

    // Second term
    epsilons[2] -= (3*epsilons_Q[5] + 2*epsilons_Q[6] + 5*epsilons_Q[7] + 3*epsilons_Q[8] - epsilons_Q[10] + 3*epsilons_Q[11] - epsilons_Q[13]) * pow(a_eval[ind], 6) * nu * pow(sigma, 2);

    // Third term
    epsilons[2] -= (11*epsilons_Q[14] + 8*epsilons_Q[15] + 5*epsilons_Q[16] + 9*epsilons_Q[17] + 6*epsilons_Q[18] + 3*epsilons_Q[19] + 2*epsilons_Q[20] + 5*epsilons_Q[21] + 9*epsilons_Q[22] + 6*epsilons_Q[23] + 3*epsilons_Q[24] + 2*epsilons_Q[25] + 5*epsilons_Q[26]) * pow(a_eval[ind], 5) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    epsilons[2] += (14*epsilons_Q[28] + 9*epsilons_Q[29] + 12*epsilons_Q[30] + 11*epsilons_Q[31] + 8*epsilons_Q[32] + 5*epsilons_Q[33] + 9*epsilons_Q[34] + 12*epsilons_Q[35] + 11*epsilons_Q[36] + 8*epsilons_Q[37] + 5*epsilons_Q[38] + 9*epsilons_Q[39] + 6*epsilons_Q[40] + 3*epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    epsilons[2] += (15*epsilons_Q[42] + 14*epsilons_Q[43] + 11*epsilons_Q[44] + 15*epsilons_Q[45] + 14*epsilons_Q[46] + 11*epsilons_Q[47] + 9*epsilons_Q[48] + 12*epsilons_Q[49] + 11*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    epsilons[2] += (14*epsilons_Q[51] + 14*epsilons_Q[52] + 15*epsilons_Q[53] + 14*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 5) * pow(sigma, 2);

    // Seventh term
    epsilons[2] += 14*epsilons_Q[55] * a_eval[ind] * pow(nu, 6) * pow(sigma, 2);

    // Eighth term
    epsilons[2] += 6*deltas_Q[0] * pow(a_eval[ind], 6) * beta;

    // Ninth term
    epsilons[2] -= (deltas_Q[1] + 2*deltas_Q[2] - deltas_Q[3]) * pow(a_eval[ind], 5) * sigma * alpha;

    // Tenth term
    epsilons[2] += 6*(deltas_Q[1] + deltas_Q[2] + deltas_Q[3]) * pow(a_eval[ind], 5) * nu * beta;

    // Eleventh term
    epsilons[2] -= (deltas_Q[4] + 2*deltas_Q[5] + 3*deltas_Q[6] - deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * alpha;

    // Twelfth term
    epsilons[2] += 6*(deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(nu, 2) * beta;

    // Thirteenth term
    epsilons[2] += (3*deltas_Q[9] + 2*deltas_Q[10] - deltas_Q[11] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * alpha;

    // Fourteenth term
    epsilons[2] -= 6*(deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 3) * beta;

    // Fifteenth term
    epsilons[2] += (3*deltas_Q[14] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * alpha;

    // Sixteenth term
    epsilons[2] -= 6*(deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 4) * beta;

    // Seventeenth term
    epsilons[2] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * alpha;

    // Eighteenth term
    epsilons[2] -= 6*deltas_Q[17] * a_eval[ind] * pow(nu, 5) * beta;

    // Nineteenth term
    epsilons[2] -= 6*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[1];

    // Twentieth term
    epsilons[2] -= 3*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[2];

    // Twenty-first term
    epsilons[2] += (2*gammas_Q[1] + 3*gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[0];

    // Twenty-second term
    epsilons[2] -= 6*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[1];

    // Twenty-third term
    epsilons[2] -= 3*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[2];

    // Twenty-fourth term
    epsilons[2] -= (4*gammas_Q[3] + 5*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[0];

    // Twenty-fifth term
    epsilons[2] += 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[1];

    // Twenty-sixth term
    epsilons[2] += 3*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[2];

    // Twenty-seventh term
    epsilons[2] -= 7*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[0];

    // Twenty-eighth term
    epsilons[2] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[1];

    // Twenty-ninth term
    epsilons[2] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[2];

    // Thirtieth term
    epsilons[2] -= (2*gammas_Q[1] - gammas_Q[2]) * pow(a_eval[ind], 3) * pow(alpha, 2);

    // Thirty-first term
    epsilons[2] += (4*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * pow(alpha, 2);

    // Thirty-second term
    epsilons[2] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(alpha, 2);

    // Thirty-third term
    epsilons[2] -= 3*alpha_Q * pow(a_eval[ind], 2) * deltas[1];

    // Thirty-fourth term
    epsilons[2] += alpha_Q * pow(a_eval[ind], 2) * deltas[3];

    // Thirty-fifth term
    epsilons[2] += 3*beta_Q * a_eval[ind] * nu * deltas[1];

    // Thirty-sixth term
    epsilons[2] -= beta_Q * a_eval[ind] * nu * deltas[3];

    // Thirty-seventh term
    epsilons[2] -= 4*beta_Q * a_eval[ind] * alpha * gammas[0];
    break;
  }
  case 4:
  {
    // First term
    epsilons[3] += epsilons_Q[3] * pow(a_eval[ind], 7) * pow(sigma, 2);

    // Second term
    epsilons[3] += (epsilons_Q[5] + epsilons_Q[6] + 2*epsilons_Q[7] + 2*epsilons_Q[8] + epsilons_Q[9] + epsilons_Q[10] + epsilons_Q[11]) * pow(a_eval[ind], 6) * nu * pow(sigma, 2);

    // Third term
    epsilons[3] += (4*epsilons_Q[14] + 3*epsilons_Q[15] + 2*epsilons_Q[16] + 4*epsilons_Q[17] + 3*epsilons_Q[18] + 2*epsilons_Q[19] + 2*epsilons_Q[20] + 3*epsilons_Q[21] + 3*epsilons_Q[22] + 2*epsilons_Q[23] + epsilons_Q[24] + epsilons_Q[25] + 2*epsilons_Q[26] + epsilons_Q[27]) * pow(a_eval[ind], 5) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    epsilons[3] -= (5*epsilons_Q[28] + 4*epsilons_Q[29] + 5*epsilons_Q[30] + 5*epsilons_Q[31] + 4*epsilons_Q[32] + 3*epsilons_Q[33] + 3*epsilons_Q[34] + 4*epsilons_Q[35] + 4*epsilons_Q[36] + 3*epsilons_Q[37] + 2*epsilons_Q[38] + 4*epsilons_Q[39] + 3*epsilons_Q[40] + 2*epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    epsilons[3] -= (6*epsilons_Q[42] + 6*epsilons_Q[43] + 5*epsilons_Q[44] + 5*epsilons_Q[45] + 5*epsilons_Q[46] + 4*epsilons_Q[47] + 4*epsilons_Q[48] + 5*epsilons_Q[49] + 5*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    epsilons[3] -= (6*epsilons_Q[51] + 5*epsilons_Q[52] + 6*epsilons_Q[53] + 6*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 5) * pow(sigma, 2);

    // Seventh term
    epsilons[3] -= 6*epsilons_Q[55] * a_eval[ind] * pow(nu, 6) * pow(sigma, 2);

    // Eighth term
    epsilons[3] -= 2*deltas_Q[0] * pow(a_eval[ind], 6) * beta;

    // Ninth term
    epsilons[3] += (deltas_Q[1] + 2*deltas_Q[2] - 2*deltas_Q[3]) * pow(a_eval[ind], 5) * sigma * alpha;

    // Tenth term
    epsilons[3] -= 2*(deltas_Q[1] + deltas_Q[2] + deltas_Q[3]) * pow(a_eval[ind], 5) * nu * beta;

    // Eleventh term
    epsilons[3] += (deltas_Q[4] + 2*deltas_Q[5] + 3*deltas_Q[6] - deltas_Q[7] - 2*deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * alpha;

    // Twelfth term
    epsilons[3] -= 2*(deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(nu, 2) * beta;

    // Thirteenth term
    epsilons[3] -= (3*deltas_Q[9] + 2*deltas_Q[10] - 2*deltas_Q[11] - deltas_Q[12]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * alpha;

    // Fourteenth term
    epsilons[3] += 2*(deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 3) * beta;

    // Fifteenth term
    epsilons[3] -= (3*deltas_Q[14] - deltas_Q[15]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * alpha;

    // Sixteenth term
    epsilons[3] += 2*(deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 4) * beta;

    // Seventeenth term
    epsilons[3] -= deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * alpha;

    // Eighteenth term
    epsilons[3] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 5) * beta;

    // Nineteenth term
    epsilons[3] += 2*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[1];

    // Twentieth term
    epsilons[3] += gammas_Q[0] * pow(a_eval[ind], 4) * gammas[2];

    // Twenty-first term
    epsilons[3] -= 3*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[0];

    // Twenty-second term
    epsilons[3] += 2*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[1];

    // Twenty-third term
    epsilons[3] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[2];

    // Twenty-fourth term
    epsilons[3] += 3*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[0];

    // Twenty-fifth term
    epsilons[3] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[1];

    // Twenty-sixth term
    epsilons[3] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[2];

    // Twenty-seventh term
    epsilons[3] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[0];

    // Twenty-eighth term
    epsilons[3] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[1];

    // Twenty-ninth term
    epsilons[3] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[2];

    // Thirtieth term
    epsilons[3] += gammas_Q[1] * pow(a_eval[ind], 3) * pow(alpha, 2);

    // Thirty-first term
    epsilons[3] -= 2*gammas_Q[2] * pow(a_eval[ind], 3) * pow(alpha, 2);

    // Thirty-second term
    epsilons[3] -= (2*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * nu * pow(alpha, 2);

    // Thirty-third term
    epsilons[3] += alpha_Q * pow(a_eval[ind], 2) * deltas[1];

    // Thirty-fourth term
    epsilons[3] += beta_Q * a_eval[ind] * sigma * deltas[0];

    // Thirty-fifth term
    epsilons[3] -= beta_Q * a_eval[ind] * nu * deltas[1];

    // Thirty-sixth term
    epsilons[3] += 5*beta_Q * a_eval[ind] * alpha * gammas[0];
    break;
  }
  case 5:
  {
    // First term
    epsilons[4] += epsilons_Q[4] * pow(a_eval[ind], 7) * pow(sigma, 2);

    // Second term
    epsilons[4] += (epsilons_Q[11] + epsilons_Q[12] + epsilons_Q[13]) * pow(a_eval[ind], 6) * nu * pow(sigma, 2);

    // Third term
    epsilons[4] += (epsilons_Q[22] + epsilons_Q[23] + epsilons_Q[24] + epsilons_Q[25] + epsilons_Q[26] + epsilons_Q[27]) * pow(a_eval[ind], 5) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    epsilons[4] -= (epsilons_Q[34] + epsilons_Q[35] + epsilons_Q[36] + epsilons_Q[37] + epsilons_Q[38] + epsilons_Q[39] + epsilons_Q[40] + epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    epsilons[4] -= (epsilons_Q[45] + epsilons_Q[46] + epsilons_Q[47] + epsilons_Q[48] + epsilons_Q[49] + epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    epsilons[4] -= (epsilons_Q[52] + epsilons_Q[53] + epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 5) * pow(sigma, 2);

    // Seventh term
    epsilons[4] -= epsilons_Q[55] * a_eval[ind] * pow(nu, 6) * pow(sigma, 2);

    // Eighth term
    epsilons[4] += 2*deltas_Q[3] * pow(a_eval[ind], 5) * sigma * alpha;

    // Ninth term
    epsilons[4] += 2*(deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * alpha;

    // Tenth term
    epsilons[4] -= 2*(deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * alpha;

    // Eleventh term
    epsilons[4] -= 2*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * alpha;

    // Twelfth term
    epsilons[4] -= 2*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * alpha;

    // Thirteenth term
    epsilons[4] += 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[0];

    // Fourteenth term
    epsilons[4] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[0];

    // Fifteenth term
    epsilons[4] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[0];

    // Sixteenth term
    epsilons[4] += gammas_Q[2] * pow(a_eval[ind], 3) * pow(alpha, 2);

    // Seventeenth term
    epsilons[4] -= gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(alpha, 2);

    // Eighteenth term
    epsilons[4] -= gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(alpha, 2);

    // Nineteenth term
    epsilons[4] -= 2*beta_Q * a_eval[ind] * sigma * deltas[0];

    // Twentieth term
    epsilons[4] -= 2*beta_Q * a_eval[ind] * alpha * gammas[0];
    break;
  }
  case 6:
  {
    // First term
    epsilons[5] += epsilons_Q[5] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[5] += (10*epsilons_Q[14] + 4*epsilons_Q[15] + epsilons_Q[16] + 4*epsilons_Q[17] + epsilons_Q[18] + epsilons_Q[19] + 4*epsilons_Q[22] + epsilons_Q[23] + epsilons_Q[24]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[5] -= (18*epsilons_Q[28] + 6*epsilons_Q[29] + 12*epsilons_Q[30] + 10*epsilons_Q[31] + 4*epsilons_Q[32] + epsilons_Q[33] + 6*epsilons_Q[34] + 12*epsilons_Q[35] + 10*epsilons_Q[36] + 4*epsilons_Q[37] + epsilons_Q[38] + 4*epsilons_Q[39] + epsilons_Q[40] + epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[5] -= (21*epsilons_Q[42] + 18*epsilons_Q[43] + 12*epsilons_Q[44] + 21*epsilons_Q[45] + 18*epsilons_Q[46] + 12*epsilons_Q[47] + 6*epsilons_Q[48] + 12*epsilons_Q[49] + 10*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[5] -= (21*epsilons_Q[51] + 21*epsilons_Q[52] + 21*epsilons_Q[53] + 18*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[5] -= 21*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[5] += deltas_Q[1] * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[5] += (2*deltas_Q[4] + deltas_Q[5] + deltas_Q[6]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[5] += (3*deltas_Q[4] + 2*deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + 2*deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[5] -= (4*deltas_Q[9] + 6*deltas_Q[10] + 4*deltas_Q[11] + 2*deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[5] -= (3*deltas_Q[9] + 4*deltas_Q[10] + 4*deltas_Q[11] + 3*deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[5] -= (11*deltas_Q[14] + 8*deltas_Q[15] + 6*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[5] -= (5*deltas_Q[14] + 5*deltas_Q[15] + 4*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[5] -= 11*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[5] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[5] -= 10*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[5] -= 4*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[5] -= (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[5] -= (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[2];

    // Twentieth term
    epsilons[5] -= 10*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twenty-first term
    epsilons[5] -= 4*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-second term
    epsilons[5] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-third term
    epsilons[5] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[2];

    // Twenty-fourth term
    epsilons[5] += 10*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-fifth term
    epsilons[5] += 4*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-sixth term
    epsilons[5] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-seventh term
    epsilons[5] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[2];

    // Twenty-eighth term
    epsilons[5] += 10*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Twenty-ninth term
    epsilons[5] += 4*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Thirtieth term
    epsilons[5] -= gammas_Q[1] * pow(a_eval[ind], 3) * alpha * beta;

    // Thirty-first term
    epsilons[5] -= gammas_Q[3] * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-second term
    epsilons[5] += (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-third term
    epsilons[5] -= gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-fourth term
    epsilons[5] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-fifth term
    epsilons[5] += 4*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-sixth term
    epsilons[5] += alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-seventh term
    epsilons[5] += alpha_Q * pow(a_eval[ind], 2) * deltas[6];

    // Thirty-eigth term
    epsilons[5] -= 4*beta_Q * a_eval[ind] * nu * deltas[4];

    // Thirty-ninth term
    epsilons[5] -= beta_Q * a_eval[ind] * nu * deltas[5];

    // Fourtieth term
    epsilons[5] -= beta_Q * a_eval[ind] * nu * deltas[6];

    // Forty-first term
    epsilons[5] += beta_Q * a_eval[ind] * alpha * gammas[1];

    // Forty-second term
    epsilons[5] += beta_Q * a_eval[ind] * alpha * gammas[2];
    break;
  }
  case 7:
  {
    // First term
    epsilons[6] += epsilons_Q[6] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[6] -= (30*epsilons_Q[14] + 13*epsilons_Q[15] + 2*epsilons_Q[16] + 14*epsilons_Q[17] + 4*epsilons_Q[18] - epsilons_Q[20] + 14*epsilons_Q[22] + 4*epsilons_Q[23] - epsilons_Q[25]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[6] += (57*epsilons_Q[28] + 18*epsilons_Q[29] + 38*epsilons_Q[30] + 30*epsilons_Q[31] + 13*epsilons_Q[32] + 2*epsilons_Q[33] + 18*epsilons_Q[34] + 38*epsilons_Q[35] + 30*epsilons_Q[36] + 13*epsilons_Q[37] + 2*epsilons_Q[38] + 14*epsilons_Q[39] + 4*epsilons_Q[40]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[6] += (66*epsilons_Q[42] + 57*epsilons_Q[43] + 36*epsilons_Q[44] + 66*epsilons_Q[45] + 57*epsilons_Q[46] + 36*epsilons_Q[47] + 18*epsilons_Q[48] + 38*epsilons_Q[49] + 30*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[6] += (63*epsilons_Q[51] + 63*epsilons_Q[52] + 66*epsilons_Q[53] + 57*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[6] += 63*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[6] += 2*deltas_Q[2] * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[6] -= (4*deltas_Q[4] + 4*deltas_Q[5] + 2*deltas_Q[6] - deltas_Q[7]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[6] -= 2*(deltas_Q[4] - deltas_Q[6] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[6] += (10*deltas_Q[9] + 10*deltas_Q[10] + 2*deltas_Q[11] + 3*deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[6] += 2*(deltas_Q[10] + 2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[6] += (20*deltas_Q[14] + 9*deltas_Q[15] + 8*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[6] += 2*(deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[6] += 17*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[6] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[6] += 30*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[6] += 14*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[6] += 4*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[6] += (2*gammas_Q[1] + 3*gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[2];

    // Twentieth term
    epsilons[6] += 30*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twenty-first term
    epsilons[6] += 14*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-second term
    epsilons[6] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[0];

    // Twenty-third term
    epsilons[6] -= 8*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-fourth term
    epsilons[6] -= (4*gammas_Q[3] + 5*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[2];

    // Twenty-fifth term
    epsilons[6] -= 30*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-sixth term
    epsilons[6] -= 14*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-seventh term
    epsilons[6] -= 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[0];

    // Twenty-eighth term
    epsilons[6] -= 12*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-ninth term
    epsilons[6] -= 7*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[2];

    // Thirtieth term
    epsilons[6] -= 30*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Thirty-first term
    epsilons[6] -= 14*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Thirty-second term
    epsilons[6] += 6*gammas_Q[1] * pow(a_eval[ind], 3) * alpha * beta;

    // Thirty-third term
    epsilons[6] -= 6*(2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-fourth term
    epsilons[6] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-fifth term
    epsilons[6] -= 12*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-sixth term
    epsilons[6] -= 4*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-seventh term
    epsilons[6] -= 14*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-eigth term
    epsilons[6] -= 4*alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-ninth term
    epsilons[6] += alpha_Q * pow(a_eval[ind], 2) * deltas[7];

    // Fourtieth term
    epsilons[6] += 2*beta_Q * a_eval[ind] * beta * gammas[0];

    // Forty-first term
    epsilons[6] -= 4*beta_Q * a_eval[ind] * alpha * gammas[1];

    // Forty-second term
    epsilons[6] -= 4*beta_Q * a_eval[ind] * alpha * gammas[2];

    // Forty-third term
    epsilons[6] += 14*beta_Q * a_eval[ind] * nu * deltas[4];

    // Forty-fourth term
    epsilons[6] += 4*beta_Q * a_eval[ind] * nu * deltas[5];

    // Forty-fifth term
    epsilons[6] -= beta_Q * a_eval[ind] * nu * deltas[7];
    break;
  }
  case 8:
  {
    // First term
    epsilons[7] += epsilons_Q[7] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[7] += (17*epsilons_Q[14] + 10*epsilons_Q[15] + 4*epsilons_Q[16] + 9*epsilons_Q[17] + 4*epsilons_Q[18] + epsilons_Q[21] + 9*epsilons_Q[22] + 4*epsilons_Q[23] + epsilons_Q[26]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[7] -= (35*epsilons_Q[28] + 13*epsilons_Q[29] + 23*epsilons_Q[30] + 17*epsilons_Q[31] + 10*epsilons_Q[32] + 4*epsilons_Q[33] + 13*epsilons_Q[34] + 23*epsilons_Q[35] + 17*epsilons_Q[36] + 10*epsilons_Q[37] + 4*epsilons_Q[38] + 9*epsilons_Q[39] + 4*epsilons_Q[40]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[7] -= (41*epsilons_Q[42] + 35*epsilons_Q[43] + 24*epsilons_Q[44] + 41*epsilons_Q[45] + 35*epsilons_Q[46] + 24*epsilons_Q[47] + 13*epsilons_Q[48] + 23*epsilons_Q[49] + 17*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[7] -= (42*epsilons_Q[51] + 42*epsilons_Q[52] + 41*epsilons_Q[53] + 35*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[7] -= 42*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[7] -= (2*deltas_Q[1] + deltas_Q[2]) * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[7] += (deltas_Q[4] + 2*deltas_Q[5] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[7] -= (5*deltas_Q[4] + 4*deltas_Q[5] + 3*deltas_Q[6] + 2*deltas_Q[7] + 3*deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[7] -= (3*deltas_Q[9] + deltas_Q[10] - deltas_Q[11] + 2*deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[7] += (6*deltas_Q[9] + 7*deltas_Q[10] + 6*deltas_Q[11] + 5*deltas_Q[12] + 4*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[7] -= (3*deltas_Q[14] + deltas_Q[15] + 2*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[7] += (9*deltas_Q[14] + 8*deltas_Q[15] + 7*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[7] -= 4*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[7] += 9*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[7] -= 15*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[7] -= 9*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[7] -= (2*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[7] -= 15*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twentieth term
    epsilons[7] -= 9*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-first term
    epsilons[7] += (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[0];

    // Twenty-second term
    epsilons[7] += (4*gammas_Q[3] + 3*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-third term
    epsilons[7] += 15*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-fourth term
    epsilons[7] += 9*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-fifth term
    epsilons[7] += 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[0];

    // Twenty-sixth term
    epsilons[7] += 5*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-seventh term
    epsilons[7] += 15*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Twenty-eighth term
    epsilons[7] += 9*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Twenty-ninth term
    epsilons[7] -= (3*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * alpha * beta;

    // Thirtieth term
    epsilons[7] += (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-first term
    epsilons[7] += (6*gammas_Q[3] + 4*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-second term
    epsilons[7] += 4*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-third term
    epsilons[7] += 7*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-fourth term
    epsilons[7] += 9*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-fifth term
    epsilons[7] += 4*alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-sixth term
    epsilons[7] += alpha_Q * pow(a_eval[ind], 2) * deltas[8];

    // Thirty-seventh term
    epsilons[7] -= 9*beta_Q * a_eval[ind] * nu * deltas[4];

    // Thirty-eigth term
    epsilons[7] -= 4*beta_Q * a_eval[ind] * nu * deltas[5];

    // Thirty-ninth term
    epsilons[7] -= beta_Q * a_eval[ind] * nu * deltas[8];

    // Fourtieth term
    epsilons[7] -= beta_Q * a_eval[ind] * beta * gammas[0];
    break;
  }
  case 9:
  {
    // First term
    epsilons[8] += epsilons_Q[8] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[8] -= (30*epsilons_Q[14] + 16*epsilons_Q[15] + 6*epsilons_Q[16] + 13*epsilons_Q[17] + 4*epsilons_Q[18] - epsilons_Q[19] - epsilons_Q[20] - 2*epsilons_Q[21] + 16*epsilons_Q[22] + 6*epsilons_Q[23]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[8] += (62*epsilons_Q[28] + 19*epsilons_Q[29] + 38*epsilons_Q[30] + 26*epsilons_Q[31] + 13*epsilons_Q[32] + 4*epsilons_Q[33] + 22*epsilons_Q[34] + 42*epsilons_Q[35] + 30*epsilons_Q[36] + 16*epsilons_Q[37] + 6*epsilons_Q[38] + 13*epsilons_Q[39] + 4*epsilons_Q[40] - epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[8] += (69*epsilons_Q[42] + 57*epsilons_Q[43] + 38*epsilons_Q[44] + 74*epsilons_Q[45] + 62*epsilons_Q[46] + 42*epsilons_Q[47] + 19*epsilons_Q[48] + 38*epsilons_Q[49] + 26*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[8] += (69*epsilons_Q[51] + 74*epsilons_Q[52] + 69*epsilons_Q[53] + 57*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[8] += 69*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[8] -= 2*deltas_Q[2] * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[8] -= (4*deltas_Q[4] + 4*deltas_Q[5] - 2*deltas_Q[6] + 2*deltas_Q[7] + 4*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[8] += 2*(deltas_Q[4] - deltas_Q[6] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[8] += (6*deltas_Q[9] + 6*deltas_Q[10] + 6*deltas_Q[11] + 10*deltas_Q[12] + 8*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[8] -= 2*(deltas_Q[10] + 2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[8] += (12*deltas_Q[14] + 16*deltas_Q[15] + 14*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[8] -= 2*(deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[8] += 22*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[8] -= 2*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[8] += 30*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[8] += 16*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[8] += 6*gammas_Q[1] * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[8] -= 3*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[2];

    // Twentieth term
    epsilons[8] += 30*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twenty-first term
    epsilons[8] += 16*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-second term
    epsilons[8] -= 6*(2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-third term
    epsilons[8] += 3*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[2];

    // Twenty-fourth term
    epsilons[8] -= 30*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-fifth term
    epsilons[8] -= 16*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-sixth term
    epsilons[8] -= 12*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-seventh term
    epsilons[8] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[2];

    // Twenty-eighth term
    epsilons[8] -= 30*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Twenty-ninth term
    epsilons[8] -= 16*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Thirtieth term
    epsilons[8] += 4*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * alpha * beta;

    // Thirty-first term
    epsilons[8] -= 7*gammas_Q[3] * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-second term
    epsilons[8] -= 8*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-third term
    epsilons[8] -= 7*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-fourth term
    epsilons[8] -= 12*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-fifth term
    epsilons[8] -= 16*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-sixth term
    epsilons[8] -= 6*alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-seventh term
    epsilons[8] += beta_Q * a_eval[ind] * sigma * deltas[1];

    // Thirty-eigth term
    epsilons[8] += 16*beta_Q * a_eval[ind] * nu * deltas[4];

    // Thirty-ninth term
    epsilons[8] += 6*beta_Q * a_eval[ind] * nu * deltas[5];

    // Fourtieth term
    epsilons[8] += 4*beta_Q * a_eval[ind] * alpha * gammas[1];

    // Forty-first term
    epsilons[8] += 5*beta_Q * a_eval[ind] * alpha * gammas[2];
    break;
  }
  case 10:
  {
    // First term
    epsilons[9] += epsilons_Q[9] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[9] += (45*epsilons_Q[14] + 22*epsilons_Q[15] + 7*epsilons_Q[16] + 21*epsilons_Q[17] + 7*epsilons_Q[18] + epsilons_Q[19] - epsilons_Q[21] + 22*epsilons_Q[22] + 7*epsilons_Q[23] + epsilons_Q[27]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[9] -= (89*epsilons_Q[28] + 29*epsilons_Q[29] + 58*epsilons_Q[30] + 43*epsilons_Q[31] + 21*epsilons_Q[32] + 7*epsilons_Q[33] + 29*epsilons_Q[34] + 59*epsilons_Q[35] + 45*epsilons_Q[36] + 22*epsilons_Q[37] + 7*epsilons_Q[38] + 21*epsilons_Q[39] + 7*epsilons_Q[40] + epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[9] -= (102*epsilons_Q[42] + 87*epsilons_Q[43] + 58*epsilons_Q[44] + 103*epsilons_Q[45] + 89*epsilons_Q[46] + 59*epsilons_Q[47] + 29*epsilons_Q[48] + 58*epsilons_Q[49] + 43*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[9] -= (102*epsilons_Q[51] + 103*epsilons_Q[52] + 102*epsilons_Q[53] + 87*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[9] -= 102*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[9] -= deltas_Q[2] * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[9] += (8*deltas_Q[4] + 7*deltas_Q[5] - 2*deltas_Q[6] + 2*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[9] += (deltas_Q[4] - deltas_Q[6] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[9] -= (13*deltas_Q[9] + 15*deltas_Q[10] + 5*deltas_Q[11] + 10*deltas_Q[12] + 7*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[9] -= (deltas_Q[10] + 2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[9] -= (29*deltas_Q[14] + 21*deltas_Q[15] + 17*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[9] -= (deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[9] -= 31*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[9] -= deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[9] -= 45*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[9] -= 22*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[9] -= (7*gammas_Q[1] + 4*gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[9] -= 45*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twentieth term
    epsilons[9] -= 22*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-first term
    epsilons[9] += (3*gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[0];

    // Twenty-second term
    epsilons[9] += (14*gammas_Q[3] + 11*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-third term
    epsilons[9] += 45*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-fourth term
    epsilons[9] += 22*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-fifth term
    epsilons[9] += 7*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[0];

    // Twenty-sixth term
    epsilons[9] += 18*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-seventh term
    epsilons[9] += 45*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Twenty-eighth term
    epsilons[9] += 22*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Twenty-ninth term
    epsilons[9] -= 2*(4*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * alpha * beta;

    // Thirtieth term
    epsilons[9] += 2*(4*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-first term
    epsilons[9] += (16*gammas_Q[3] + 10*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-second term
    epsilons[9] += 12*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-third term
    epsilons[9] += 18*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-fourth term
    epsilons[9] += 22*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-fifth term
    epsilons[9] += 7*alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-sixth term
    epsilons[9] += beta_Q * a_eval[ind] * sigma * deltas[2];

    // Thirty-seventh term
    epsilons[9] -= 22*beta_Q * a_eval[ind] * nu * deltas[4];

    // Thirty-eigth term
    epsilons[9] -= 7*beta_Q * a_eval[ind] * nu * deltas[5];

    // Thirty-ninth term
    epsilons[9] -= 2*beta_Q * a_eval[ind] * beta * gammas[0];

    // Fourtieth term
    epsilons[9] += 2*beta_Q * a_eval[ind] * alpha * gammas[1];
    break;
  }
  case 11:
  {
    // First term
    epsilons[10] += epsilons_Q[10] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[10] -= (9*epsilons_Q[14] + 4*epsilons_Q[15] + epsilons_Q[16] + 4*epsilons_Q[17] + epsilons_Q[18] - epsilons_Q[20] - epsilons_Q[21] + 4*epsilons_Q[22] + epsilons_Q[23]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[10] += (17*epsilons_Q[28] + 5*epsilons_Q[29] + 11*epsilons_Q[30] + 8*epsilons_Q[31] + 3*epsilons_Q[32] + 5*epsilons_Q[34] + 11*epsilons_Q[35] + 9*epsilons_Q[36] + 4*epsilons_Q[37]+epsilons_Q[38] + 4*epsilons_Q[39] + epsilons_Q[40]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[10] += (19*epsilons_Q[42] + 16*epsilons_Q[43] + 10*epsilons_Q[44] + 19*epsilons_Q[45] + 17*epsilons_Q[46] + 11*epsilons_Q[47] + 5*epsilons_Q[48] + 11*epsilons_Q[49] + 8*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[10] += (18*epsilons_Q[51] + 19*epsilons_Q[52] + 19*epsilons_Q[53] + 16*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[10] += 18*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[10] += (deltas_Q[2] + deltas_Q[3]) * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[10] -= (2*deltas_Q[4] + deltas_Q[5] - 2*deltas_Q[6]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[10] -= (deltas_Q[4] - deltas_Q[6] - deltas_Q[7]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[10] += (deltas_Q[9] + 3*deltas_Q[10] + deltas_Q[11] + 2*deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[10] += (deltas_Q[10] + deltas_Q[11] - deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[10] += (5*deltas_Q[14] + 5*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[10] += (deltas_Q[14] + deltas_Q[15]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[10] += 5*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[10] += 9*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Sixteenth term
    epsilons[10] += 4*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Seventeenth term
    epsilons[10] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Eighteenth term
    epsilons[10] += 9*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Nineteenth term
    epsilons[10] += 4*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twentieth term
    epsilons[10] -= 3*gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[0];

    // Twenty-first term
    epsilons[10] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-second term
    epsilons[10] -= 9*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-third term
    epsilons[10] -= 4*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-fourth term
    epsilons[10] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[0];

    // Twenty-fifth term
    epsilons[10] -= 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-sixth term
    epsilons[10] -= 9*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Twenty-seventh term
    epsilons[10] -= 4*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Twenty-eighth term
    epsilons[10] += (2*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * alpha * beta;

    // Twenty-ninth term
    epsilons[10] -= 2*gammas_Q[3] * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirtieth term
    epsilons[10] -= (4*gammas_Q[3] + 3*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-first term
    epsilons[10] -= 2*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-second term
    epsilons[10] -= 5*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-third term
    epsilons[10] -= 4*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-fourth term
    epsilons[10] -= alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-fifth term
    epsilons[10] += beta_Q * a_eval[ind] * sigma * deltas[3];

    // Thirty-sixth term
    epsilons[10] += 4*beta_Q * a_eval[ind] * nu * deltas[4];

    // Thirty-seventh term
    epsilons[10] += beta_Q * a_eval[ind] * nu * deltas[5];

    // Thirty-eigth term
    epsilons[10] -= beta_Q * a_eval[ind] * beta * gammas[0];

    // Thirty-ninth term
    epsilons[10] -= beta_Q * a_eval[ind] * alpha * gammas[1];
    break;
  }
  case 12:
  {
    // First term
    epsilons[11] += epsilons_Q[11] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[11] += (18*epsilons_Q[14] + 9*epsilons_Q[15] + 3*epsilons_Q[16] + 9*epsilons_Q[17] + 3*epsilons_Q[18] + 12*epsilons_Q[22] + 5*epsilons_Q[23] + epsilons_Q[24] + epsilons_Q[25] + 2*epsilons_Q[26] + epsilons_Q[27]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[11] -= (36*epsilons_Q[28] + 12*epsilons_Q[29] + 24*epsilons_Q[30] + 18*epsilons_Q[31] + 9*epsilons_Q[32] + 3*epsilons_Q[33] + 15*epsilons_Q[34] + 28*epsilons_Q[35] + 22*epsilons_Q[36] + 12*epsilons_Q[37] + 5*epsilons_Q[38] + 13*epsilons_Q[39] + 6*epsilons_Q[40] + 2*epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[11] -= (42*epsilons_Q[42] + 36*epsilons_Q[43] + 24*epsilons_Q[44] + 47*epsilons_Q[45] + 41*epsilons_Q[46] + 28*epsilons_Q[47] + 16*epsilons_Q[48] + 29*epsilons_Q[49] + 23*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[11] -= (42*epsilons_Q[51] + 47*epsilons_Q[52] + 48*epsilons_Q[53] + 42*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[11] -= 48*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[11] -= deltas_Q[3] * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[11] += (3*deltas_Q[4] + 3*deltas_Q[5] + 2*deltas_Q[7] + 2*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[11] -= (deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[11] -= (6*deltas_Q[9] + 6*deltas_Q[10] + 4*deltas_Q[11] + 7*deltas_Q[12] + 7*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[11] += (deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[11] -= 12*(deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[11] += (deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[11] -= 20*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[11] += deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[11] -= 18*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[11] -= 9*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[11] -= (3*gammas_Q[1] - gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[11] += 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[2];

    // Twentieth term
    epsilons[11] -= 18*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twenty-first term
    epsilons[11] -= 9*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-second term
    epsilons[11] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[0];

    // Twenty-third term
    epsilons[11] += 2*(3*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-fourth term
    epsilons[11] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[2];

    // Twenty-fifth term
    epsilons[11] += 18*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-sixth term
    epsilons[11] += 9*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-seventh term
    epsilons[11] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[0];

    // Twenty-eighth term
    epsilons[11] += 5*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-ninth term
    epsilons[11] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[2];

    // Thirtieth term
    epsilons[11] += 18*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Thirty-first term
    epsilons[11] += 9*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Thirty-second term
    epsilons[11] -= (3*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * alpha * beta;

    // Thirty-third term
    epsilons[11] += (3*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-fourth term
    epsilons[11] += (6*gammas_Q[3] + 4*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-fifth term
    epsilons[11] += gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-sixth term
    epsilons[11] += 7*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-seventh term
    epsilons[11] += 9*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-eigth term
    epsilons[11] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-ninth term
    epsilons[11] -= 2*beta_Q * a_eval[ind] * sigma * deltas[1];

    // Fourtieth term
    epsilons[11] -= 9*beta_Q * a_eval[ind] * nu * deltas[4];

    // Forty-first term
    epsilons[11] -= 3*beta_Q * a_eval[ind] * nu * deltas[5];

    // Forty-second term
    epsilons[11] -= 2*beta_Q * a_eval[ind] * beta * gammas[0];

    // Forty-third term
    epsilons[11] -= beta_Q * a_eval[ind] * alpha * gammas[1];

    // Forty-fourth term
    epsilons[11] -= 2*beta_Q * a_eval[ind] * alpha * gammas[2];
    break;
  }
  case 13:
  {
    // First term
    epsilons[12] += epsilons_Q[12] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[12] -= (24*epsilons_Q[14] + 12*epsilons_Q[15] + 4*epsilons_Q[16] + 12*epsilons_Q[17] + 4*epsilons_Q[18] + 13*epsilons_Q[22] + 4*epsilons_Q[23] - epsilons_Q[24] + epsilons_Q[26] + 2*epsilons_Q[27]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[12] += (48*epsilons_Q[28] + 16*epsilons_Q[29] + 32*epsilons_Q[30] + 24*epsilons_Q[31] + 12*epsilons_Q[32] + 4*epsilons_Q[33] + 16*epsilons_Q[34] + 33*epsilons_Q[35] + 26*epsilons_Q[36] + 13*epsilons_Q[37] + 4*epsilons_Q[38] + 16*epsilons_Q[39] + 7*epsilons_Q[40] + 2*epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[12] += (56*epsilons_Q[42] + 48*epsilons_Q[43] + 32*epsilons_Q[44] + 57*epsilons_Q[45] + 50*epsilons_Q[46] + 33*epsilons_Q[47] + 19*epsilons_Q[48] + 36*epsilons_Q[49] + 29*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[12] += (56*epsilons_Q[51] + 57*epsilons_Q[52] + 60*epsilons_Q[53] + 53*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[12] += 60*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[12] -= 2*deltas_Q[3] * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[12] -= 2*(2*deltas_Q[4] + 2*deltas_Q[5] + deltas_Q[7]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[12] -= 2*(deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[12] += (8*deltas_Q[9] + 8*deltas_Q[10] + 6*deltas_Q[12] + 8*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[12] += 2*(deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[12] += (16*deltas_Q[14] + 10*deltas_Q[15] + 12*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[12] += 2*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[12] += 22*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[12] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[12] += 24*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[12] += 12*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[12] += 2*(2*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[12] += 24*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twentieth term
    epsilons[12] += 12*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-first term
    epsilons[12] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[0];

    // Twenty-second term
    epsilons[12] -= (8*gammas_Q[3] + 6*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-third term
    epsilons[12] -= 24*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-fourth term
    epsilons[12] -= 12*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-fifth term
    epsilons[12] -= 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[0];

    // Twenty-sixth term
    epsilons[12] -= 10*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-seventh term
    epsilons[12] -= 24*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Twenty-eighth term
    epsilons[12] -= 12*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Twenty-ninth term
    epsilons[12] += 2*(2*gammas_Q[1] - gammas_Q[2]) * pow(a_eval[ind], 3) * alpha * beta;

    // Thirtieth term
    epsilons[12] -= (4*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-first term
    epsilons[12] -= 2*(4*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-second term
    epsilons[12] -= 2*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-third term
    epsilons[12] -= 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-fourth term
    epsilons[12] -= 12*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-fifth term
    epsilons[12] -= 4*alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-sixth term
    epsilons[12] -= 2*beta_Q * a_eval[ind] * sigma * deltas[2];

    // Thirty-seventh term
    epsilons[12] += 12*beta_Q * a_eval[ind] * nu * deltas[4];

    // Thirty-eigth term
    epsilons[12] += 4*beta_Q * a_eval[ind] * nu * deltas[5];

    // Thirty-ninth term
    epsilons[12] += 6*beta_Q * a_eval[ind] * beta * gammas[0];

    // Fourtieth term
    epsilons[12] -= 2*beta_Q * a_eval[ind] * alpha * gammas[1];
    break;
  }
  case 14:
  {
    // First term
    epsilons[13] += epsilons_Q[13] * pow(a_eval[ind], 6) * pow(sigma, 3);

    // Second term
    epsilons[13] += (6*epsilons_Q[14] + 3*epsilons_Q[15] + epsilons_Q[16] + 3*epsilons_Q[17] + epsilons_Q[18] + 3*epsilons_Q[22] + epsilons_Q[23] + epsilons_Q[25] + epsilons_Q[26] + 3*epsilons_Q[27]) * pow(a_eval[ind], 5) * nu * pow(sigma, 3);

    // Third term
    epsilons[13] -= (12*epsilons_Q[28] + 4*epsilons_Q[29] + 8*epsilons_Q[30] + 6*epsilons_Q[31] + 3*epsilons_Q[32] + epsilons_Q[33] + 4*epsilons_Q[34] + 8*epsilons_Q[35] + 7*epsilons_Q[36] + 4*epsilons_Q[37] + 2*epsilons_Q[38] + 6*epsilons_Q[39] + 4*epsilons_Q[40] + 3*epsilons_Q[41]) * pow(a_eval[ind], 4) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    epsilons[13] -= (14*epsilons_Q[42] + 12*epsilons_Q[43] + 8*epsilons_Q[44] + 14*epsilons_Q[45] + 13*epsilons_Q[46] + 9*epsilons_Q[47] + 7*epsilons_Q[48] + 11*epsilons_Q[49] + 10*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    epsilons[13] -= (14*epsilons_Q[51] + 15*epsilons_Q[52] + 17*epsilons_Q[53] + 16*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 4) * pow(sigma, 3);

    // Sixth term
    epsilons[13] -= 18*epsilons_Q[55] * a_eval[ind] * pow(nu, 5) * pow(sigma, 3);

    // Seventh term
    epsilons[13] += deltas_Q[3] * pow(a_eval[ind], 5) * sigma * beta;

    // Eighth term
    epsilons[13] += (deltas_Q[4] + deltas_Q[5] + 2*deltas_Q[7]) * pow(a_eval[ind], 4) * pow(sigma, 2) * alpha;

    // Ninth term
    epsilons[13] += (deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * sigma * beta;

    // Tenth term
    epsilons[13] -= (2*deltas_Q[9] + 2*deltas_Q[10] + 3*deltas_Q[12] + 5*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * alpha;

    // Eleventh term
    epsilons[13] -= (deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * sigma * beta;

    // Twelfth term
    epsilons[13] -= (4*deltas_Q[14] + 4*deltas_Q[15] + 6*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * alpha;

    // Thirteenth term
    epsilons[13] -= (deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * sigma * beta;

    // Fourteenth term
    epsilons[13] -= 10*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * alpha;

    // Fifteenth term
    epsilons[13] -= deltas_Q[17] * a_eval[ind] * pow(nu, 4) * sigma * beta;

    // Sixteenth term
    epsilons[13] -= 6*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[3];

    // Seventeenth term
    epsilons[13] -= 3*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[4];

    // Eighteenth term
    epsilons[13] -= (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[1];

    // Nineteenth term
    epsilons[13] -= 6*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[3];

    // Twentieth term
    epsilons[13] -= 3*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[4];

    // Twenty-first term
    epsilons[13] -= gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[0];

    // Twenty-second term
    epsilons[13] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[1];

    // Twenty-third term
    epsilons[13] += 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[3];

    // Twenty-fourth term
    epsilons[13] += 3*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[4];

    // Twenty-fifth term
    epsilons[13] -= 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[0];

    // Twenty-sixth term
    epsilons[13] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[1];

    // Twenty-seventh term
    epsilons[13] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[3];

    // Twenty-eighth term
    epsilons[13] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[4];

    // Twenty-ninth term
    epsilons[13] -= (gammas_Q[1] - gammas_Q[2]) * pow(a_eval[ind], 3) * alpha * beta;

    // Thirtieth term
    epsilons[13] += (gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(alpha, 2);

    // Thirty-first term
    epsilons[13] += 2*gammas_Q[3] * pow(a_eval[ind], 2) * nu * alpha * beta;

    // Thirty-second term
    epsilons[13] -= gammas_Q[5] * a_eval[ind] * nu * sigma * pow(alpha, 2);

    // Thirty-third term
    epsilons[13] += gammas_Q[5] * a_eval[ind] * pow(nu, 2) * alpha * beta;

    // Thirty-fourth term
    epsilons[13] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[4];

    // Thirty-fifth term
    epsilons[13] += alpha_Q * pow(a_eval[ind], 2) * deltas[5];

    // Thirty-sixth term
    epsilons[13] -= 2*beta_Q * a_eval[ind] * sigma * deltas[3];

    // Thirty-seventh term
    epsilons[13] -= 3*beta_Q * a_eval[ind] * nu * deltas[4];

    // Thirty-eigth term
    epsilons[13] -= beta_Q * a_eval[ind] * nu * deltas[5];

    // Thirty-ninth term
    epsilons[13] -= 2*beta_Q * a_eval[ind] * beta * gammas[0];

    // Fourtieth term
    epsilons[13] += beta_Q * a_eval[ind] * alpha * gammas[1];
    break;
  }
  case 15:
  {
    // First term
    epsilons[14] += epsilons_Q[14] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[14] -= (9*epsilons_Q[28] + 4*epsilons_Q[30] + epsilons_Q[31] + 4*epsilons_Q[35] + epsilons_Q[36] - 3*epsilons_Q[39] - 2*epsilons_Q[40] - epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[14] -= (12*epsilons_Q[42] + 9*epsilons_Q[43] + 5*epsilons_Q[44] + 12*epsilons_Q[45] + 9*epsilons_Q[46] + 5*epsilons_Q[47] - 3*epsilons_Q[48] - 5*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[14] -= (14*epsilons_Q[51] + 14*epsilons_Q[52] + 7*epsilons_Q[53] + 2*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[14] -= 7*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[14] += (deltas_Q[4] + deltas_Q[6] - deltas_Q[7] - 2*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[14] -= (2*deltas_Q[10] + deltas_Q[11] - 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[14] -= 2*(deltas_Q[9] + 2*deltas_Q[10] - deltas_Q[12] - deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[14] -= (6*deltas_Q[14] + 5*deltas_Q[15] - deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[14] -= (7*deltas_Q[14] + deltas_Q[15]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[14] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[14] -= 2*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[14] -= 8*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[14] += (2*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fifteenth term
    epsilons[14] -= 8*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Sixteenth term
    epsilons[14] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Seventeenth term
    epsilons[14] -= (gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Eighteenth term
    epsilons[14] -= (4*gammas_Q[3] + 3*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Nineteenth term
    epsilons[14] += 8*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Twentieth term
    epsilons[14] += 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twenty-first term
    epsilons[14] -= 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-second term
    epsilons[14] -= 5*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-third term
    epsilons[14] -= 5*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-fourth term
    epsilons[14] += 8*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-fifth term
    epsilons[14] += 2*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-sixth term
    epsilons[14] -= (4*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twenty-seventh term
    epsilons[14] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-eighth term
    epsilons[14] -= 6*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Twenty-ninth term
    epsilons[14] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Thirtieth term
    epsilons[14] -= 4*alpha_Q * pow(a_eval[ind], 2) * deltas[10];

    // Thirty-first term
    epsilons[14] -= alpha_Q * pow(a_eval[ind], 2) * deltas[11];

    // Thirty-second term
    epsilons[14] += 4*beta_Q * a_eval[ind] * nu * deltas[10];

    // Thirty-third term
    epsilons[14] += beta_Q * a_eval[ind] * nu * deltas[11];

    // Thirty-fourth term
    epsilons[14] += 4*beta_Q * a_eval[ind] * beta * gammas[1];

    // Thirty-fifth term
    epsilons[14] += 3*beta_Q * a_eval[ind] * beta * gammas[2];
    break;
  }
  case 16:
  {
    // First term
    epsilons[15] += epsilons_Q[15] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[15] += (10*epsilons_Q[28] - 2*epsilons_Q[29] + 5*epsilons_Q[30] - epsilons_Q[32] - 2*epsilons_Q[34] + 5*epsilons_Q[35] - epsilons_Q[37] - 12*epsilons_Q[39] - 8*epsilons_Q[40] - 4*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[15] += (13*epsilons_Q[42] + 10*epsilons_Q[43] + 5*epsilons_Q[44] + 13*epsilons_Q[45] + 10*epsilons_Q[46] + 5*epsilons_Q[47] - 14*epsilons_Q[48] - 11*epsilons_Q[49] - 24*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[15] += (14*epsilons_Q[51] + 14*epsilons_Q[52] - 7*epsilons_Q[53] - 18*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[15] -= 14*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[15] -= (2*deltas_Q[4] + 3*deltas_Q[6] - 4*deltas_Q[7] - 8*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[15] += (5*deltas_Q[10] - deltas_Q[12] - 8*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[15] += (5*deltas_Q[9] + 9*deltas_Q[10] - 7*deltas_Q[11] - 10*deltas_Q[12] - 8*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[15] += (15*deltas_Q[14] + 8*deltas_Q[15] - 11*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[15] += (16*deltas_Q[14] - 7*deltas_Q[15] - 7*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[15] -= 24*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[15] -= 4*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[15] += 12*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[15] -= 5*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fifteenth term
    epsilons[15] -= gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[4];

    // Sixteenth term
    epsilons[15] += 12*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Seventeenth term
    epsilons[15] += 8*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Eighteenth term
    epsilons[15] += (3*gammas_Q[3] + 7*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Nineteenth term
    epsilons[15] += 10*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Twentieth term
    epsilons[15] += gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[4];

    // Twenty-first term
    epsilons[15] -= 12*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Twenty-second term
    epsilons[15] -= 32*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twenty-third term
    epsilons[15] += 16*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-fourth term
    epsilons[15] += 17*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-fifth term
    epsilons[15] += 15*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-sixth term
    epsilons[15] += gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[4];

    // Twenty-seventh term
    epsilons[15] -= 12*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-eighth term
    epsilons[15] -= 8*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-ninth term
    epsilons[15] += (11*gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Thirtieth term
    epsilons[15] += 8*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Thirty-first term
    epsilons[15] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Thirty-second term
    epsilons[15] += 15*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Thirty-third term
    epsilons[15] += 8*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Thirty-fourth term
    epsilons[15] -= 2*alpha_Q * pow(a_eval[ind], 2) * deltas[9];

    // Thirty-fifth term
    epsilons[15] += 5*alpha_Q * pow(a_eval[ind], 2) * deltas[10];

    // Thirty-sixth term
    epsilons[15] -= alpha_Q * pow(a_eval[ind], 2) * deltas[12];

    // Thirty-seventh term
    epsilons[15] += 2*beta_Q * a_eval[ind] * nu * deltas[9];

    // Thirty-eigth term
    epsilons[15] -= 5*beta_Q * a_eval[ind] * nu * deltas[10];

    // Thirty-ninth term
    epsilons[15] += beta_Q * a_eval[ind] * nu * deltas[12];

    // Fourtieth term
    epsilons[15] -= 16*beta_Q * a_eval[ind] * beta * gammas[1];

    // Forty-first term
    epsilons[15] -= 11*beta_Q * a_eval[ind] * beta * gammas[2];

    // Forty-second term
    epsilons[15] += 5*beta_Q * a_eval[ind] * alpha * gammas[3];

    // Forty-third term
    epsilons[15] += 2*beta_Q * a_eval[ind] * alpha * gammas[4];
    break;
  }
  case 17:
  {
    // First term
    epsilons[16] += epsilons_Q[16] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[16] -= (6*epsilons_Q[28] + 3*epsilons_Q[30] + epsilons_Q[33] + 3*epsilons_Q[35] + epsilons_Q[38] - 9*epsilons_Q[39] - 6*epsilons_Q[40] - 3*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[16] -= (9*epsilons_Q[42] + 6*epsilons_Q[43] + 5*epsilons_Q[44] + 9*epsilons_Q[45] + 6*epsilons_Q[46] + 5*epsilons_Q[47] - 9*epsilons_Q[48] - 9*epsilons_Q[49] - 18*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[16] -= (14*epsilons_Q[51] + 14*epsilons_Q[52] - 6*epsilons_Q[53] - 15*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[16] += 7*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[16] += (deltas_Q[4] + deltas_Q[5] + 3*deltas_Q[6] - 3*deltas_Q[7] - 6*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[16] -= (3*deltas_Q[10] - 5*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[16] -= (5*deltas_Q[9] + 5*deltas_Q[10] - 9*deltas_Q[11] - 8*deltas_Q[12] - 5*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[16] -= (9*deltas_Q[14] + 6*deltas_Q[15] - 7*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[16] -= (10*deltas_Q[14] - 10*deltas_Q[15] - 7*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[16] += 16*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[16] += 5*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[16] -= 6*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[16] += 3*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fifteenth term
    epsilons[16] -= 6*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Sixteenth term
    epsilons[16] -= (gammas_Q[3] + 7*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Seventeenth term
    epsilons[16] -= 3*(gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Eighteenth term
    epsilons[16] -= 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Nineteenth term
    epsilons[16] += 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Twentieth term
    epsilons[16] += 24*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twenty-first term
    epsilons[16] -= 15*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-second term
    epsilons[16] -= 15*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-third term
    epsilons[16] -= 9*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-fourth term
    epsilons[16] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-fifth term
    epsilons[16] -= (gammas_Q[1] - 6*gammas_Q[2]) * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-sixth term
    epsilons[16] -= 8*gammas_Q[3] * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twenty-seventh term
    epsilons[16] += (2*gammas_Q[3] - 5*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-eighth term
    epsilons[16] += 3*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-ninth term
    epsilons[16] -= 8*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Thirtieth term
    epsilons[16] -= 4*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Thirty-first term
    epsilons[16] -= 3*alpha_Q * pow(a_eval[ind], 2) * deltas[10];

    // Thirty-second term
    epsilons[16] -= alpha_Q * pow(a_eval[ind], 2) * deltas[13];

    // Thirty-third term
    epsilons[16] += 3*beta_Q * a_eval[ind] * nu * deltas[10];

    // Thirty-fourth term
    epsilons[16] += beta_Q * a_eval[ind] * nu * deltas[13];

    // Thirty-fifth term
    epsilons[16] += 13*beta_Q * a_eval[ind] * beta * gammas[1];

    // Thirty-sixth term
    epsilons[16] += 9*beta_Q * a_eval[ind] * beta * gammas[2];

    // Thirty-seventh term
    epsilons[16] -= 3*beta_Q * a_eval[ind] * alpha * gammas[3];
    break;
  }
  case 18:
  {
    // First term
    epsilons[17] += epsilons_Q[17] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[17] += (6*epsilons_Q[28] - epsilons_Q[29] - 3*epsilons_Q[31] - epsilons_Q[32] + 3*epsilons_Q[35] + 8*epsilons_Q[39] + 6*epsilons_Q[40] + 3*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[17] += (4*epsilons_Q[42] + epsilons_Q[43] + 9*epsilons_Q[45] + 6*epsilons_Q[46] + 3*epsilons_Q[47] + 8*epsilons_Q[48] + 12*epsilons_Q[49] + 15*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[17] += (4*epsilons_Q[51] + 9*epsilons_Q[52] + 19*epsilons_Q[53] + 22*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[17] += 25*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[17] += (deltas_Q[5] - 3*deltas_Q[7] - 6*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[17] -= (2*deltas_Q[9] + deltas_Q[10] - 6*deltas_Q[11] - 2*deltas_Q[12] - 6*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[17] -= (deltas_Q[9] - deltas_Q[10] - 15*deltas_Q[11] - 9*deltas_Q[12] - 5*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[17] -= (5*deltas_Q[14] - 10*deltas_Q[15] - 17*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[17] += (deltas_Q[14] + 18*deltas_Q[15] + 13*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[17] += 35*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[17] += 16*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[17] += 6*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[17] -= 3*(gammas_Q[1] - 2 * gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fifteenth term
    epsilons[17] += 3*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[4];

    // Sixteenth term
    epsilons[17] += 6*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Seventeenth term
    epsilons[17] -= 6*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Eighteenth term
    epsilons[17] -= 3*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Nineteenth term
    epsilons[17] += 3*(2*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Twentieth term
    epsilons[17] -= 3*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[4];

    // Twenty-first term
    epsilons[17] -= 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Twenty-second term
    epsilons[17] += 24*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twenty-third term
    epsilons[17] -= 12*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-fourth term
    epsilons[17] -= 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-fifth term
    epsilons[17] -= 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[4];

    // Twenty-sixth term
    epsilons[17] -= 6*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-seventh term
    epsilons[17] += 6*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-eighth term
    epsilons[17] += (3*gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twenty-ninth term
    epsilons[17] -= 6*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Thirtieth term
    epsilons[17] += 11*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Thirty-first term
    epsilons[17] += 7*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Thirty-second term
    epsilons[17] -= 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Thirty-third term
    epsilons[17] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[10];

    // Thirty-fourth term
    epsilons[17] += beta_Q * a_eval[ind] * sigma * deltas[4];

    // Thirty-fifth term
    epsilons[17] -= 3*beta_Q * a_eval[ind] * nu * deltas[10];

    // Thirty-sixth term
    epsilons[17] += 12*beta_Q * a_eval[ind] * beta * gammas[1];

    // Thirty-seventh term
    epsilons[17] += 6*beta_Q * a_eval[ind] * beta * gammas[2];

    // Thirty-eigth term
    epsilons[17] -= 12*beta_Q * a_eval[ind] * alpha * gammas[3];

    // Thirty-ninth term
    epsilons[17] -= 5*beta_Q * a_eval[ind] * alpha * gammas[4];
    break;
  }
  case 19:
  {
    // First term
    epsilons[18] += epsilons_Q[18] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[18] -= (12*epsilons_Q[28] - epsilons_Q[29] + 4*epsilons_Q[30] - 3*epsilons_Q[31] + epsilons_Q[33] - 2*epsilons_Q[34] + 5*epsilons_Q[35] + 6*epsilons_Q[39] + 5*epsilons_Q[40] + 2*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[18] -= (12*epsilons_Q[42] + 9*epsilons_Q[43] + 4*epsilons_Q[44] + 13*epsilons_Q[45] + 12*epsilons_Q[46] + 5*epsilons_Q[47] + 5*epsilons_Q[48] + 12*epsilons_Q[49] + 9*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[18] -= (12*epsilons_Q[51] + 13*epsilons_Q[52] + 22*epsilons_Q[53] + 23*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[18] -= 26*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[18] -= (2*deltas_Q[4] + 5*deltas_Q[5] - 2*deltas_Q[7] - 4*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[18] += (4*deltas_Q[9] + deltas_Q[10] - 6*deltas_Q[11] - 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[18] += (7*deltas_Q[9] + deltas_Q[10] - 17*deltas_Q[11] - 4*deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[18] += (7*deltas_Q[14] - 8*deltas_Q[15] - 9*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[18] += (5*deltas_Q[14] - 15*deltas_Q[15] - 7*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[18] -= 15*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[18] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[18] -= 12*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[18] += (5*gammas_Q[1] - 4*gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fifteenth term
    epsilons[18] -= 2*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[4];

    // Sixteenth term
    epsilons[18] -= 12*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Seventeenth term
    epsilons[18] += 3*(gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Eighteenth term
    epsilons[18] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Nineteenth term
    epsilons[18] -= (10*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Twentieth term
    epsilons[18] += 4*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[4];

    // Twenty-first term
    epsilons[18] += 12*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Twenty-second term
    epsilons[18] -= 16*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twenty-third term
    epsilons[18] += 15*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-fourth term
    epsilons[18] += 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-fifth term
    epsilons[18] -= 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-sixth term
    epsilons[18] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[4];

    // Twenty-seventh term
    epsilons[18] += 12*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-eighth term
    epsilons[18] += 2*(gammas_Q[1] - 2*gammas_Q[2]) * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-ninth term
    epsilons[18] -= (10*gammas_Q[3] + 6*gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Thirtieth term
    epsilons[18] -= 2*(2*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Thirty-first term
    epsilons[18] -= 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Thirty-second term
    epsilons[18] -= 22*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Thirty-third term
    epsilons[18] += 2*alpha_Q * pow(a_eval[ind], 2) * deltas[9];

    // Thirty-fourth term
    epsilons[18] -= 5*alpha_Q * pow(a_eval[ind], 2) * deltas[10];

    // Thirty-fifth term
    epsilons[18] += beta_Q * a_eval[ind] * sigma * deltas[5];

    // Thirty-sixth term
    epsilons[18] -= 2*beta_Q * a_eval[ind] * nu * deltas[9];

    // Thirty-seventh term
    epsilons[18] += 5*beta_Q * a_eval[ind] * nu * deltas[10];

    // Thirty-eigth term
    epsilons[18] -= 10*beta_Q * a_eval[ind] * beta * gammas[1];

    // Thirty-ninth term
    epsilons[18] -= 4*beta_Q * a_eval[ind] * beta * gammas[2];

    // Fourtieth term
    epsilons[18] += 10*beta_Q * a_eval[ind] * alpha * gammas[3];

    // Forty-first term
    epsilons[18] += 2*beta_Q * a_eval[ind] * alpha * gammas[4];
    break;
  }
  case 20:
  {
    // First term
    epsilons[19] += epsilons_Q[19] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[19] += (8*epsilons_Q[28] - epsilons_Q[29] + 3*epsilons_Q[30] - epsilons_Q[31] + 4*epsilons_Q[35] - 6*epsilons_Q[39] - 4*epsilons_Q[40] - 3*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[19] += (10*epsilons_Q[42] + 7*epsilons_Q[43] + 3*epsilons_Q[44] + 12*epsilons_Q[45] + 8*epsilons_Q[46] + 4*epsilons_Q[47] - 7*epsilons_Q[48] - 5*epsilons_Q[49] - 13*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[19] += (10*epsilons_Q[51] + 12*epsilons_Q[52] - 7*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[19] -= 4*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[19] += (deltas_Q[5] - deltas_Q[6] + 2*deltas_Q[7] + 4*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[19] += 2*(deltas_Q[11] - 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[19] += (2*deltas_Q[10] - 4*deltas_Q[11] - 6*deltas_Q[12] - 5*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[19] += 2*(deltas_Q[15] - 3*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[19] += 3*(deltas_Q[14] - 2*deltas_Q[15] - 2*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[19] -= 16*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[19] -= 7*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[19] += 8*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[19] -= (4*gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fifteenth term
    epsilons[19] += 8*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Sixteenth term
    epsilons[19] += 4*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Seventeenth term
    epsilons[19] += (3*gammas_Q[3] + 4*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Eighteenth term
    epsilons[19] += (8*gammas_Q[3] + 5*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Nineteenth term
    epsilons[19] -= 8*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Twentieth term
    epsilons[19] -= 16*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twenty-first term
    epsilons[19] += 8*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-second term
    epsilons[19] += 11*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-third term
    epsilons[19] += 9*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-fourth term
    epsilons[19] -= 8*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-fifth term
    epsilons[19] -= 4*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-sixth term
    epsilons[19] += 2*(3*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twenty-seventh term
    epsilons[19] += 4*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-eighth term
    epsilons[19] += 10*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Twenty-ninth term
    epsilons[19] += 4*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Thirtieth term
    epsilons[19] += 4*alpha_Q * pow(a_eval[ind], 2) * deltas[10];

    // Thirty-first term
    epsilons[19] += beta_Q * a_eval[ind] * sigma * deltas[6];

    // Thirty-second term
    epsilons[19] -= 4*beta_Q * a_eval[ind] * nu * deltas[10];

    // Thirty-third term
    epsilons[19] -= 8*beta_Q * a_eval[ind] * beta * gammas[1];

    // Thirty-fourth term
    epsilons[19] -= 6*beta_Q * a_eval[ind] * beta * gammas[2];

    // Thirty-fifth term
    epsilons[19] -= beta_Q * a_eval[ind] * alpha * gammas[3];
    break;
  }
  case 21:
  {
    // First term
    epsilons[20] += epsilons_Q[20] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[20] -= (2*epsilons_Q[28] + epsilons_Q[30] + epsilons_Q[32] + 2*epsilons_Q[33] + epsilons_Q[35]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[20] -= (3*epsilons_Q[42] + 3*epsilons_Q[43] + 3*epsilons_Q[44] + 3*epsilons_Q[45] + 2*epsilons_Q[46] + epsilons_Q[47] + epsilons_Q[49]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[20] -= 3*(2*epsilons_Q[51] + epsilons_Q[52] + epsilons_Q[53] + epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[20] -= 6*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[20] += (deltas_Q[6] + deltas_Q[7]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[20] -= (2*deltas_Q[9] + deltas_Q[10]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[20] -= (deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[20] -= (5*deltas_Q[14] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[20] -= (2*deltas_Q[14] + 2*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[20] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[20] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[20] -= 2*gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[20] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fifteenth term
    epsilons[20] -= 2*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Sixteenth term
    epsilons[20] -= 3*gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Seventeenth term
    epsilons[20] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Eighteenth term
    epsilons[20] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Nineteenth term
    epsilons[20] += 2*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twentieth term
    epsilons[20] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-first term
    epsilons[20] -= 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-second term
    epsilons[20] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-third term
    epsilons[20] -= (3*gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twenty-fourth term
    epsilons[20] -= 7*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Twenty-fifth term
    epsilons[20] -= alpha_Q * pow(a_eval[ind], 2) * deltas[10];

    // Twenty-sixth term
    epsilons[20] += beta_Q * a_eval[ind] * sigma * deltas[7];

    // Twenty-seventh term
    epsilons[20] += beta_Q * a_eval[ind] * nu * deltas[10];

    // Twenty-eighth term
    epsilons[20] -= beta_Q * a_eval[ind] * beta * gammas[2];

    // Twenty-ninth term
    epsilons[20] -= beta_Q * a_eval[ind] * alpha * gammas[3];
    break;
  }
  case 22:
  {
    // First term
    epsilons[21] += epsilons_Q[21] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[21] += (epsilons_Q[28] - epsilons_Q[29] - 2*epsilons_Q[31] - epsilons_Q[32] - epsilons_Q[34] + 3*epsilons_Q[39] + 2*epsilons_Q[40] + epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[21] -= (epsilons_Q[42] + epsilons_Q[43] + epsilons_Q[44] + epsilons_Q[45] - epsilons_Q[46] - 2*epsilons_Q[48] - 4*epsilons_Q[49] - 4*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[21] -= (2*epsilons_Q[51] + epsilons_Q[52] - 4*epsilons_Q[53] - 6*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[21] += 5*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[21] += (deltas_Q[4] + deltas_Q[5] - 2*deltas_Q[6] - deltas_Q[7] - deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[21] -= (deltas_Q[9] - 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[21] += (4*deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[21] -= (deltas_Q[14] - deltas_Q[15] - 4*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[21] += 3*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[21] += 9*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[21] += 4*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[21] += gammas_Q[0] * pow(a_eval[ind], 4) * gammas[5];

    // Fourteenth term
    epsilons[21] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[4];

    // Fifteenth term
    epsilons[21] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * gammas[5];

    // Sixteenth term
    epsilons[21] -= (3*gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Seventeenth term
    epsilons[21] -= gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Eighteenth term
    epsilons[21] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[4];

    // Nineteenth term
    epsilons[21] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * gammas[5];

    // Twentieth term
    epsilons[21] += 7*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twenty-first term
    epsilons[21] -= 7*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-second term
    epsilons[21] -= 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-third term
    epsilons[21] -= 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[4];

    // Twenty-fourth term
    epsilons[21] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * gammas[5];

    // Twenty-fifth term
    epsilons[21] -= (gammas_Q[1] - gammas_Q[2]) * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-sixth term
    epsilons[21] += (5*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twenty-seventh term
    epsilons[21] += 2*gammas_Q[3] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-eighth term
    epsilons[21] += 3*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-ninth term
    epsilons[21] += 7*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Thirtieth term
    epsilons[21] += gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Thirty-first term
    epsilons[21] -= alpha_Q * pow(a_eval[ind], 2) * deltas[9];

    // Thirty-second term
    epsilons[21] += beta_Q * a_eval[ind] * sigma * deltas[8];

    // Thirty-third term
    epsilons[21] += beta_Q * a_eval[ind] * nu * deltas[9];

    // Thirty-fourth term
    epsilons[21] += 3*beta_Q * a_eval[ind] * beta * gammas[1];

    // Thirty-fifth term
    epsilons[21] += 2*beta_Q * a_eval[ind] * beta * gammas[2];

    // Thirty-sixth term
    epsilons[21] -= beta_Q * a_eval[ind] * alpha * gammas[4];
    break;
  }
  case 23:
  {
    // First term
    epsilons[22] += epsilons_Q[22] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[22] -= (epsilons_Q[34] + 3*epsilons_Q[35] + 3*epsilons_Q[36] + epsilons_Q[37] + 7*epsilons_Q[39] + 4*epsilons_Q[40] + 2*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[22] -= (5*epsilons_Q[45] + 5*epsilons_Q[46] + 3*epsilons_Q[47] + 7*epsilons_Q[48] + 11*epsilons_Q[49] + 14*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[22] -= (5*epsilons_Q[52] + 15*epsilons_Q[53] + 18*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[22] -= 18*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[22] += (2*deltas_Q[7] + 5*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[22] -= 2*(2*deltas_Q[11] + deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[22] -= (10*deltas_Q[11] + 7*deltas_Q[12] + 4*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[22] -= (8*deltas_Q[15] + 12*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[22] -= (12*deltas_Q[15] + 9*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[22] -= 26*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[22] -= 11*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[22] -= 6*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fourteenth term
    epsilons[22] -= 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[4];

    // Fifteenth term
    epsilons[22] += 4*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Sixteenth term
    epsilons[22] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Seventeenth term
    epsilons[22] += 6*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Eighteenth term
    epsilons[22] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[4];

    // Nineteenth term
    epsilons[22] -= 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Twentieth term
    epsilons[22] += 8*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Twenty-first term
    epsilons[22] += 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twenty-second term
    epsilons[22] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-third term
    epsilons[22] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[4];

    // Twenty-fourth term
    epsilons[22] -= 4*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-fifth term
    epsilons[22] += 4*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-sixth term
    epsilons[22] -= 6*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-seventh term
    epsilons[22] += 4*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Twenty-eighth term
    epsilons[22] -= 2*beta_Q * a_eval[ind] * sigma * deltas[4];

    // Twenty-ninth term
    epsilons[22] -= 8*beta_Q * a_eval[ind] * beta * gammas[1];

    // Thirtieth term
    epsilons[22] -= 4*beta_Q * a_eval[ind] * beta * gammas[2];

    // Thirty-first term
    epsilons[22] += 6*beta_Q * a_eval[ind] * alpha * gammas[3];

    // Thirty-second term
    epsilons[22] += 2*beta_Q * a_eval[ind] * alpha * gammas[4];
    break;
  }
  case 24:
  {
    // First term
    epsilons[23] += epsilons_Q[23] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[23] -= (epsilons_Q[34] - epsilons_Q[35] - 3*epsilons_Q[36] + epsilons_Q[38] - 9*epsilons_Q[39] - 5*epsilons_Q[40] - 3*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[23] += (epsilons_Q[45] + 3*epsilons_Q[46] + epsilons_Q[47] + 8*epsilons_Q[48] + 13*epsilons_Q[49] + 18*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[23] += (epsilons_Q[52] + 16*epsilons_Q[53] + 21*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[23] += 19*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[23] -= 3*(deltas_Q[7] + 3*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[23] += (6*deltas_Q[11] + 2*deltas_Q[12] + 4*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[23] += (18*deltas_Q[11] + 12*deltas_Q[12] + 6*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[23] += (10*deltas_Q[15] + 14*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[23] += (21*deltas_Q[15] + 15*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[23] += 30*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[23] += 18*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[23] += 6*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fourteenth term
    epsilons[23] -= 5*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Fifteenth term
    epsilons[23] -= 3*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Sixteenth term
    epsilons[23] -= 6*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Seventeenth term
    epsilons[23] += 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Eighteenth term
    epsilons[23] -= 10*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Nineteenth term
    epsilons[23] -= 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Twentieth term
    epsilons[23] -= 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Twenty-first term
    epsilons[23] += 6*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Twenty-second term
    epsilons[23] += gammas_Q[4] * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twenty-third term
    epsilons[23] -= 6*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-fourth term
    epsilons[23] += 10*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-fifth term
    epsilons[23] += 2*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Twenty-sixth term
    epsilons[23] -= 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Twenty-seventh term
    epsilons[23] -= 2*beta_Q * a_eval[ind] * sigma * deltas[5];

    // Twenty-eighth term
    epsilons[23] += 12*beta_Q * a_eval[ind] * beta * gammas[1];

    // Twenty-ninth term
    epsilons[23] += 6*beta_Q * a_eval[ind] * beta * gammas[2];

    // Thirtieth term
    epsilons[23] -= 6*beta_Q * a_eval[ind] * alpha * gammas[3];
    break;
  }
  case 25:
  {
    // First term
    epsilons[24] += epsilons_Q[24] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[24] -= (epsilons_Q[34] + epsilons_Q[35] + epsilons_Q[36] + epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[24] -= (2*epsilons_Q[45] + epsilons_Q[46] + epsilons_Q[47] + epsilons_Q[48] + epsilons_Q[49] + 2*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[24] -= 2*(epsilons_Q[52] + epsilons_Q[53] + epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[24] -= 3*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[24] += (deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[24] -= 4*deltas_Q[11] * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[24] -= 2*(deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[24] -= 4*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[24] -= 3*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[24] -= 8*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[24] -= 4*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[24] -= 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[3];

    // Fourteenth term
    epsilons[24] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Fifteenth term
    epsilons[24] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[3];

    // Sixteenth term
    epsilons[24] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Seventeenth term
    epsilons[24] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Eighteenth term
    epsilons[24] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[3];

    // Nineteenth term
    epsilons[24] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twentieth term
    epsilons[24] -= 6*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-first term
    epsilons[24] -= 4*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Twenty-second term
    epsilons[24] -= 2*beta_Q * a_eval[ind] * sigma * deltas[6];

    // Twenty-third term
    epsilons[24] += 2*beta_Q * a_eval[ind] * alpha * gammas[3];
    break;
  }
  case 26:
  {
    // First term
    epsilons[25] += epsilons_Q[25] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[25] -= (epsilons_Q[37] + 2*epsilons_Q[38] - 6*epsilons_Q[39] - 4*epsilons_Q[40] - 2*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[25] -= (epsilons_Q[46] + 2*epsilons_Q[47] - 6*epsilons_Q[48] - 8*epsilons_Q[49] - 14*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[25] -= (3*epsilons_Q[52] - 10*epsilons_Q[53] - 15*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[25] += 13*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[25] -= 4*(deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[25] -= 2*(deltas_Q[12] - 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[25] += 8*(deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[25] -= 4*(deltas_Q[15] - 2*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[25] += 12*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[25] += 22*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[25] += 16*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[25] -= 4*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Fourteenth term
    epsilons[25] -= 6*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Fifteenth term
    epsilons[25] += 16*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Sixteenth term
    epsilons[25] -= 8*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Seventeenth term
    epsilons[25] -= 12*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Eighteenth term
    epsilons[25] += 4*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Nineteenth term
    epsilons[25] += 4*gammas_Q[4] * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twentieth term
    epsilons[25] -= 4*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-first term
    epsilons[25] += 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-second term
    epsilons[25] += 8*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Twenty-third term
    epsilons[25] -= 4*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Twenty-fourth term
    epsilons[25] -= 2*beta_Q * a_eval[ind] * sigma * deltas[7];

    // Twenty-fifth term
    epsilons[25] += 8*beta_Q * a_eval[ind] * beta * gammas[1];

    // Twenty-sixth term
    epsilons[25] += 8*beta_Q * a_eval[ind] * beta * gammas[2];
    break;
  }
  case 27:
  {
    // First term
    epsilons[26] += epsilons_Q[26] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[26] -= (2*epsilons_Q[36] + epsilons_Q[37] + 9*epsilons_Q[39] + 6*epsilons_Q[40] + 3*epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[26] -= (2*epsilons_Q[46] + epsilons_Q[47] + 9*epsilons_Q[48] + 12*epsilons_Q[49] + 17*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[26] -= (epsilons_Q[52] + 15*epsilons_Q[53] + 20*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[26] -= 19*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[26] += (deltas_Q[7] + 4*deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[26] -= 4*deltas_Q[13] * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[26] -= (8*deltas_Q[11] + 5*deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[26] -= 8*deltas_Q[16] * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[26] -= (9*deltas_Q[15] + 6*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[26] -= 20*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[26] -= 7*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[26] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Fourteenth term
    epsilons[26] += 3*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Fifteenth term
    epsilons[26] -= 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Sixteenth term
    epsilons[26] += 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Seventeenth term
    epsilons[26] += 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Eighteenth term
    epsilons[26] -= 4*gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Nineteenth term
    epsilons[26] += 4*gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twentieth term
    epsilons[26] -= 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-first term
    epsilons[26] += 4*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Twenty-second term
    epsilons[26] -= 2*beta_Q * a_eval[ind] * sigma * deltas[8];

    // Twenty-third term
    epsilons[26] -= 8*beta_Q * a_eval[ind] * beta * gammas[1];

    // Twenty-fourth term
    epsilons[26] -= 6*beta_Q * a_eval[ind] * beta * gammas[2];
    break;
  }
  case 28:
  {
    // First term
    epsilons[27] += epsilons_Q[27] * pow(a_eval[ind], 5) * pow(sigma, 4);

    // Second term
    epsilons[27] -= 2*(epsilons_Q[39] + epsilons_Q[40] + epsilons_Q[41]) * pow(a_eval[ind], 4) * nu * pow(sigma, 4);

    // Third term
    epsilons[27] -= (3*epsilons_Q[48] + 3*epsilons_Q[49] + 5*epsilons_Q[50]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    epsilons[27] -= (4*epsilons_Q[53] + 6*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 4);

    // Fifth term
    epsilons[27] -= 7*epsilons_Q[55] * a_eval[ind] * pow(nu, 4) * pow(sigma, 4);

    // Sixth term
    epsilons[27] += (deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(sigma, 2) * beta;

    // Seventh term
    epsilons[27] -= 2*deltas_Q[13] * pow(a_eval[ind], 3) * pow(sigma, 3) * alpha;

    // Eighth term
    epsilons[27] -= 2*(deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 2) * beta;

    // Ninth term
    epsilons[27] -= 4*deltas_Q[16] * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * alpha;

    // Tenth term
    epsilons[27] -= 3*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 2) * beta;

    // Eleventh term
    epsilons[27] -= 10*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * alpha;

    // Twelfth term
    epsilons[27] -= 4*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 2) * beta;

    // Thirteenth term
    epsilons[27] -= gammas_Q[2] * pow(a_eval[ind], 3) * pow(beta, 2);

    // Fourteenth term
    epsilons[27] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[1];

    // Fifteenth term
    epsilons[27] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[2];

    // Sixteenth term
    epsilons[27] -= 6*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[0];

    // Seventeenth term
    epsilons[27] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[1];

    // Eighteenth term
    epsilons[27] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[2];

    // Nineteenth term
    epsilons[27] -= gammas_Q[4] * pow(a_eval[ind], 2) * sigma * alpha * beta;

    // Twentieth term
    epsilons[27] += gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(beta, 2);

    // Twenty-first term
    epsilons[27] -= 3*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(alpha, 2);

    // Twenty-second term
    epsilons[27] -= 2*gammas_Q[5] * a_eval[ind] * nu * sigma * alpha * beta;

    // Twenty-third term
    epsilons[27] += gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(beta, 2);

    // Twenty-fourth term
    epsilons[27] -= 2*beta_Q * a_eval[ind] * beta * gammas[1];

    // Twenty-fifth term
    epsilons[27] -= 2*beta_Q * a_eval[ind] * beta * gammas[2];
    break;
  }
  case 29:
  {
    // First term
    epsilons[28] += epsilons_Q[28] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[28] += (3*epsilons_Q[42] + epsilons_Q[43] + epsilons_Q[44] + 3*epsilons_Q[45] + epsilons_Q[46] + epsilons_Q[47]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[28] += (7*epsilons_Q[51] + 7*epsilons_Q[52] + 3*epsilons_Q[53] + epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[28] += 7*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[28] += (deltas_Q[9] + deltas_Q[10]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[28] += (deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[28] += (4*deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[28] += 5*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[28] += 4*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[28] -= gammas_Q[1] * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[28] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[28] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[28] += (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[28] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Fifteenth term
    epsilons[28] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Sixteenth term
    epsilons[28] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Seventeenth term
    epsilons[28] -= gammas_Q[3] * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Eighteenth term
    epsilons[28] -= gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Nineteenth term
    epsilons[28] -= gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twentieth term
    epsilons[28] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-first term
    epsilons[28] += alpha_Q * pow(a_eval[ind], 2) * deltas[15];

    // Twenty-second term
    epsilons[28] += alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-third term
    epsilons[28] -= 3*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fourth term
    epsilons[28] -= beta_Q * a_eval[ind] * nu * deltas[15];

    // Twenty-fifth term
    epsilons[28] -= beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-sixth term
    epsilons[28] += beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-seventh term
    epsilons[28] += beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-eighth term
    epsilons[28] -= beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 30:
  {
    // First term
    epsilons[29] += epsilons_Q[29] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[29] += (4*epsilons_Q[42] - epsilons_Q[43] + 3*epsilons_Q[44] + 3*epsilons_Q[45] + 3*epsilons_Q[47] + epsilons_Q[48]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[29] += (19*epsilons_Q[51] + 18*epsilons_Q[52] + 4*epsilons_Q[53] - epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[29] += 19*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[29] += (5*deltas_Q[9] + deltas_Q[10]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[29] += (deltas_Q[14] + 2*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[29] += (8*deltas_Q[14] - 3*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[29] += 14*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[29] += 8*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[29] -= 3*gammas_Q[1] * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[29] -= 3*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[29] -= (3*gammas_Q[3] + 4*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[29] += 3*(2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[29] -= 9*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Fifteenth term
    epsilons[29] -= 11*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Sixteenth term
    epsilons[29] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Seventeenth term
    epsilons[29] -= 6*gammas_Q[3] * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Eighteenth term
    epsilons[29] -= 5*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Nineteenth term
    epsilons[29] -= 6*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twentieth term
    epsilons[29] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-first term
    epsilons[29] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-second term
    epsilons[29] += beta_Q * a_eval[ind] * sigma * deltas[9];

    // Twenty-third term
    epsilons[29] -= 3*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fourth term
    epsilons[29] -= 3*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-fifth term
    epsilons[29] += 3*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-sixth term
    epsilons[29] += 4*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-seventh term
    epsilons[29] -= 2*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 31:
  {
    // First term
    epsilons[30] += epsilons_Q[30] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[30] -= (3*epsilons_Q[42] - 2*epsilons_Q[43] + epsilons_Q[44] + 6*epsilons_Q[45] + 2*epsilons_Q[47] - epsilons_Q[49]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[30] -= (13*epsilons_Q[51] + 16*epsilons_Q[52] + 3*epsilons_Q[53] - 2*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[30] -= 13*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[30] -= (4*deltas_Q[9] + 3*deltas_Q[10]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[30] -= 4*(deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[30] -= (13*deltas_Q[14] + 2*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[30] -= 20*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[30] -= 13*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[30] += 6*gammas_Q[1] * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[30] += (5*gammas_Q[3] + 4*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[30] += 4*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[30] -= 6*(2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[30] += 13*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Fifteenth term
    epsilons[30] += 12*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Sixteenth term
    epsilons[30] -= 12*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Seventeenth term
    epsilons[30] += 4*gammas_Q[3] * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Eighteenth term
    epsilons[30] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Nineteenth term
    epsilons[30] += 4*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twentieth term
    epsilons[30] -= 6*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-first term
    epsilons[30] -= 2*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-second term
    epsilons[30] += beta_Q * a_eval[ind] * sigma * deltas[10];

    // Twenty-third term
    epsilons[30] += 6*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fourth term
    epsilons[30] += 2*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-fifth term
    epsilons[30] -= 4*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-sixth term
    epsilons[30] -= 4*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-seventh term
    epsilons[30] += 4*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 32:
  {
    // First term
    epsilons[31] += epsilons_Q[31] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[31] += (epsilons_Q[42] + epsilons_Q[43] + epsilons_Q[44] + epsilons_Q[45] + epsilons_Q[50]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[31] += (3*epsilons_Q[51] + epsilons_Q[52] + epsilons_Q[53] + epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[31] += 3*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[31] += (deltas_Q[9] + deltas_Q[10] + deltas_Q[11]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[31] += 5*deltas_Q[14] * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[31] += 2*(2*deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[31] += 5*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[31] += 6*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[31] -= (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[31] -= 3*gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[31] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Thirteenth term
    epsilons[31] -= gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Fourteenth term
    epsilons[31] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Fifteenth term
    epsilons[31] += 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Sixteenth term
    epsilons[31] += (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Seventeenth term
    epsilons[31] += 3*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Eighteenth term
    epsilons[31] += 3*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Nineteenth term
    epsilons[31] += alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twentieth term
    epsilons[31] += beta_Q * a_eval[ind] * sigma * deltas[11];

    // Twenty-first term
    epsilons[31] -= beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-second term
    epsilons[31] -= beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-third term
    epsilons[31] += beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 33:
  {
    // First term
    epsilons[32] += epsilons_Q[32] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[32] -= (3*epsilons_Q[42] - 2*epsilons_Q[43] + 3*epsilons_Q[44] + 3*epsilons_Q[45] + 3*epsilons_Q[47]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[32] -= (19*epsilons_Q[51] + 18*epsilons_Q[52] + 3*epsilons_Q[53] - 2*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[32] -= 19*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[32] -= (9*deltas_Q[9] + 3*deltas_Q[10] - deltas_Q[12]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[32] -= (11*deltas_Q[14] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[32] -= (18*deltas_Q[14] - 5*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[32] -= 26*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[32] -= 19*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[32] += 3*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[32] += 3*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[32] += 3*(gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[32] -= 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[32] += 3*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Fifteenth term
    epsilons[32] += 9*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Sixteenth term
    epsilons[32] += 15*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Seventeenth term
    epsilons[32] -= 9*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Eighteenth term
    epsilons[32] += (9*gammas_Q[3] - 2*gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Nineteenth term
    epsilons[32] += gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Twentieth term
    epsilons[32] += 5*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-first term
    epsilons[32] -= 3*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-second term
    epsilons[32] -= 3*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-third term
    epsilons[32] += beta_Q * a_eval[ind] * sigma * deltas[12];

    // Twenty-fourth term
    epsilons[32] += 3*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fifth term
    epsilons[32] += 3*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-sixth term
    epsilons[32] -= 3*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-seventh term
    epsilons[32] -= 7*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-eighth term
    epsilons[32] -= 3*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 34:
  {
    // First term
    epsilons[33] += epsilons_Q[33] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[33] += (3*epsilons_Q[42] + 4*epsilons_Q[44] + 3*epsilons_Q[45] + 2*epsilons_Q[47]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[33] += (18*epsilons_Q[51] + 13*epsilons_Q[52] + 3*epsilons_Q[53]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[33] += 18*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[33] += (5*deltas_Q[9] + 2*deltas_Q[10] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[33] += (9*deltas_Q[14] + 2*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[33] += (11*deltas_Q[14] - deltas_Q[15] + 4*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[33] += 19*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[33] += 16*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[33] -= 3*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[33] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[33] -= 4*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[33] += 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[33] -= gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[1];

    // Fifteenth term
    epsilons[33] -= 3*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Sixteenth term
    epsilons[33] -= 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Seventeenth term
    epsilons[33] -= 12*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Eighteenth term
    epsilons[33] += 9*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Nineteenth term
    epsilons[33] -= 7*gammas_Q[3] * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Twentieth term
    epsilons[33] -= 7*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-first term
    epsilons[33] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-second term
    epsilons[33] += 2*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-third term
    epsilons[33] += beta_Q * a_eval[ind] * sigma * deltas[13];

    // Twenty-fourth term
    epsilons[33] -= 3*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fifth term
    epsilons[33] -= 2*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-sixth term
    epsilons[33] += 2*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-seventh term
    epsilons[33] += 4*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-eighth term
    epsilons[33] += 3*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 35:
  {
    // First term
    epsilons[34] += epsilons_Q[34] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[34] -= (2*epsilons_Q[42] + epsilons_Q[44] + epsilons_Q[45] + epsilons_Q[46] + epsilons_Q[47] - epsilons_Q[48] - 2*epsilons_Q[50]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[34] -= (7*epsilons_Q[51] + 6*epsilons_Q[52] + epsilons_Q[53] - epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[34] -= 5*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[34] -= (2*deltas_Q[9] + deltas_Q[10] + 2*deltas_Q[11] - deltas_Q[12] - deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[34] -= (2*deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[34] -= (5*deltas_Q[14] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[34] -= 5*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[34] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[34] += 2*gammas_Q[1] * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[34] += (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[34] += (2*gammas_Q[3] + 3*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[34] -= 2*(2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[34] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[1];

    // Fifteenth term
    epsilons[34] += 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Sixteenth term
    epsilons[34] += 8*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Seventeenth term
    epsilons[34] -= 4*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Eighteenth term
    epsilons[34] += (3*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Nineteenth term
    epsilons[34] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Twentieth term
    epsilons[34] += gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-first term
    epsilons[34] -= 2*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-second term
    epsilons[34] -= alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-third term
    epsilons[34] -= 2*beta_Q * a_eval[ind] * sigma * deltas[9];

    // Twenty-fourth term
    epsilons[34] += 2*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fifth term
    epsilons[34] += beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-sixth term
    epsilons[34] -= beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-seventh term
    epsilons[34] -= 2*beta_Q * a_eval[ind] * beta * gammas[4];
    break;
  }
  case 36:
  {
    // First term
    epsilons[35] += epsilons_Q[35] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[35] += (3*epsilons_Q[45] + 2*epsilons_Q[46] + epsilons_Q[47] + epsilons_Q[49]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[35] += 3*(epsilons_Q[52] + epsilons_Q[53] + epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[35] += 6*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[35] += deltas_Q[11] * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[35] += 2*(deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[35] += (deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[35] += 8*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[35] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[35] += 4*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[35] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[35] -= 4*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Thirteenth term
    epsilons[35] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Fourteenth term
    epsilons[35] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Fifteenth term
    epsilons[35] -= 4*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Sixteenth term
    epsilons[35] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Seventeenth term
    epsilons[35] -= 2*beta_Q * a_eval[ind] * sigma * deltas[10];

    // Eighteenth term
    epsilons[35] -= 4*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 37:
  {
    // First term
    epsilons[36] += epsilons_Q[36] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[36] += (4*epsilons_Q[42] + 2*epsilons_Q[44] + 4*epsilons_Q[45] + epsilons_Q[46] + 3*epsilons_Q[47] + epsilons_Q[50]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[36] += (14*epsilons_Q[51] + 16*epsilons_Q[52] + 4*epsilons_Q[53] - 3*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[36] += 4*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[36] += (4*deltas_Q[9] + 2*deltas_Q[10] - 2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[36] += 4*(deltas_Q[14] + deltas_Q[15]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[36] += 10*deltas_Q[14] * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[36] += 5*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Ninth term
    epsilons[36] -= 4*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Tenth term
    epsilons[36] -= 2*(gammas_Q[3] + 3*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Eleventh term
    epsilons[36] -= 4*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Twelfth term
    epsilons[36] += 8*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Thirteenth term
    epsilons[36] -= 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Fourteenth term
    epsilons[36] -= 14*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Fifteenth term
    epsilons[36] -= 12*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Sixteenth term
    epsilons[36] += 12*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Seventeenth term
    epsilons[36] -= 2*(3*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Eighteenth term
    epsilons[36] -= 20*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Nineteenth term
    epsilons[36] -= 10*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twentieth term
    epsilons[36] += 4*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-first term
    epsilons[36] += 2*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-second term
    epsilons[36] -= 2*beta_Q * a_eval[ind] * sigma * deltas[11];

    // Twenty-third term
    epsilons[36] -= 4*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fourth term
    epsilons[36] -= 2*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-fifth term
    epsilons[36] += 6*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-sixth term
    epsilons[36] += 4*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-seventh term
    epsilons[36] += 4*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 38:
  {
    // First term
    epsilons[37] += epsilons_Q[37] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[37] += 2*(2*epsilons_Q[42] + epsilons_Q[44] + 2*epsilons_Q[45] + epsilons_Q[46] + epsilons_Q[47] - 2*epsilons_Q[50]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[37] += (14*epsilons_Q[51] + 13*epsilons_Q[52] + 4*epsilons_Q[53] + 3*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[37] += 26*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[37] += 2*(2*deltas_Q[9] + deltas_Q[10] - 3*deltas_Q[12] - deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[37] += 2*(2*deltas_Q[14] - 4*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[37] += 2*(5*deltas_Q[14] - 6*deltas_Q[15] - deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[37] += 22*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[37] += 6*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[37] -= 4*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[37] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[37] -= 4*(gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[37] += 8*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[37] -= 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[1];

    // Fifteenth term
    epsilons[37] += 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Sixteenth term
    epsilons[37] -= 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Seventeenth term
    epsilons[37] -= 20*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Eighteenth term
    epsilons[37] += 12*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Nineteenth term
    epsilons[37] -= 6*(gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Twentieth term
    epsilons[37] += 24*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Twenty-first term
    epsilons[37] += 6*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-second term
    epsilons[37] += 4*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-third term
    epsilons[37] += 2*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-fourth term
    epsilons[37] -= 2*beta_Q * a_eval[ind] * sigma * deltas[12];

    // Twenty-fifth term
    epsilons[37] -= 4*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-sixth term
    epsilons[37] -= 2*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-seventh term
    epsilons[37] += 2*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-eighth term
    epsilons[37] += 8*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-ninth term
    epsilons[37] += 4*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 39:
  {
    // First term
    epsilons[38] += epsilons_Q[38] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[38] -= (6*epsilons_Q[42] + 3*epsilons_Q[44] + 6*epsilons_Q[45] + epsilons_Q[47]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[38] -= (21*epsilons_Q[51] + 16*epsilons_Q[52] + 6*epsilons_Q[53] + 3*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[38] -= 25*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[38] -= (6*deltas_Q[9] + 3*deltas_Q[10] - deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[38] -= 3*(2*deltas_Q[14] - 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[38] -= (15*deltas_Q[14] - 2*deltas_Q[15] + 5*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[38] -= 27*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[38] -= 21*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[38] += 6*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[38] += 3*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[38] += 6*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[38] -= 12*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[38] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[1];

    // Fifteenth term
    epsilons[38] += 9*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Sixteenth term
    epsilons[38] += 18*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Seventeenth term
    epsilons[38] -= 18*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Eighteenth term
    epsilons[38] += (9*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Nineteenth term
    epsilons[38] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Twentieth term
    epsilons[38] += 7*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-first term
    epsilons[38] -= 6*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-second term
    epsilons[38] -= 3*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-third term
    epsilons[38] -= 2*beta_Q * a_eval[ind] * sigma * deltas[13];

    // Twenty-fourth term
    epsilons[38] += 6*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fifth term
    epsilons[38] += 3*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-sixth term
    epsilons[38] -= 3*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-seventh term
    epsilons[38] -= 6*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-eighth term
    epsilons[38] -= 6*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 40:
  {
    // First term
    epsilons[39] += epsilons_Q[39] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[39] -= (6*epsilons_Q[42] + 3*epsilons_Q[44] + 6*epsilons_Q[45] + 3*epsilons_Q[47] - epsilons_Q[48] - 3*epsilons_Q[49] - 3*epsilons_Q[50]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[39] -= (21*epsilons_Q[51] + 21*epsilons_Q[52] + epsilons_Q[53] - 8*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[39] -= 7*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[39] -= (6*deltas_Q[9] + 3*deltas_Q[10] - 3*deltas_Q[11] - deltas_Q[12] - deltas_Q[13]) * pow(a_eval[ind], 3) * beta * pow(sigma, 3);

    // Sixth term
    epsilons[39] -= (6*deltas_Q[14] - deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[39] -= (15*deltas_Q[14] - 5*deltas_Q[15] - 2*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[39] -= deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[39] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[39] += 6*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[39] += 3*(gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[39] += (6*gammas_Q[3] + 7*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[39] -= 12*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[39] += 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Fifteenth term
    epsilons[39] += 15*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Sixteenth term
    epsilons[39] += 20*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Seventeenth term
    epsilons[39] -= 18*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Eighteenth term
    epsilons[39] += (9*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Nineteenth term
    epsilons[39] += 16*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Twentieth term
    epsilons[39] += 11*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-first term
    epsilons[39] -= 6*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-second term
    epsilons[39] -= 3*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-third term
    epsilons[39] += 6*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fourth term
    epsilons[39] += 3*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-fifth term
    epsilons[39] -= 9*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-sixth term
    epsilons[39] -= 8*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-seventh term
    epsilons[39] -= 6*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 41:
  {
    // First term
    epsilons[40] += epsilons_Q[40] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[40] += (8*epsilons_Q[42] + 4*epsilons_Q[44] + 8*epsilons_Q[45] + 4*epsilons_Q[47] + epsilons_Q[48] - epsilons_Q[49] + epsilons_Q[50]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[40] += (28*epsilons_Q[51] + 28*epsilons_Q[52] + 7*epsilons_Q[53] - 4*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[40] += 14*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[40] += (8*deltas_Q[9] + 4*deltas_Q[10] - 3*deltas_Q[11] - 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[40] += 2*(4*deltas_Q[14] - deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[40] += (20*deltas_Q[14] - 3*deltas_Q[15] - 3*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[40] -= 2*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[40] += 4*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[40] -= 8*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[40] -= (4*gammas_Q[3] + 7*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[40] -= 8*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[40] += 16*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[40] += 2*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[1];

    // Fifteenth term
    epsilons[40] -= 20*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Sixteenth term
    epsilons[40] -= 18*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Seventeenth term
    epsilons[40] -= 24*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Eighteenth term
    epsilons[40] += 24*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Nineteenth term
    epsilons[40] -= 3*(4*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Twentieth term
    epsilons[40] -= 30*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Twenty-first term
    epsilons[40] -= 18*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-second term
    epsilons[40] += 8*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-third term
    epsilons[40] += 4*alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-fourth term
    epsilons[40] -= 8*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fifth term
    epsilons[40] -= 4*beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-sixth term
    epsilons[40] += 10*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-seventh term
    epsilons[40] += 8*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-eighth term
    epsilons[40] += 8*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 42:
  {
    // First term
    epsilons[41] += epsilons_Q[41] * pow(a_eval[ind], 4) * pow(sigma, 5);

    // Second term
    epsilons[41] -= (2*epsilons_Q[42] + epsilons_Q[44] + 2*epsilons_Q[45] + epsilons_Q[47] - epsilons_Q[48] - epsilons_Q[49] - epsilons_Q[50]) * pow(a_eval[ind], 3) * nu * pow(sigma, 5);

    // Third term
    epsilons[41] -= (7*epsilons_Q[51] + 7*epsilons_Q[52] - 5*epsilons_Q[54]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 5);

    // Fourth term
    epsilons[41] += 7*epsilons_Q[55] * a_eval[ind] * pow(nu, 3) * pow(sigma, 5);

    // Fifth term
    epsilons[41] -= (2*deltas_Q[9] + deltas_Q[10] - deltas_Q[11] - deltas_Q[13]) * pow(a_eval[ind], 3) * pow(sigma, 3) * beta;

    // Sixth term
    epsilons[41] += (3*deltas_Q[16] - 2*deltas_Q[14]) * pow(a_eval[ind], 2) * pow(sigma, 4) * alpha;

    // Seventh term
    epsilons[41] -= (5*deltas_Q[14] - deltas_Q[15] - 2*deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 3) * beta;

    // Eighth term
    epsilons[41] += 13*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * alpha;

    // Ninth term
    epsilons[41] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 3) * beta;

    // Tenth term
    epsilons[41] += 2*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * gammas[5];

    // Eleventh term
    epsilons[41] += (gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[3];

    // Twelfth term
    epsilons[41] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[4];

    // Thirteenth term
    epsilons[41] -= 4*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * gammas[5];

    // Fourteenth term
    epsilons[41] += 10*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[2];

    // Fifteenth term
    epsilons[41] += 5*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[3];

    // Sixteenth term
    epsilons[41] += 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[4];

    // Seventeenth term
    epsilons[41] -= 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * gammas[5];

    // Eighteenth term
    epsilons[41] += (3*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * sigma * pow(beta, 2);

    // Nineteenth term
    epsilons[41] += 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * alpha * beta;

    // Twentieth term
    epsilons[41] += 5*gammas_Q[5] * a_eval[ind] * nu * sigma * pow(beta, 2);

    // Twenty-first term
    epsilons[41] -= 2*alpha_Q * pow(a_eval[ind], 2) * deltas[14];

    // Twenty-second term
    epsilons[41] -= alpha_Q * pow(a_eval[ind], 2) * deltas[16];

    // Twenty-third term
    epsilons[41] += 2*beta_Q * a_eval[ind] * nu * deltas[14];

    // Twenty-fourth term
    epsilons[41] += beta_Q * a_eval[ind] * nu * deltas[16];

    // Twenty-fifth term
    epsilons[41] -= 3*beta_Q * a_eval[ind] * beta * gammas[3];

    // Twenty-sixth term
    epsilons[41] -= 2*beta_Q * a_eval[ind] * beta * gammas[4];

    // Twenty-seventh term
    epsilons[41] -= 2*beta_Q * a_eval[ind] * alpha * gammas[5];
    break;
  }
  case 43:
  {
    // First term
    epsilons[42] += epsilons_Q[42] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[42] += (4*epsilons_Q[51] + 3*epsilons_Q[52] +epsilons_Q[53] + 2*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[42] += 18*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[42] += 2*deltas_Q[14] * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[42] += 11*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[42] += 2*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[42] -= gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[5];

    // Eighth term
    epsilons[42] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Ninth term
    epsilons[42] -= 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[5];

    // Tenth term
    epsilons[42] -= 3*gammas_Q[5] * a_eval[ind] * pow(beta, 2) * pow(sigma, 2);

    // Eleventh term
    epsilons[42] += 3*alpha_Q * pow(a_eval[ind], 2) * deltas[17];

    // Twelfth term
    epsilons[42] += beta_Q * a_eval[ind] * sigma * deltas[14];

    // Thirteenth term
    epsilons[42] -= 3*beta_Q * a_eval[ind] * nu * deltas[17];

    // Fourteenth term
    epsilons[42] += beta_Q * a_eval[ind] * beta * gammas[5];
    break;
  }
  case 44:
  {
    // First term
    epsilons[43] += epsilons_Q[43] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[43] -= (4*epsilons_Q[51] + 3*epsilons_Q[52] + 7*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[43] -= 60*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[43] -= (2*deltas_Q[14] - deltas_Q[15]) * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[43] -= 43*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[43] -= 3*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[43] += 3*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[5];

    // Eighth term
    epsilons[43] += 15*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Ninth term
    epsilons[43] += 6*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[5];

    // Tenth term
    epsilons[43] += 14*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);

    // Eleventh term
    epsilons[43] -= 3*alpha_Q * pow(a_eval[ind], 2) * deltas[17];

    // Twelfth term
    epsilons[43] += beta_Q * a_eval[ind] * sigma * deltas[15];

    // Thirteenth term
    epsilons[43] += 3*beta_Q * a_eval[ind] * nu * deltas[17];

    // Fourteenth term
    epsilons[43] -= 4*beta_Q * a_eval[ind] * gammas[5] * beta;
    break;
  }
  case 45:
  {
    // First term
    epsilons[44] += epsilons_Q[44] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[44] += (6*epsilons_Q[51] + epsilons_Q[52] + 6*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[44] += 48*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[44] -= (deltas_Q[14] - deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[44] += 31*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[44] += 4*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[44] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[5];

    // Eighth term
    epsilons[44] += gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[3];

    // Ninth term
    epsilons[44] -= 10*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Tenth term
    epsilons[44] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[5];

    // Eleventh term
    epsilons[44] -= 12*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);

    // Twelfth term
    epsilons[44] += alpha_Q * pow(a_eval[ind], 2) * deltas[17];

    // Thirteenth term
    epsilons[44] += beta_Q * a_eval[ind] * sigma * deltas[16];

    // Fourteenth term
    epsilons[44] -= beta_Q * a_eval[ind] * nu * deltas[17];

    // Fifteenth term
    epsilons[44] += beta_Q * a_eval[ind] * beta * gammas[5];
    break;
  }
  case 46:
  {
    // First term
    epsilons[45] += epsilons_Q[45] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[45] += (epsilons_Q[52] + epsilons_Q[53] - 2*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[45] -= 18*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[45] += deltas_Q[15] * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[45] -= 14*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[45] -= deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[45] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[5];

    // Eighth term
    epsilons[45] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Ninth term
    epsilons[45] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[5];

    // Tenth term
    epsilons[45] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);

    // Eleventh term
    epsilons[45] -= 2*beta_Q * a_eval[ind] * sigma * deltas[14];
    break;
  }
  case 47:
  {
    // First term
    epsilons[46] += epsilons_Q[46] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[46] -= (epsilons_Q[52] - 15*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[46] += 102*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[46] -= 2*deltas_Q[15] * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[46] += 74*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[46] += 2*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[46] -= 4*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[5];

    // Eighth term
    epsilons[46] -= 24*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Ninth term
    epsilons[46] -= 8*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[5];

    // Tenth term
    epsilons[46] -= 24*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);

    // Eleventh term
    epsilons[46] -= 2*beta_Q * a_eval[ind] * sigma * deltas[15];

    // Twelfth term
    epsilons[46] += 4*beta_Q * a_eval[ind] * beta * gammas[5];
    break;
  }
  case 48:
  {
    // First term
    epsilons[47] += epsilons_Q[47] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[47] += (5*epsilons_Q[52] - 12*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[47] -= 69*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[47] -= (2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[47] -= 48*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[47] -= 3*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[47] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[3];

    // Eighth term
    epsilons[47] += 16*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Ninth term
    epsilons[47] += 19*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);

    // Tenth term
    epsilons[47] -= 2*beta_Q * a_eval[ind] * sigma * deltas[16];
    break;
  }
  case 49:
  {
    // First term
    epsilons[48] += epsilons_Q[48] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[48] += (epsilons_Q[53] + 7*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[48] += 42*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[48] -= (deltas_Q[15] + 2*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[48] += 20*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[48] -= 9*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[48] -= gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[5];

    // Eighth term
    epsilons[48] -= 10*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Ninth term
    epsilons[48] -= 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[5];

    // Tenth term
    epsilons[48] -= 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);

    // Eleventh term
    epsilons[48] += 2*beta_Q * a_eval[ind] * beta * gammas[5];
    break;
  }
  case 50:
  {
    // First term
    epsilons[49] += epsilons_Q[49] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[49] += (3*epsilons_Q[53] - 6*epsilons_Q[54]) * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[49] -= 63*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[49] += 2*deltas_Q[15] * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[49] -= 50*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[49] -= 2*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[49] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * gammas[5];

    // Eighth term
    epsilons[49] += 6*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[3];

    // Ninth term
    epsilons[49] += 16*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Tenth term
    epsilons[49] += 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * gammas[5];

    // Eleventh term
    epsilons[49] -= 4*beta_Q * a_eval[ind] * beta * gammas[5];

    // Twelfth term
    epsilons[49] += 13*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);
    break;
  }
  case 51:
  {
    // First term
    epsilons[50] += epsilons_Q[50] * pow(a_eval[ind], 3) * pow(sigma, 6);

    // Second term
    epsilons[50] += 3*epsilons_Q[54] * pow(a_eval[ind], 2) * nu * pow(sigma, 6);

    // Third term
    epsilons[50] += 21*epsilons_Q[55] * a_eval[ind] * pow(nu, 2) * pow(sigma, 6);

    // Fourth term
    epsilons[50] += deltas_Q[16] * pow(a_eval[ind], 2) * pow(sigma, 4) * beta;

    // Fifth term
    epsilons[50] += 20*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * alpha;

    // Sixth term
    epsilons[50] += 5*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 4) * beta;

    // Seventh term
    epsilons[50] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[3];

    // Eighth term
    epsilons[50] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[4];

    // Ninth term
    epsilons[50] -= 3*gammas_Q[5] * a_eval[ind] * pow(sigma, 2) * pow(beta, 2);
    break;
  }
  case 52:
  {
    // First term
    epsilons[51] += epsilons_Q[51] * pow(a_eval[ind], 2) * pow(sigma, 7);

    // Second term
    epsilons[51] += epsilons_Q[55] * a_eval[ind] * nu * pow(sigma, 7);

    // Third term
    epsilons[51] += deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * beta;

    // Fourth term
    epsilons[51] += gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[5];

    // Fifth term
    epsilons[51] += beta_Q * a_eval[ind] * sigma * deltas[17];
    break;
  }
  case 53:
  {
    // First term
    epsilons[52] += epsilons_Q[52] * pow(a_eval[ind], 2) * pow(sigma, 7);

    // Second term
    epsilons[52] += 6*epsilons_Q[55] * a_eval[ind] * nu * pow(sigma, 7);

    // Third term
    epsilons[52] -= 3*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * beta;

    // Fourth term
    epsilons[52] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[5];

    // Fifth term
    epsilons[52] -= 2*beta_Q * a_eval[ind] * sigma * deltas[17];
    break;
  }
  case 54:
  {
    // First term
    epsilons[53] += epsilons_Q[53] * pow(a_eval[ind], 2) * pow(sigma, 7);

    // Second term
    epsilons[53] -= 14*epsilons_Q[55] * a_eval[ind] * nu * pow(sigma, 7);

    // Third term
    epsilons[53] += 6*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * beta;

    // Fourth term
    epsilons[53] += 6*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[5];
    break;
  }
  case 55:
  {
    // First term
    epsilons[54] += epsilons_Q[54] * pow(a_eval[ind], 2) * pow(sigma, 7);

    // Second term
    epsilons[54] += 14*epsilons_Q[55] * a_eval[ind] * nu * pow(sigma, 7);

    // Third term
    epsilons[54] -= 5*deltas_Q[17] * a_eval[ind] * pow(sigma, 5) * beta;

    // Fourth term
    epsilons[54] -= 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * gammas[5];
    break;
  }
  case 56:
  {
    // First term
    epsilons[55] += epsilons_Q[55] * a_eval[ind] * pow(sigma, 8);
    break;
  }
  default:
    throw invalid_argument("Invalid epsilon index");
  }
  return nullptr;
}


// Explicit instantiation for the template class
template class NCoefficients<VectorXd>;
template class NCoefficients<VectorXcd>;
template class NCoefficients<VectorXld>;
template class NCoefficients<VectorXcld>;
template class NCoefficients<VectorXQ>;
template class NCoefficients<VectorXcQ>;
