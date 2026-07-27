#include "NScheme.h"

// Constructor for NScheme
template <typename Vec>
NScheme<Vec>::NScheme(const int &order, const int &no_cycles,
                      const Vec &a_eval, const Vec &b_eval, const char &verbose)
    : coefs(a_eval, b_eval), n(order), q(no_cycles), a_eval(a_eval), b_eval(b_eval), verbose(verbose)
{
  // If the number of cycles is even start with ind_A = 1 and add a1 to nu,
  // otherwise ind_B = 1 and add b1 to sigma
  if (q % 2 == 0)
  {
    ind_A += 1;
    coefs.nu += a_eval[0];
  }
  else
  {
    ind_B += 1;
    coefs.sigma += b_eval[0];
  }
}

// Destructor for NScheme
template <typename Vec>
NScheme<Vec>::~NScheme()
{
}

// Methods which updates the NScheme with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename Vec>
void NScheme<Vec>::step_A()
{
  // Go down the order and take A steps, because higher orders need coefficients of smaller order
  if (n >= 8)
  {
    for (int i{1}; i <= coefs.epsilons.size(); ++i)
    {
      coefs.epsilons_step(i, ind_A);
      if (verbose == 2)
      {
        cout << "step A: epsilon" << i << ": " << coefs.epsilons[i - 1] << endl;
      }
    }
  }
  if (n >= 6)
  {
    for (int i{1}; i <= coefs.deltas.size(); ++i)
    {
      coefs.deltas_step(i, ind_A);
      if (verbose == 2)
      {
        cout << "step A: delta" << i << ": " << coefs.deltas[i - 1] << endl;
      }
    }
  }
  if (n >= 4)
  {
    for (int i{1}; i <= coefs.gammas.size(); ++i)
    {
      coefs.gammas_step(i, ind_A);
      if (verbose == 2)
      {
        cout << "step A: gamma" << i << ": " << coefs.gammas[i - 1] << endl;
      }
    }
  }

  coefs.alpha_stepA(ind_A);
  coefs.beta_stepA(ind_A);

  if (verbose == 2)
  {
    cout << "step A: alpha: " << coefs.alpha << endl;
    cout << "step A: beta: " << coefs.beta << endl;
  }

  // Step for the nu coefficent
  coefs.nu += a_eval[ind_A];

  // Update the A index
  ind_A += 1;
}

// Method which updates the NScheme with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename Vec>
void NScheme<Vec>::step_B()
{
  // Go down the order and take B steps, because higher orders need coefficients of smaller order
  if (n >= 8)
  {
    // Reverse the coefficients
    reverse(coefs.epsilons.begin(), coefs.epsilons.end());
    reverse(coefs.deltas.begin(), coefs.deltas.end());
    reverse(coefs.gammas.begin(), coefs.gammas.end());
    swap(coefs.alpha, coefs.beta);
    swap(coefs.nu, coefs.sigma);
    swap(coefs.a_eval, coefs.b_eval);
    for (int i{1}; i <= coefs.epsilons.size(); ++i)
    {
      coefs.epsilons_step(i, ind_B);
      if (verbose == 2)
      {
        cout << "step B: epsilon" << i << ": " << coefs.epsilons[i - 1] << endl;
      }
    }
    // Un-reverse the coefficients
    swap(coefs.a_eval, coefs.b_eval);
    swap(coefs.nu, coefs.sigma);
    swap(coefs.alpha, coefs.beta);
    reverse(coefs.gammas.begin(), coefs.gammas.end());
    reverse(coefs.deltas.begin(), coefs.deltas.end());
    reverse(coefs.epsilons.begin(), coefs.epsilons.end());
  }
  if (n >= 6)
  {
    // Reverse the coefficients
    reverse(coefs.deltas.begin(), coefs.deltas.end());
    reverse(coefs.gammas.begin(), coefs.gammas.end());
    swap(coefs.alpha, coefs.beta);
    swap(coefs.nu, coefs.sigma);
    swap(coefs.a_eval, coefs.b_eval);
    for (int i{1}; i <= coefs.deltas.size(); ++i)
    {
      coefs.deltas_step(i, ind_B);
      if (verbose == 2)
      {
        cout << "step B: delta" << i << ": " << coefs.deltas[i - 1] << endl;
      }
    }
    // Un-reverse the coefficients
    swap(coefs.a_eval, coefs.b_eval);
    swap(coefs.nu, coefs.sigma);
    swap(coefs.alpha, coefs.beta);
    reverse(coefs.gammas.begin(), coefs.gammas.end());
    reverse(coefs.deltas.begin(), coefs.deltas.end());
  }
  if (n >= 4)
  {
    // Reverse the coefficients
    reverse(coefs.gammas.begin(), coefs.gammas.end());
    swap(coefs.alpha, coefs.beta);
    swap(coefs.nu, coefs.sigma);
    swap(coefs.a_eval, coefs.b_eval);
    for (int i{1}; i <= coefs.gammas.size(); ++i)
    {
      coefs.gammas_step(i, ind_B);
      if (verbose == 2)
      {
        cout << "step B: gamma" << i << ": " << coefs.gammas[i - 1] << endl;
      }
    }
    // Un-reverse the coefficients
    swap(coefs.a_eval, coefs.b_eval);
    swap(coefs.nu, coefs.sigma);
    swap(coefs.alpha, coefs.beta);
    reverse(coefs.gammas.begin(), coefs.gammas.end());
  }

  coefs.alpha_stepB(ind_B);
  coefs.beta_stepB(ind_B);

  if (verbose == 2)
  {
    cout << "step B: alpha: " << coefs.alpha << endl;
    cout << "step B: beta: " << coefs.beta << endl;
  }

  // Step for the nu coefficent
  coefs.sigma += b_eval[ind_B];

  // Update the B index
  ind_B += 1;
}

// Method which iterates over the number of cycles
template <typename Vec>
void NScheme<Vec>::iterate()
{
  // Determine the number of iterations
  int total_steps = (q % 2 == 0) ? q / 2 : (q + 1) / 2;
  for (int i = 0; i < total_steps; ++i)
  {
    // For even q, start with step_B
    if (q % 2 == 0)
    {
      step_B();
      step_A();
    }
    else
    { // For odd q, start with step_A
      step_A();
      if (i < total_steps - 1)
      { // Only call step_B if there's another step
        step_B();
      }
    }
  }
}

// Method, which resets the scheme, introduces new evaluation vectors a
// and b and iterates the scheme with them
template <typename Vec>
void NScheme<Vec>::reiterate(const Vec &a_vec, const Vec &b_vec)
{
  a_eval = a_vec;
  b_eval = b_vec;
  coefs = NCoefficients(a_vec, b_vec);
  if (q % 2 == 0)
  {
    ind_A = 1;
    ind_B = 0;
    coefs.nu += a_eval[0];
  }
  else
  {
    ind_A = 0;
    ind_B = 1;
    coefs.sigma += b_eval[0];
  }
  iterate();
}

// Method, which calculates the squared error for the desired order
// Specialization for Vec = VectorXd
template <>
double NScheme<VectorXd>::err2(const int &order)
{
  double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += coefs.alpha * coefs.alpha;
    err2 += coefs.beta * coefs.beta;
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += coefs.gammas[i] * coefs.gammas[i];
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += coefs.deltas[i] * coefs.deltas[i];
    }
    break;
  }
  case 8:
  {
    for (size_t i{0}; i < coefs.epsilons.size(); ++i)
    {
      err2 += coefs.epsilons[i] * coefs.epsilons[i];
    }
    break;
  }
  }
  return err2;
}

// Specialization for Vec = VectorXcd
template <>
double NScheme<VectorXcd>::err2(const int &order)
{
  double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += abs(coefs.alpha) * abs(coefs.alpha);
    err2 += abs(coefs.beta) * abs(coefs.beta);
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += abs(coefs.gammas[i]) * abs(coefs.gammas[i]);
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += abs(coefs.deltas[i]) * abs(coefs.deltas[i]);
    }
    break;
  }
  case 8:
  {
    for (size_t i{0}; i < coefs.epsilons.size(); ++i)
    {
      err2 += abs(coefs.epsilons[i]) * abs(coefs.epsilons[i]);
    }
    break;
  }
  }
  return err2;
}

// Specialization for Vec = VectorXld
template <>
long double NScheme<VectorXld>::err2(const int &order)
{
  long double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += coefs.alpha * coefs.alpha;
    err2 += coefs.beta * coefs.beta;
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += coefs.gammas[i] * coefs.gammas[i];
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += coefs.deltas[i] * coefs.deltas[i];
    }
    break;
  }
  case 8:
  {
    for (size_t i{0}; i < coefs.epsilons.size(); ++i)
    {
      err2 += coefs.epsilons[i] * coefs.epsilons[i];
    }
    break;
  }
  }
  return err2;
}

// Specialization for Vec = VectorXcld
template <>
long double NScheme<VectorXcld>::err2(const int &order)
{
  long double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += abs(coefs.alpha) * abs(coefs.alpha);
    err2 += abs(coefs.beta) * abs(coefs.beta);
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += abs(coefs.gammas[i]) * abs(coefs.gammas[i]);
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += abs(coefs.deltas[i]) * abs(coefs.deltas[i]);
    }
    break;
  }
  case 8:
  {
    for (size_t i{0}; i < coefs.epsilons.size(); ++i)
    {
      err2 += abs(coefs.epsilons[i]) * abs(coefs.epsilons[i]);
    }
    break;
  }
  }
  return err2;
}

// Specialization for Vec = VectorXQ
template <>
quad NScheme<VectorXQ>::err2(const int &order)
{
  quad err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += coefs.alpha * coefs.alpha;
    err2 += coefs.beta * coefs.beta;
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += coefs.gammas[i] * coefs.gammas[i];
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += coefs.deltas[i] * coefs.deltas[i];
    }
    break;
  }
  case 8:
  {
    for (size_t i{0}; i < coefs.epsilons.size(); ++i)
    {
      err2 += coefs.epsilons[i] * coefs.epsilons[i];
    }
    break;
  }
  }
  return err2;
}

// Specialization for Vec = VectorXcQ
template <>
quad NScheme<VectorXcQ>::err2(const int &order)
{
  quad err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += abs(coefs.alpha) * abs(coefs.alpha);
    err2 += abs(coefs.beta) * abs(coefs.beta);
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += abs(coefs.gammas[i]) * abs(coefs.gammas[i]);
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += abs(coefs.deltas[i]) * abs(coefs.deltas[i]);
    }
    break;
  }
  case 8:
  {
    for (size_t i{0}; i < coefs.epsilons.size(); ++i)
    {
      err2 += abs(coefs.epsilons[i]) * abs(coefs.epsilons[i]);
    }
    break;
  }
  }
  return err2;
}


// Method, which calculates the efficiency for the desired order
template <typename Vec>
typename Eigen::NumTraits<typename Vec::Scalar>::Real NScheme<Vec>::eff(const int &order)
{
  return 1 / (pow(q, order) * sqrt(err2(order)));
}

// Method, which displays squared errors for errors of different order (if known)
// Specialization for Vec = VectorXd
template <>
void NScheme<VectorXd>::display_err2()
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6) << endl;
  }
  if (n >= 8)
  {
    cout << string(22, ' ') << setprecision(10) << "Err8 = " << err2(8) << endl;
  }
}

// Specialization for Vec = VectorXcd
template <>
void NScheme<VectorXcd>::display_err2()
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6) << endl;
  }
  if (n >= 8)
  {
    cout << string(22, ' ') << setprecision(10) << "Err8 = " << err2(8) << endl;
  }
}

// Specialization for Vec = VectorXld
template <>
void NScheme<VectorXld>::display_err2()
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6) << endl;
  }
  if (n >= 8)
  {
    cout << string(22, ' ') << setprecision(10) << "Err8 = " << err2(8) << endl;
  }
}

// Specialization for Vec = VectorXcld
template <>
void NScheme<VectorXcld>::display_err2()
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6) << endl;
  }
  if (n >= 8)
  {
    cout << string(22, ' ') << setprecision(10) << "Err8 = " << err2(8) << endl;
  }
}

// Specialization for Vec = VectorXQ
template <>
void NScheme<VectorXQ>::display_err2()
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6) << endl;
  }
  if (n >= 8)
  {
    cout << string(22, ' ') << setprecision(10) << "Err8 = " << err2(8) << endl;
  }
}

// Specialization for Vec = VectorXcQ
template <>
void NScheme<VectorXcQ>::display_err2()
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6) << endl;
  }
  if (n >= 8)
  {
    cout << string(22, ' ') << setprecision(10) << "Err8 = " << err2(8) << endl;
  }
}

// Method, which displays squared errors for errors of different order (if known)
// Specialization for Vec = VectorXd
template <>
void NScheme<VectorXd>::display_eff()
{
  double eff2 = eff(2);
  double eff4 = eff(4);
  double eff6 = eff(6);
  double eff8 = eff(8);

  if (n == 2)
  {
    cout << "Efficiency value: ";
  }
  else
  {
    cout << "Efficiency values: ";
  }
  if (eff2 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  cout << "Eff2 = " << eff2 << endl;

  if (eff4 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 4)
  {
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (eff6 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 6)
  {
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }

  if (eff8 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 8)
  {
    cout << string(19, ' ') << "Eff8 = " << eff8 << endl;
  }
}

// Specialization for Vec = VectorXcd
template <>
void NScheme<VectorXcd>::display_eff()
{
  double eff2 = eff(2);
  double eff4 = eff(4);
  double eff6 = eff(6);
  double eff8 = eff(8);

  if (n == 2)
  {
    cout << "Efficiency value: ";
  }
  else
  {
    cout << "Efficiency values: ";
  }
  if (eff2 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  cout << "Eff2 = " << eff2 << endl;

  if (eff4 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 4)
  {
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (eff6 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 6)
  {
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }

  if (eff8 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 8)
  {
    cout << string(19, ' ') << "Eff8 = " << eff8 << endl;
  }
}

// Specialization for Vec = VectorXld
template <>
void NScheme<VectorXld>::display_eff()
{
  long double eff2 = eff(2);
  long double eff4 = eff(4);
  long double eff6 = eff(6);
  long double eff8 = eff(8);

  if (n == 2)
  {
    cout << "Efficiency value: ";
  }
  else
  {
    cout << "Efficiency values: ";
  }
  if (eff2 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  cout << "Eff2 = " << eff2 << endl;

  if (eff4 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 4)
  {
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (eff6 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 6)
  {
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }
  
  if (eff8 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 8)
  {
    cout << string(19, ' ') << "Eff8 = " << eff8 << endl;
  }
}

// Specialization for Vec = VectorXcld
template <>
void NScheme<VectorXcld>::display_eff()
{
  long double eff2 = eff(2);
  long double eff4 = eff(4);
  long double eff6 = eff(6);
  long double eff8 = eff(8);

  if (n == 2)
  {
    cout << "Efficiency value: ";
  }
  else
  {
    cout << "Efficiency values: ";
  }
  if (eff2 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  cout << "Eff2 = " << eff2 << endl;

  if (eff4 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 4)
  {
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (eff6 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 6)
  {
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }

  if (eff8 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 8)
  {
    cout << string(19, ' ') << "Eff8 = " << eff8 << endl;
  }
}

// Specialization for Vec = VectorXQ
template <>
void NScheme<VectorXQ>::display_eff()
{
  quad eff2 = eff(2);
  quad eff4 = eff(4);
  quad eff6 = eff(6);
  quad eff8 = eff(8);

  if (n == 2)
  {
    cout << "Efficiency value: ";
  }
  else
  {
    cout << "Efficiency values: ";
  }
  if (eff2 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  cout << "Eff2 = " << eff2 << endl;

  if (eff4 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 4)
  {
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (eff6 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 6)
  {
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }

  if (eff8 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 8)
  {
    cout << string(19, ' ') << "Eff8 = " << eff8 << endl;
  }
}

// Specialization for Vec = VectorXcQ
template <>
void NScheme<VectorXcQ>::display_eff()
{
  quad eff2 = eff(2);
  quad eff4 = eff(4);
  quad eff6 = eff(6);
  quad eff8 = eff(8);

  if (n == 2)
  {
    cout << "Efficiency value: ";
  }
  else
  {
    cout << "Efficiency values: ";
  }
  if (eff2 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  cout << "Eff2 = " << eff2 << endl;

  if (eff4 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 4)
  {
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (eff6 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 6)
  {
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }

  if (eff8 > 1e4)
  {
    cout << scientific << setprecision(10);
  }
  else
  {
    cout << fixed << setprecision(10);
  }
  if (n >= 8)
  {
    cout << string(19, ' ') << "Eff8 = " << eff8 << endl;
  }
}


// Explicit instantiation for the template class
template class NScheme<VectorXd>;
template class NScheme<VectorXcd>;
template class NScheme<VectorXld>;
template class NScheme<VectorXcld>;
template class NScheme<VectorXQ>;
template class NScheme<VectorXcQ>;
