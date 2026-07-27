#include "Scheme.h"

// Constructor for Scheme
template <typename RealT>
Scheme<RealT>::Scheme(const int &order, const int &no_cycles, const char &verbose)
    : coefs(order, no_cycles, ind_A, ind_B, swtch), n(order), q(no_cycles), verbose(verbose)
{
  // If the number of cycles is even start with ind_A = 1 and the switch bool = true,
  // otherwise ind_B = 1 and switch = false
  (q % 2 == 0) ? ind_A += 1 : ind_B += 1;
  (q % 2 == 0) ? swtch = true : swtch = false;
}

// Load constructor for Scheme
template <typename RealT>
Scheme<RealT>::Scheme(const string &directory, const int &order, const int &no_cycles, const char &verbose)
    : coefs(directory, order, no_cycles, ind_A, ind_B, swtch), n(order), q(no_cycles), verbose(verbose)
{
  if (q % 2 == 0)
  {
    ind_A += no_cycles / 2 + 1;
    ind_B += no_cycles / 2;
    swtch = true;
  }
  else
  {
    ind_A += (no_cycles + 1) / 2;
    ind_B += (no_cycles + 1) / 2;
    swtch = true;
  }
}

// Default constructor for Scheme
template <typename RealT>
Scheme<RealT>::Scheme()
    : coefs(0, 0, ind_A, ind_B, swtch), n(0), q(0), verbose(0)
{
}

// Copy constructor
template <typename RealT>
Scheme<RealT>::Scheme(const Scheme<RealT> &other)
    : coefs(other.coefs), n(other.n), q(other.q), verbose(other.verbose),
      ind_A(other.ind_A), ind_B(other.ind_B), swtch(other.swtch)
{
}

// Copy assignment operator
template <typename RealT>
Scheme<RealT> &Scheme<RealT>::operator=(const Scheme<RealT> &other)
{
  if (this != &other)
  {
    coefs = other.coefs;
    n = other.n;
    q = other.q;
    verbose = other.verbose;
    ind_A = other.ind_A;
    ind_B = other.ind_B;
    swtch = other.swtch;
  }
  return *this;
}

// Move assignment operator
template <typename RealT>
Scheme<RealT> &Scheme<RealT>::operator=(Scheme &&other) noexcept
{
  if (this != &other)
  {
    coefs = std::move(other.coefs);
    n = other.n;
    q = other.q;
    verbose = other.verbose;
    ind_A = other.ind_A;
    ind_B = other.ind_B;
    swtch = other.swtch;
  }
  return *this;
}

// Destructor for Scheme
template <typename RealT>
Scheme<RealT>::~Scheme()
{
}

// Methods which updates the Scheme with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename RealT>
void Scheme<RealT>::step_A()
{
  // Go down the order and take A steps, because higher orders need coefficients of smaller order
  if (n >= 8)
  {
    for (int i{1}; i <= coefs.epsilons.size(); ++i)
    {
      coefs.epsilons_step(i);
      if (verbose == 2)
      {
        cout << "step A: epsilon" << i << endl;
        coefs.epsilons[i - 1].display_tensor();
      }
    }
  }
  if (n >= 6)
  {
    for (int i{1}; i <= coefs.deltas.size(); ++i)
    {
      coefs.deltas_step(i);
      if (verbose == 2)
      {
        cout << "step A: delta" << i << endl;
        coefs.deltas[i - 1].display_tensor();
      }
    }
  }
  if (n >= 4)
  {
    for (int i{1}; i <= coefs.gammas.size(); ++i)
    {
      coefs.gammas_step(i);
      if (verbose == 2)
      {
        cout << "step A: gamma" << i << endl;
        coefs.gammas[i - 1].display_tensor();
      }
    }
  }

  coefs.alpha_stepA();
  coefs.beta_stepA();

  if (verbose == 2)
  {
    cout << "step A: alpha" << endl;
    coefs.alpha.display_tensor();
    cout << "step A: beta" << endl;
    coefs.beta.display_tensor();
  }

  // Update the A index
  ind_A += 1;
  // Switch to step B
  swtch = true;
}

// Method which updates the Scheme with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename RealT>
void Scheme<RealT>::step_B()
{
  // Go down the order and take B steps, because higher orders need coefficients of smaller order
  if (n >= 8)
  {
    // Reverse the coefficients of order 7
    reverse(coefs.epsilons.begin(), coefs.epsilons.end());
    for (int i{1}; i <= coefs.epsilons.size(); ++i)
    {
      coefs.epsilons_step(i);
      if (verbose == 2)
      {
        cout << "step B: epsilon" << i << endl;
        coefs.epsilons[i - 1].display_tensor();
      }
    }
    reverse(coefs.epsilons.begin(), coefs.epsilons.end());
  }
  if (n >= 6)
  {
    // Reverse the coefficients of order 7
    reverse(coefs.deltas.begin(), coefs.deltas.end());
    for (int i{1}; i <= coefs.deltas.size(); ++i)
    {
      coefs.deltas_step(i);
      if (verbose == 2)
      {
        cout << "step B: delta" << i << endl;
        coefs.deltas[i - 1].display_tensor();
      }
    }
    reverse(coefs.deltas.begin(), coefs.deltas.end());
  }
  if (n >= 4)
  {
    // Reverse the coefficients of order 5
    reverse(coefs.gammas.begin(), coefs.gammas.end());
    for (int i{1}; i <= coefs.gammas.size(); ++i)
    {
      coefs.gammas_step(i);
      if (verbose == 2)
      {
        cout << "step B: gamma" << i << endl;
        coefs.gammas[i - 1].display_tensor();
      }
    }
    reverse(coefs.gammas.begin(), coefs.gammas.end());
  }

  coefs.alpha_stepB();
  coefs.beta_stepB();

  if (verbose == 2)
  {
    cout << "step B: alpha" << endl;
    coefs.alpha.display_tensor();
    cout << "step B: beta" << endl;
    coefs.beta.display_tensor();
  }

  // Update the B index
  ind_B += 1;
  // Switch to step A
  swtch = false;
}

// Method which iterates over the number of cycles
template <typename RealT>
void Scheme<RealT>::iterate(optional<int> total_steps_opt)
{
  // Determine the number of iterations if total_steps is not None
  int total_steps = total_steps_opt.value_or((q % 2 == 0) ? q / 2 : (q + 1) / 2);
  for (int i{0}; i < total_steps; ++i)
  {
    // For even q, start with step_B
    if (swtch)
    {
      if (verbose > 0)
      {
        cout << "Cycle " << ind_A + ind_B << " / " << q << endl;
      }
      step_B();
      if (verbose > 0)
      {
        cout << "Cycle " << ind_A + ind_B << " / " << q << endl;
      }
      step_A();
    }
    else
    { // For odd q, start with step_A
      if (verbose > 0)
      {
        cout << "Cycle " << ind_A + ind_B << " / " << q << endl;
      }
      step_A();
      if (i < total_steps - 1)
      { // Only call step_B if there's another step
        if (verbose > 0)
        {
          cout << "Cycle " << ind_A + ind_B << " / " << q << endl;
        }
        step_B();
      }
    }
  }
}

// Method to transform the Scheme into a higher number of cycles
template <typename RealT>
void Scheme<RealT>::transform(const int &nq)
{

  if (q < 7 && n < 4)
  {
    runtime_error("Cannot transform to higher q if the n < 4\nIncrese the desired order");
  }
  else if (q < 17 && n < 6)
  {
    runtime_error("Cannot transform to higher q if the n < 6\nIncrese the desired order");
  }

  coefs.alpha = coefs.alpha.transform(nq);
  coefs.beta = coefs.beta.transform(nq);
  if (n >= 4)
  {
    for (int i{1}; i <= coefs.gammas.size(); ++i)
    {
      coefs.gammas[i - 1] = coefs.gammas[i - 1].transform(nq);
    }
  }
  if (n >= 6)
  {
    for (int i{1}; i <= coefs.deltas.size(); ++i)
    {
      coefs.deltas[i - 1] = coefs.deltas[i - 1].transform(nq);
    }
  }
}

// Method to save the computed Scheme
template <typename RealT>
void Scheme<RealT>::save(const string &directory)
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
  if (fs::exists(dir))
  {
    std::cout << "Directory already exists: " << dir << std::endl;
  }
  else
  {
    // Create the directory
    if (fs::create_directory(dir))
    {
      std::cout << "Directory created: " << dir << std::endl;
    }
    else
    {
      std::cerr << "Failed to create directory: " << dir << std::endl;
    }
  }

  coefs.alpha.save(dir, "alpha");
  coefs.beta.save(dir, "beta");

  if (n >= 4)
  {
    for (int i{1}; i <= coefs.gammas.size(); ++i)
    {
      coefs.gammas[i - 1].save(dir, "gamma" + to_string(i));
    }
  }
  if (n >= 6)
  {
    for (int i{1}; i <= coefs.deltas.size(); ++i)
    {
      coefs.deltas[i - 1].save(dir, "delta" + to_string(i));
    }
  }
}

// Method, which calculates the squared error for the desired order
// Specialization for RealT = double, Vec = VectorXd
template <>
template <>
double Scheme<double>::err2(const int &order, const VectorXd &a_eval, const VectorXd &b_eval)
{
  double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += pow(coefs.alpha.evaluate(a_eval, b_eval), 2);
    err2 += pow(coefs.beta.evaluate(a_eval, b_eval), 2);
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += pow(coefs.gammas[i].evaluate(a_eval, b_eval), 2);
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += pow(coefs.deltas[i].evaluate(a_eval, b_eval), 2);
    }
    break;
  }
  }
  return err2;
}

// Specialization for RealT = double, Vec = VectorXcd
template <>
template <>
double Scheme<double>::err2(const int &order, const VectorXcd &a_eval, const VectorXcd &b_eval)
{
  double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += pow(abs(coefs.alpha.evaluate(a_eval, b_eval)), 2);
    err2 += pow(abs(coefs.beta.evaluate(a_eval, b_eval)), 2);
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += pow(abs(coefs.gammas[i].evaluate(a_eval, b_eval)), 2);
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += pow(abs(coefs.deltas[i].evaluate(a_eval, b_eval)), 2);
    }
    break;
  }
  }
  return err2;
}

// Specialization for RealT = long double, Vec = VectorXld
template <>
template <>
long double Scheme<long double>::err2(const int &order, const VectorXld &a_eval, const VectorXld &b_eval)
{
  long double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += pow(coefs.alpha.evaluate(a_eval, b_eval), 2);
    err2 += pow(coefs.beta.evaluate(a_eval, b_eval), 2);
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += pow(coefs.gammas[i].evaluate(a_eval, b_eval), 2);
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += pow(coefs.deltas[i].evaluate(a_eval, b_eval), 2);
    }
    break;
  }
  }
  return err2;
}

// Specialization for RealT = long double, Vec = VectorXcld
template <>
template <>
long double Scheme<long double>::err2(const int &order, const VectorXcld &a_eval, const VectorXcld &b_eval)
{
  long double err2{0};
  switch (order)
  {
  case 2:
  {
    err2 += pow(abs(coefs.alpha.evaluate(a_eval, b_eval)), 2);
    err2 += pow(abs(coefs.beta.evaluate(a_eval, b_eval)), 2);
    break;
  }
  case 4:
  {
    for (size_t i{0}; i < coefs.gammas.size(); ++i)
    {
      err2 += pow(abs(coefs.gammas[i].evaluate(a_eval, b_eval)), 2);
    }
    break;
  }
  case 6:
  {
    for (size_t i{0}; i < coefs.deltas.size(); ++i)
    {
      err2 += pow(abs(coefs.deltas[i].evaluate(a_eval, b_eval)), 2);
    }
    break;
  }
  }
  return err2;
}

// Method, which calculates the squared error for the desired order
// Specialization for RealT = double, Vec = VectorXd
template <typename RealT>
template <typename Vec>
RealT Scheme<RealT>::eff(const int &order, const Vec &a_eval, const Vec &b_eval)
{
  return 1 / (pow(q, order) * sqrt(err2(order, a_eval, b_eval)));
}

// Method, which displays squared errors for errors of different order (if known)
// Specialization for RealT = double, Vec = VectorXd
template <>
template <>
void Scheme<double>::display_err2(const VectorXd &a_eval, const VectorXd &b_eval)
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2, a_eval, b_eval) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4, a_eval, b_eval) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6, a_eval, b_eval) << endl;
  }
}

// Specialization for RealT = double, Vec = VectorXcd
template <>
template <>
void Scheme<double>::display_err2(const VectorXcd &a_eval, const VectorXcd &b_eval)
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2, a_eval, b_eval) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4, a_eval, b_eval) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6, a_eval, b_eval) << endl;
  }
}

// Specialization for RealT = long double, Vec = VectorXld
template <>
template <>
void Scheme<long double>::display_err2(const VectorXld &a_eval, const VectorXld &b_eval)
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2, a_eval, b_eval) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4, a_eval, b_eval) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6, a_eval, b_eval) << endl;
  }
}

// Specialization for RealT = long double, Vec = VectorXcld
template <>
template <>
void Scheme<long double>::display_err2(const VectorXcld &a_eval, const VectorXcld &b_eval)
{
  if (n == 2)
  {
    cout << "Squared error value: ";
  }
  else
  {
    cout << "Squared error values: ";
  }
  cout << scientific << setprecision(10) << "Err2 = " << err2(2, a_eval, b_eval) << endl;
  if (n >= 4)
  {
    cout << string(22, ' ') << setprecision(10) << "Err4 = " << err2(4, a_eval, b_eval) << endl;
  }
  if (n >= 6)
  {
    cout << string(22, ' ') << setprecision(10) << "Err6 = " << err2(6, a_eval, b_eval) << endl;
  }
}

// Method, which displays squared errors for errors of different order (if known)
// Specialization for RealT = double, Vec = VectorXd
template <>
template <>
void Scheme<double>::display_eff(const VectorXd &a_eval, const VectorXd &b_eval)
{
  double eff2 = 1 / (pow(q, 2) * sqrt(err2(2, a_eval, b_eval)));
  double eff4;
  double eff6;

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

  if (n >= 4)
  {
    eff4 = 1 / (pow(q, 4) * sqrt(err2(4, a_eval, b_eval)));
    if (eff4 > 1e4)
    {
      cout << scientific << setprecision(10);
    }
    else
    {
      cout << fixed << setprecision(10);
    }
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (n >= 6)
  {
    eff6 = 1 / (pow(q, 6) * sqrt(err2(6, a_eval, b_eval)));
    if (eff6 > 1e4)
    {
      cout << scientific << setprecision(10);
    }
    else
    {
      cout << fixed << setprecision(10);
    }
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }
}

// Specialization for RealT = double, Vec = VectorXcd
template <>
template <>
void Scheme<double>::display_eff(const VectorXcd &a_eval, const VectorXcd &b_eval)
{
  double eff2 = 1 / (pow(q, 2) * sqrt(err2(2, a_eval, b_eval)));
  double eff4;
  double eff6;

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

  if (n >= 4)
  {
    eff4 = 1 / (pow(q, 4) * sqrt(err2(4, a_eval, b_eval)));
    if (eff4 > 1e4)
    {
      cout << scientific << setprecision(10);
    }
    else
    {
      cout << fixed << setprecision(10);
    }
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (n >= 6)
  {
    eff6 = 1 / (pow(q, 6) * sqrt(err2(6, a_eval, b_eval)));
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
  }
}

// Specialization for RealT = long double, Vec = VectorXld
template <>
template <>
void Scheme<long double>::display_eff(const VectorXld &a_eval, const VectorXld &b_eval)
{
  long double eff2 = 1 / (pow(q, 2) * sqrt(err2(2, a_eval, b_eval)));
  long double eff4;
  long double eff6;

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

  if (n >= 4)
  {
    eff4 = 1 / (pow(q, 4) * sqrt(err2(4, a_eval, b_eval)));
    if (eff4 > 1e4)
    {
      cout << scientific << setprecision(10);
    }
    else
    {
      cout << fixed << setprecision(10);
    }
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (n >= 6)
  {
    eff6 = 1 / (pow(q, 6) * sqrt(err2(6, a_eval, b_eval)));
    if (eff6 > 1e4)
    {
      cout << scientific << setprecision(10);
    }
    else
    {
      cout << fixed << setprecision(10);
    }
    cout << string(19, ' ') << "Eff6 = " << eff6 << endl;
  }
}

// Specialization for RealT = long double, Vec = VectorXcld
template <>
template <>
void Scheme<long double>::display_eff(const VectorXcld &a_eval, const VectorXcld &b_eval)
{
  long double eff2 = 1 / (pow(q, 2) * sqrt(err2(2, a_eval, b_eval)));
  long double eff4;
  long double eff6;

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

  if (n >= 4)
  {
    eff4 = 1 / (pow(q, 4) * sqrt(err2(4, a_eval, b_eval)));
    if (eff4 > 1e4)
    {
      cout << scientific << setprecision(10);
    }
    else
    {
      cout << fixed << setprecision(10);
    }
    cout << string(19, ' ') << "Eff4 = " << eff4 << endl;
  }

  if (n >= 6)
  {
    eff6 = 1 / (pow(q, 6) * sqrt(err2(6, a_eval, b_eval)));
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
  }
}

// Explicit instantiation of eff for required types
template double Scheme<double>::eff<VectorXd>(const int &, const VectorXd &, const VectorXd &);
template double Scheme<double>::eff<VectorXcd>(const int &, const VectorXcd &, const VectorXcd &);
template long double Scheme<long double>::eff<VectorXld>(const int &, const VectorXld &, const VectorXld &);
template long double Scheme<long double>::eff<VectorXcld>(const int &, const VectorXcld &, const VectorXcld &);

// Explicit instantiation for the template class
template class Scheme<double>;
template class Scheme<long double>;
