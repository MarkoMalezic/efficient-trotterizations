#include "NMinimization.h"

// Constructor for NMinimization
template <typename Vec>
NMinimization<Vec>::NMinimization(const Vec &a_init, const Vec &b_init, const RealT step,
                                  const int &n, const int &q, const vector<RealT> &W_vec,
                                  const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls)
    : a_vec(a_init), b_vec(b_init), n(n), q(q), step(step), n_iter(n_iter), eps(eps), Ls(Ls)
{
  // Construct the diagonal weigths matrix
  if constexpr (is_same_v<Vec, VectorXd> || is_same_v<Vec, VectorXcd>)
  {
    W = DiagonalMatrix<RealT, Dynamic>(VectorXd::Map(W_vec.data(), W_vec.size()));
  }
  else if constexpr (is_same_v<Vec, VectorXld> || is_same_v<Vec, VectorXcld>)
  {
    W = DiagonalMatrix<RealT, Dynamic>(VectorXld::Map(W_vec.data(), W_vec.size()));
  }
  else if constexpr (is_same_v<Vec, VectorXQ> || is_same_v<Vec, VectorXcQ>)
  {
    W = DiagonalMatrix<RealT, Dynamic>(VectorXQ::Map(W_vec.data(), W_vec.size()));
  }

  // The number of parameters is dependent on even/odd number of cycles
  if (q % 2 == 0)
  {
    nps = {q / 2, q / 2 + 1};
  }
  else
  {
    nps = {(q + 1) / 2, (q + 1) / 2};
  }

  if (a_vec.size() != nps[1] || b_vec.size() != nps[0])
  {
    throw invalid_argument("Size of input vectors does not match expected dimensions.");
  }

  // The number of polynomials we need depends on the given order
  switch (n)
  {
  case 2:
  {
    m = 2;
    break;
  }
  case 4:
  {
    m = 8;
    break;
  }
  case 6:
  {
    m = 26;
    break;
  }
  case 8:
  {
    m = 82;
    break;
  }
  default:
    throw invalid_argument("Invalid order");
  }
}

// Overloaded constructor for NMinimization
template <typename Vec>
NMinimization<Vec>::NMinimization(const RealT step, const int &n, const int &q, const vector<RealT> &W_vec,
                                  const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls)
    : n(n), q(q), step(step), n_iter(n_iter), eps(eps), Ls(Ls)
{
  // Construct the diagonal weigths matrix
  if constexpr (is_same_v<Vec, VectorXd> || is_same_v<Vec, VectorXcd>)
  {
    W = DiagonalMatrix<RealT, Dynamic>(VectorXd::Map(W_vec.data(), W_vec.size()));
  }
  else if constexpr (is_same_v<Vec, VectorXld> || is_same_v<Vec, VectorXcld>)
  {
    W = DiagonalMatrix<RealT, Dynamic>(VectorXld::Map(W_vec.data(), W_vec.size()));
  }
  else if constexpr (is_same_v<Vec, VectorXQ> || is_same_v<Vec, VectorXcQ>)
  {
    W = DiagonalMatrix<RealT, Dynamic>(VectorXQ::Map(W_vec.data(), W_vec.size()));
  }

  // The number of parameters is dependent on even/odd number of cycles
  if (q % 2 == 0)
  {
    nps = {q / 2, q / 2 + 1};
  }
  else
  {
    nps = {(q + 1) / 2, (q + 1) / 2};
  }

  a_vec = Vec(nps[1]);
  b_vec = Vec(nps[0]);

  // The number of polynomials we need depends on the given order
  switch (n)
  {
  case 2:
  {
    m = 2;
    break;
  }
  case 4:
  {
    m = 8;
    break;
  }
  case 6:
  {
    m = 26;
    break;
  }
  case 8:
  {
    m = 82;
    break;
  }
  default:
    throw invalid_argument("Invalid order");
  }
}

// Destructor for NMinimization
template <typename Vec>
NMinimization<Vec>::~NMinimization()
{
}

// Method to construct the vector y
template <typename Vec>
Vec NMinimization<Vec>::y_eval(const Vec &a_eval, const Vec &b_eval)
{
  NScheme scheme(n, q, a_eval, b_eval, 0);
  scheme.iterate();
  Vec y(m);
  y(0) = scheme.coefs.alpha;
  y(1) = scheme.coefs.beta;
  if (n >= 4)
  {
    for (size_t i{0}; i < scheme.coefs.gammas.size(); ++i)
    {
      y(i + 2) = scheme.coefs.gammas[i];
    }
  }
  if (n >= 6)
  {
    for (size_t i{0}; i < scheme.coefs.deltas.size(); ++i)
    {
      y(i + 8) = scheme.coefs.deltas[i];
    }
  }
  if (n >= 8)
  {
    for (size_t i{0}; i < scheme.coefs.epsilons.size(); ++i)
    {
      y(i + 26) = scheme.coefs.epsilons[i];
    }
  }
  return y;
}

// Method to construct the vector y with the additional term (distance from the origin)
template <typename Vec>
Vec NMinimization<Vec>::y_eval_origin(const Vec &a_eval, const Vec &b_eval, const RealT &ratio)
{
  NScheme scheme(n, q, a_eval, b_eval, 0);
  scheme.iterate();
  Vec y(m+1);
  y(0) = scheme.coefs.alpha;
  y(1) = scheme.coefs.beta;
  // Fill up the y_eval vector with the coefficents (multiplied by the 1-ratio factor in the final order)
  if (n == 4)
  {
    for (size_t i{0}; i < scheme.coefs.gammas.size(); ++i)
    {
      y(i + 2) = (1 - ratio) * scheme.coefs.gammas[i];
    }
  }
  else if (n == 6)
  {
    for (size_t i{0}; i < scheme.coefs.gammas.size(); ++i)
    {
      y(i + 2) = scheme.coefs.gammas[i];
    }
    for (size_t i{0}; i < scheme.coefs.deltas.size(); ++i)
    {
      y(i + 8) = (1 - ratio) * scheme.coefs.deltas[i];
    }
  }
  else if (n == 8)
  {
    for (size_t i{0}; i < scheme.coefs.gammas.size(); ++i)
    {
      y(i + 2) = scheme.coefs.gammas[i];
    }
    for (size_t i{0}; i < scheme.coefs.deltas.size(); ++i)
    {
      y(i + 8) = scheme.coefs.deltas[i];
    }
    for (size_t i{0}; i < scheme.coefs.epsilons.size(); ++i)
    {
      y(i + 26) = (1 - ratio) * scheme.coefs.epsilons[i];
    }
  }
  // Add the extra term to the end of the y_eval vector
  pair<Vec, Vec> cd_vec{to_ramp(q, a_eval, b_eval)};
  Scalar x_bar2 = Scalar(0.0);
  for(int i{0}; i < q; i++)
  {
    x_bar2 += pow(cd_vec.first[i] - Scalar(1.0 / (2*q)), 2);
    x_bar2 += pow(cd_vec.second[i] - Scalar(1.0 / (2*q)), 2);
  }
  y(m) = ratio * pow(x_bar2, 0.5);
  return y;
}

// Method to construct the Jacobian
template <typename Vec>
typename NMinimization<Vec>::Mat NMinimization<Vec>::jacobian()
{
  // Firstly calculate the derivatives of aq (second term in the chain rule over ai derivatives)
  Vec aqp{a_vec};
  aqp(nps[1] - 1) += step;
  Vec aqm{a_vec};
  aqm(nps[1] - 1) -= step;
  Vec y_aq{(y_eval(aqp, b_vec) - y_eval(aqm, b_vec)) / (2 * step)};

  // Secondly calculate the derivatives of bq (second term in the chain rule over bi derivatives)
  Vec bqp{b_vec};
  bqp(nps[0] - 1) += step;
  Vec bqm{b_vec};
  bqm(nps[0] - 1) -= step;
  Vec y_bq{(y_eval(a_vec, bqp) - y_eval(a_vec, bqm)) / (2 * step)};

  // Initialize the Jacobian matrix of size m x n, n = nA + nB - 2 (because of the constraint on ai/bi)
  Mat J(m, nps[0] + nps[1] - 2);
  for (int j = 0; j < nps[0] + nps[1] - 2; ++j)
  {
    Vec step_jp;
    Vec step_jm;
    Vec y_j;
    if (j < nps[1] - 1)
    {
      step_jp = a_vec;
      step_jp(j) += step;
      step_jm = a_vec;
      step_jm(j) -= step;
      y_j = (y_eval(step_jp, b_vec) - y_eval(step_jm, b_vec)) / (2 * step);
      for (int i = 0; i < m; ++i)
      {
        J(i, j) = y_j[i] - y_aq[i];
      }
    }
    else
    {
      step_jp = b_vec;
      step_jp(j - nps[1] + 1) += step;
      step_jm = b_vec;
      step_jm(j - nps[1] + 1) -= step;
      y_j = (y_eval(a_vec, step_jp) - y_eval(a_vec, step_jm)) / (2 * step);
      for (int i = 0; i < m; ++i)
      {
        J(i, j) = y_j[i] - y_bq[i];
      }
    }
  }
  return J;
}

// Method to construct the Jacobian with the additional term (distance from the origin)
template <typename Vec>
typename NMinimization<Vec>::Mat NMinimization<Vec>::jacobian_origin(const RealT &ratio)
{
  // Firstly calculate the derivatives of aq (second term in the chain rule over ai derivatives)
  Vec aqp{a_vec};
  aqp(nps[1] - 1) += step;
  Vec aqm{a_vec};
  aqm(nps[1] - 1) -= step;
  Vec y_aq{(y_eval_origin(aqp, b_vec, ratio) - y_eval_origin(aqm, b_vec, ratio)) / (2 * step)};

  // Secondly calculate the derivatives of bq (second term in the chain rule over bi derivatives)
  Vec bqp{b_vec};
  bqp(nps[0] - 1) += step;
  Vec bqm{b_vec};
  bqm(nps[0] - 1) -= step;
  Vec y_bq{(y_eval_origin(a_vec, bqp, ratio) - y_eval_origin(a_vec, bqm, ratio)) / (2 * step)};

  // Initialize the Jacobian matrix of size m x n, n = nA + nB - 2 (because of the constraint on ai/bi)
  Mat J(m+1, nps[0] + nps[1] - 2);
  for (int j = 0; j < nps[0] + nps[1] - 2; ++j)
  {
    Vec step_jp;
    Vec step_jm;
    Vec y_j;
    if (j < nps[1] - 1)
    {
      step_jp = a_vec;
      step_jp(j) += step;
      step_jm = a_vec;
      step_jm(j) -= step;
      y_j = (y_eval_origin(step_jp, b_vec, ratio) - y_eval_origin(step_jm, b_vec, ratio)) / (2 * step);
      for (int i = 0; i < m+1; ++i)
      {
        J(i, j) = y_j[i] - y_aq[i];
      }
    }
    else
    {
      step_jp = b_vec;
      step_jp(j - nps[1] + 1) += step;
      step_jm = b_vec;
      step_jm(j - nps[1] + 1) -= step;
      y_j = (y_eval_origin(a_vec, step_jp, ratio) - y_eval_origin(a_vec, step_jm, ratio)) / (2 * step);
      for (int i = 0; i < m+1; ++i)
      {
        J(i, j) = y_j[i] - y_bq[i];
      }
    }
  }
  return J;
}

// Method to minimize the polynomial manifold
// Input: - lambda, initial value for the damping coefficent
//        - hessian, fixed hessian matrix (if provided it is used throughout the minimization)
//        - verbose, if true the method displays the intermediate results of the minimization
template <typename Vec>
MinResult<Vec> NMinimization<Vec>::minimize(RealT lambda, const HessMatrixForVector_t<Vec> *hess, const bool &verbose)
{
  // Declare minimization results
  MinResult<Vec> mini(nps);

  // Declare variables and initialize some of them in the first step
  Vec y{y_eval(a_vec, b_vec)}; // Vector of evaluated polynomials
  Vec yn{};                    // Vector of evaluated polynomials for the new vectors a_vecn, b_vecn

  Mat J{jacobian()};               // The Jacobian matrix
  Mat JtW{compute_JtW(J, W)};      // Compute the transpose/adjoint J * weights
  Vec JtWy{-JtW * y};              // The righthand side of the step h equation
  HessMatrixForVector_t<Vec> JtWJ; // The approximate Hessian matrix
  if (hess)
  {
    JtWJ = *hess; // Use the provided Hessian if available
  }
  else
  {
    JtWJ = compute_hess<Mat, Vec>(J, JtW); // Compute the approximate Hessian if not provided
  }
  Mat diag{JtWJ.diagonal().asDiagonal()}; // The diagonal elements of the approximate Hessian

  Mat left(nps[0] + nps[1] - 2, nps[0] + nps[1] - 2); // The lefthand side of the step h equation

  Vec h(nps[0] + nps[1] - 2); // The step vector
  Vec hovera{h.size()};       // Vector h / a (needed for the convergence condition)

  // The new vectors a_vecn, b_vecn
  Vec a_vecn(nps[1]);
  Vec b_vecn(nps[0]);

  // The Chi squared values
  RealT chi2{compute_chi2(y, W)}; // Initial Chi2 values
  RealT chi2n;                    // The new Chi2 value, computed from a_vecn, b_vecn

  // The metric used to accept the step or not
  RealT metric;

  // Iterate the minimization until convergence or until the maximal number of iterations is passed
  for (int iter{0}; iter < n_iter; ++iter)
  {
    // Print the values of the current iteration
    if (verbose)
    {
      cout << "iteration: " << iter + 1 << endl;
      cout << fixed << setprecision(12);
      cout << "a_vec:" << endl
           << a_vec << endl;
      cout << "b_vec:" << endl
           << b_vec << endl;
      cout << "Chi2 = " << setprecision(5) << scientific << chi2 << endl;
      cout << string(30, '-') << endl;
    }

    // Check if the gradient has converged
    if (JtWy.array().abs().maxCoeff() < eps[0])
    {
      mini.conv = 1;
      mini.iter = iter + 1;
      break;
    }

    // Evaluate the left-hand side of the step equation
    left = JtWJ + lambda * diag;

    // Solve the step equation
    h = left.ldlt().solve(JtWy);

    // Compute h / a and check for convergence in the variables
    hovera.head(nps[1] - 1) = a_vec.head(nps[1] - 1);
    hovera.tail(nps[0] - 1) = b_vec.head(nps[0] - 1);
    hovera = h.array() / hovera.array();
    if (hovera.array().abs().maxCoeff() < eps[1])
    {
      mini.conv = 2;
      mini.iter = iter + 1;
      break;
    }

    // Add the step to the new vectors a_vecn
    a_vecn.head(nps[1] - 1) = h.head(nps[1] - 1) + a_vec.head(nps[1] - 1);
    a_vecn[nps[1] - 1] = static_cast<Scalar>(1.0) - a_vecn.head(nps[1] - 1).sum();
    // and b_vecn
    b_vecn.head(nps[0] - 1) = h.tail(nps[0] - 1) + b_vec.head(nps[0] - 1);
    b_vecn[nps[0] - 1] = static_cast<Scalar>(1.0) - b_vecn.head(nps[0] - 1).sum();

    // Evaluate the polynomials with the new vectors a_vecn, b_vecn
    yn = y_eval(a_vecn, b_vecn);
    // Compute the new Chi2 value
    chi2n = compute_chi2(yn, W);

    // Check for convergence in the reduced Chi2 function
    if (chi2n / (m - nps[0] - nps[1] + 2) < eps[2] && m > nps[0] + nps[1] - 2)
    {
      mini.conv = 3;
      mini.iter = iter + 1;
      break;
    }

    // Accept the step if the difference is larger than the tolerance (always go downhill with some numerical instability tolerance)
    if (chi2 - chi2n > -eps[3])
    {
      // Calculate the difference in change for a_vec and b_vec and store in MinResult
      mini.da_vec.head(nps[1] - 1) = h.head(nps[1] - 1);
      mini.da_vec[nps[1] - 1] = a_vecn[nps[1] - 1] - a_vec[nps[1] - 1];

      mini.db_vec.head(nps[0] - 1) = h.tail(nps[0] - 1);
      mini.db_vec[nps[0] - 1] = b_vecn[nps[0] - 1] - b_vec[nps[0] - 1];

      // Increase the damping coefficient
      lambda = max(lambda / Ls[0], static_cast<RealT>(1e-12));

      // Replace the evaluated polynomials and the Chi2
      y = move(yn);
      chi2 = move(chi2n);
      // Replace the a_vec, b_vec with the new vectors a_vecn, b_vecn
      a_vec = move(a_vecn);
      b_vec = move(b_vecn);

      // Compute the new Jacobian
      J = jacobian();
      JtW = compute_JtW(J, W);
      JtWy = -JtW * y;
      if (!hess)
      {
        JtWJ = compute_hess<Mat, Vec>(J, JtW);
        diag = JtWJ.diagonal().asDiagonal();
      }
    }
    else
    {
      // Otherwise, decrease the damping coefficient
      lambda = min(lambda * Ls[1], static_cast<RealT>(1e7));
    }
  }

  // Store the minimization results in MinResult mini
  if (mini.conv == 0)
  {
    mini.iter = n_iter;
  }

  if (n == 2)
  {
    mini.err2 = pow(abs(y[0]), 2) + pow(abs(y[1]), 2);
  }
  else if (n == 4)
  {
    for (int i{0}; i < 6; ++i)
    {
      mini.err2 += pow(abs(y[i+2]), 2);
    }
  }
  else if (n == 6)
  {
    for (int i{0}; i < 18; ++i)
    {
      mini.err2 += pow(abs(y[i+8]), 2);
    }
  }
  else if (n == 8)
  {
    for (int i{0}; i < 56; ++i)
    {
      mini.err2 += pow(abs(y[i+26]), 2);
    }
  }

  mini.hess = JtWJ;
  mini.a_vec = a_vec;
  mini.b_vec = b_vec;

  return mini;
}

// Method to minimize the polynomial manifold with the additional term (distance from the origin)
// Input: - ratio, ratio a between the theoretical error at given order and the distance from the mean
//        - lambda, initial value for the damping coefficent
//        - hessian, fixed hessian matrix (if provided it is used throughout the minimization)
//        - verbose, if true the method displays the intermediate results of the minimization
template <typename Vec>
MinResult<Vec> NMinimization<Vec>::minimize_origin(RealT lambda, const RealT ratio, const HessMatrixForVector_t<Vec> *hess, const bool &verbose)
{
  // Declare minimization results
  MinResult<Vec> mini(nps);

  // Declare variables and initialize some of them in the first step
  Vec y{y_eval_origin(a_vec, b_vec, ratio)}; // Vector of evaluated polynomials
  Vec yn{};                    // Vector of evaluated polynomials for the new vectors a_vecn, b_vecn

  Mat J{jacobian_origin(ratio)};               // The Jacobian matrix
  Mat JtW{compute_JtW(J, W)};      // Compute the transpose/adjoint J * weights
  Vec JtWy{-JtW * y};              // The righthand side of the step h equation
  HessMatrixForVector_t<Vec> JtWJ; // The approximate Hessian matrix
  if (hess)
  {
    JtWJ = *hess; // Use the provided Hessian if available
  }
  else
  {
    JtWJ = compute_hess<Mat, Vec>(J, JtW); // Compute the approximate Hessian if not provided
  }
  Mat diag{JtWJ.diagonal().asDiagonal()}; // The diagonal elements of the approximate Hessian

  Mat left(nps[0] + nps[1] - 2, nps[0] + nps[1] - 2); // The lefthand side of the step h equation

  Vec h(nps[0] + nps[1] - 2); // The step vector
  Vec hovera{h.size()};       // Vector h / a (needed for the convergence condition)

  // The new vectors a_vecn, b_vecn
  Vec a_vecn(nps[1]);
  Vec b_vecn(nps[0]);

  // The Chi squared values
  RealT chi2{compute_chi2(y, W)}; // Initial Chi2 values
  RealT chi2n;                    // The new Chi2 value, computed from a_vecn, b_vecn

  // The metric used to accept the step or not
  RealT metric;

  // Iterate the minimization until convergence or until the maximal number of iterations is passed
  for (int iter{0}; iter < n_iter; ++iter)
  {
    // Print the values of the current iteration
    if (verbose)
    {
      cout << "iteration: " << iter + 1 << endl;
      cout << fixed << setprecision(12);
      cout << "a_vec:" << endl
           << a_vec << endl;
      cout << "b_vec:" << endl
           << b_vec << endl;
      cout << "Chi2 = " << setprecision(5) << scientific << chi2 << endl;
      cout << string(30, '-') << endl;
    }

    // Check if the gradient has converged
    if (JtWy.array().abs().maxCoeff() < eps[0])
    {
      mini.conv = 1;
      mini.iter = iter + 1;
      break;
    }

    // cout << "Gradient: \n" << JtWy << endl;

    // Eigen::SelfAdjointEigenSolver<MatrixXd> solver(JtWJ);
    // VectorXd eigs = solver.eigenvalues();
    // std::cout << "Eigenvalues (Hessian): \n" << eigs << std::endl;

    // Evaluate the left-hand side of the step equation
    left = JtWJ + lambda * diag;

    // Solve the step equation
    h = left.ldlt().solve(JtWy);

    // Compute h / a and check for convergence in the variables
    hovera.head(nps[1] - 1) = a_vec.head(nps[1] - 1);
    hovera.tail(nps[0] - 1) = b_vec.head(nps[0] - 1);
    hovera = h.array() / hovera.array();
    if (hovera.array().abs().maxCoeff() < eps[1])
    {
      mini.conv = 2;
      mini.iter = iter + 1;
      break;
    }

    // Add the step to the new vectors a_vecn
    a_vecn.head(nps[1] - 1) = h.head(nps[1] - 1) + a_vec.head(nps[1] - 1);
    a_vecn[nps[1] - 1] = static_cast<Scalar>(1.0) - a_vecn.head(nps[1] - 1).sum();
    // and b_vecn
    b_vecn.head(nps[0] - 1) = h.tail(nps[0] - 1) + b_vec.head(nps[0] - 1);
    b_vecn[nps[0] - 1] = static_cast<Scalar>(1.0) - b_vecn.head(nps[0] - 1).sum();

    // Evaluate the polynomials with the new vectors a_vecn, b_vecn
    yn = y_eval_origin(a_vecn, b_vecn, ratio);
    // Compute the new Chi2 value
    chi2n = compute_chi2(yn, W);

    // Check for convergence in the reduced Chi2 function
    if (chi2n / (m+1 - nps[0] - nps[1] + 2) < eps[2] && m+1 > nps[0] + nps[1] - 2)
    {
      mini.conv = 3;
      mini.iter = iter + 1;
      break;
    }

    // Accept the step if the difference is larger than the tolerance (always go downhill with some numerical instability tolerance)
    if (chi2 - chi2n > -eps[3])
    {
      // Calculate the difference in change for a_vec and b_vec and store in MinResult
      mini.da_vec.head(nps[1] - 1) = h.head(nps[1] - 1);
      mini.da_vec[nps[1] - 1] = a_vecn[nps[1] - 1] - a_vec[nps[1] - 1];

      mini.db_vec.head(nps[0] - 1) = h.tail(nps[0] - 1);
      mini.db_vec[nps[0] - 1] = b_vecn[nps[0] - 1] - b_vec[nps[0] - 1];

      // Increase the damping coefficient
      lambda = max(lambda / Ls[0], static_cast<RealT>(1e-12));

      // Replace the evaluated polynomials and the Chi2
      y = move(yn);
      chi2 = move(chi2n);
      // Replace the a_vec, b_vec with the new vectors a_vecn, b_vecn
      a_vec = move(a_vecn);
      b_vec = move(b_vecn);

      // Compute the new Jacobian
      J = jacobian_origin(ratio);
      JtW = compute_JtW(J, W);
      JtWy = -JtW * y;
      if (!hess)
      {
        JtWJ = compute_hess<Mat, Vec>(J, JtW);
        diag = JtWJ.diagonal().asDiagonal();
      }
    }
    else
    {
      // Otherwise, decrease the damping coefficient
      lambda = min(lambda * Ls[1], static_cast<RealT>(1e7));
    }
  }

  // Store the minimization results in MinResult mini
  if (mini.conv == 0)
  {
    mini.iter = n_iter;
  }

  if (n == 2)
  {
    mini.err2 = pow(abs(y[0]), 2) + pow(abs(y[1]), 2);
  }
  else if (n == 4)
  {
    for (int i{0}; i < 6; ++i)
    {
      mini.err2 += pow(abs(y[i+2]), 2);
    }
  }
  else if (n == 6)
  {
    for (int i{0}; i < 18; ++i)
    {
      mini.err2 += pow(abs(y[i+8]), 2);
    }
  }
  else if (n == 8)
  {
    for (int i{0}; i < 56; ++i)
    {
      mini.err2 += pow(abs(y[i+26]), 2);
    }
  }

  mini.hess = JtWJ;
  mini.a_vec = a_vec;
  mini.b_vec = b_vec;

  return mini;
}


// Method to minimize the polynomial manifold in two steps
// First step: minimize with inital weights to enter the region of the minima
// Second step: Exactly ensure the constraint by adjusting weights and making use of the hessian from the 1. step
// Input: - lambda, initial value for the damping coefficent
//        - verbose, if true the method displays the intermediate results of the minimization
template <typename Vec>
MinResult<Vec> NMinimization<Vec>::min_twostep(RealT lambda, const array<RealT, 4> *eps2, const bool &verbose, const bool &freeze)
{
  // Minimize normally in the first step
  RealT lambda0{lambda};
  MinResult<Vec> mini1{minimize(lambda0, nullptr, verbose)};

  // Display results after 1. step if desired
  if (verbose)
  {
    mini1.display();
  }


  // Lower the order and adjust the weights to impose the constraint
  switch (n)
  {
  case 2:
  {
    throw invalid_argument("Double step minimization does not make sense for order 2 schemes");
    break;
  }
  case 4:
  {
    n = 2;
    m = 2;
    W.diagonal().conservativeResize(2);
    break;
  }
  case 6:
  {
    n = 4;
    m = 8;
    W.diagonal().conservativeResize(8);
    break;
  }
  case 8:
  {
    n = 6;
    m = 26;
    W.diagonal().conservativeResize(26);
    break;
  }
  default:
    throw invalid_argument("Invalid order");
  }

  // Change the convergence criteria if provided
  if (eps2)
  {
    eps = *eps2;
  }

  MinResult<Vec> mini2;
  // Minimize with the new weights and the hessian from the previous step
  if (freeze)
  {
    mini2 = minimize(lambda, &mini1.hess, verbose);
  }
  // Minimize with the new weights without freezing the hessian
  else
  {
    mini2 = minimize(lambda, nullptr, verbose);
  }

  if (verbose)
  {
    mini2.display();
  }

  return mini2;
}

// Method to minimize the polynomial manifold in two steps
// First step: minimize with inital weights to enter the region of the minima
// Second step: Exactly ensure the constraint by adjusting weights and making use of the hessian from the 1. step
// Input: - lambda, initial value for the damping coefficent
//        - verbose, if true the method displays the intermediate results of the minimization
template <typename Vec>
MinResult<Vec> NMinimization<Vec>::min_twostep_origin(RealT lambda, const RealT ratio, const array<RealT, 4> *eps2, const bool &verbose, const bool &freeze)
{
  // Minimize normally in the first step
  RealT lambda0{lambda};
  MinResult<Vec> mini1{minimize_origin(lambda0, ratio, nullptr, verbose)};

  // Display results after 1. step if desired
  if (verbose)
  {
    mini1.display();
  }

  // Lower the order and adjust the weights to impose the constraint
  switch (n)
  {
  case 2:
  {
    throw invalid_argument("Double step minimization does not make sense for order 2 schemes");
    break;
  }
  case 4:
  {
    n = 2;
    m = 2;
    W.diagonal().conservativeResize(2);
    break;
  }
  case 6:
  {
    n = 4;
    m = 8;
    W.diagonal().conservativeResize(8);
    break;
  }
  case 8:
  {
    n = 6;
    m = 26;
    W.diagonal().conservativeResize(26);
    break;
  }
  default:
    throw invalid_argument("Invalid order");
  }

  // Change the convergence criteria if provided
  if (eps2)
  {
    eps = *eps2;
  }

  // Minimize without the origin point and just impose the constraints
  MinResult<Vec> mini2;
  // Minimize with the new weights and the hessian from the previous step
  if (freeze)
  {
    mini2 = minimize(lambda, &mini1.hess, verbose);
  }
  // Minimize with the new weights without freezing the hessian
  else
  {
    mini2 = minimize(lambda, nullptr, verbose);
  }

  if (verbose)
  {
    mini2.display();
  }

  return mini2;
}

// Method to find as many minima of the polynomial manifold
// Input: - N, Number of maximal initial values
//        - lambda, initial value for the damping coefficent
//        - steps, number of minimization steps to take (not the number of iterations with a minimization process)
//        - mu, the average of the initial value distribution
//        - sigma, the variance of the initial value distribution
//        - eps2 - Convergence criteria in the second step if needed
//        - bsort, if true the method sorts the minima by err2 values low->high
template <typename Vec>
pair<array<vector<Vec>, 4>, vector<int>> NMinimization<Vec>::find(const int &N, RealT &lambda, const int &steps,
                                                                  const Scalar &mu, const Scalar &sigma,
                                                                  const array<RealT, 4> *eps2, const bool &bsort, const bool &verbose,
                                                                  const RealT &tol, const bool &freeze)
{
  // Declare the vectors of found minima and their "errors"
  vector<Vec> a_vecs;
  vector<Vec> b_vecs;
  vector<Vec> da_vecs;
  vector<Vec> db_vecs;
  vector<int> min_convs;
  vector<RealT> err2s;

  // Declare the helper vectors to compare the found minima vectors
  Vec v1(nps[0] + nps[1]);
  Vec v2(nps[0] + nps[1]);
  
  // Declare the MinResult
  MinResult<Vec> mini;

  array<Vec, 2> ab_vec;
  int convs{0};
  // Iterate until the maximal initial condition
  for (int i{0}; i < N; ++i)
  {
    cout << "\r" << "Initial condition: " << i + 1 << flush;

    // Sample the initial conditions
    ab_vec = sample<Vec, Scalar>(mu, sigma, nps);
    a_vec = ab_vec[0];
    b_vec = ab_vec[1];

    // Take either 1 minimization step or multiple
    if (steps == 1)
    {
      mini = minimize(lambda);
    }
    else
    {
      // Copy the weights and dimensions to restore them after the two-step minimization
      DiagonalMatrix<RealT, Dynamic> W_copy{W};
      int n_copy{n};
      int m_copy{m};
      mini = min_twostep(lambda, eps2, false, freeze);
      W = W_copy;
      n = n_copy;
      m = m_copy;
    }

    if (verbose)
    {
      cout << endl;
      mini.display();
    }

    // Only continue if the minimization has converged
    if (mini.conv)
    {
      convs += 1;
      // Append the first found minima
      if (a_vecs.empty())
      {
        a_vecs.push_back(mini.a_vec);
        b_vecs.push_back(mini.b_vec);
        da_vecs.push_back(mini.da_vec);
        db_vecs.push_back(mini.db_vec);
        min_convs.push_back(1);
        err2s.push_back(mini.err2);
      }
      else
      {
        // For the next minima check if the new one is already included
        size_t ab_size{a_vecs.size()};
        bool included{false};
        for (size_t j{0}; j < ab_size; ++j)
        {
          v1.head(nps[1]) = a_vecs[j];
          v1.tail(nps[0]) = b_vecs[j];
          v2.head(nps[1]) = mini.a_vec;
          v2.tail(nps[0]) = mini.b_vec;
          // Compare the minima and break the loop if they are the same
          if (compare_min(v1, v2, tol))
          {
            included = true;
            min_convs[j] += 1;
            // Also replace the new minimum with the old one if the Chi2 value is smaller
            if (mini.err2 < err2s[j])
            {
              a_vecs[j] = mini.a_vec;
              b_vecs[j] = mini.b_vec;
              da_vecs[j] = mini.da_vec;
              db_vecs[j] = mini.db_vec;
              err2s[j] = mini.err2;
            }
            break;
          }
        }
        // If the new minimum isn't included, then append it
        if (!included)
        {
          a_vecs.push_back(mini.a_vec);
          b_vecs.push_back(mini.b_vec);
          da_vecs.push_back(mini.da_vec);
          db_vecs.push_back(mini.db_vec);
          min_convs.push_back(1);
          err2s.push_back(mini.err2);
        }
      }
    }
  }
  cout << endl;

  // Sort the vectors by err2 (low -> high) if desired
  if (bsort)
  {
    // Find the leading order error for each minimum
    vector<RealT> lead_error(err2s.size());
    for (size_t i{0}; i < err2s.size(); ++i)
    {
      NScheme scheme(n, q, a_vecs[i], b_vecs[i], 0);
      scheme.iterate();
      lead_error[i] = scheme.err2(n);
    }
    // Find the indices of the sorted array and sort the Chi2 values
    vector<size_t> inds(lead_error.size());
    iota(inds.begin(), inds.end(), 0);
    sort(inds.begin(), inds.end(), [&lead_error](size_t i1, size_t i2)
         { return lead_error[i1] < lead_error[i2]; });

    // Then sort the minima as well
    vector<Vec> sorted_a_vecs(a_vecs.size());
    vector<Vec> sorted_b_vecs(b_vecs.size());
    vector<Vec> sorted_da_vecs(da_vecs.size());
    vector<Vec> sorted_db_vecs(db_vecs.size());
    vector<int> sorted_min_convs(min_convs.size());
    for (size_t i{0}; i < inds.size(); ++i)
    {
      sorted_a_vecs[i] = a_vecs[inds[i]];
      sorted_b_vecs[i] = b_vecs[inds[i]];
      sorted_da_vecs[i] = da_vecs[inds[i]];
      sorted_db_vecs[i] = db_vecs[inds[i]];
      sorted_min_convs[i] = min_convs[inds[i]];
    }

    array<vector<Vec>, 4> sorted_ab_vecs{{sorted_a_vecs, sorted_da_vecs,
                                          sorted_b_vecs, sorted_db_vecs}};
    pair<array<vector<Vec>, 4>, vector<int>> result;
    result.first = sorted_ab_vecs;
    result.second = sorted_min_convs;
    return result;
  }
  else
  {
    array<vector<Vec>, 4> ab_vecs{{a_vecs, da_vecs,
                                   b_vecs, db_vecs}};
    pair<array<vector<Vec>, 4>, vector<int>> result;
    result.first = ab_vecs;
    result.second = min_convs;
    return result;
  }
}


// Explicit instantiation for the NMinimization class
template class NMinimization<VectorXd>;
template class NMinimization<VectorXcd>;
template class NMinimization<VectorXld>;
template class NMinimization<VectorXcld>;
template class NMinimization<VectorXQ>;
template class NMinimization<VectorXcQ>;
