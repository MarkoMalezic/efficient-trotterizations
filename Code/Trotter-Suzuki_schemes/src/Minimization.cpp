#include "Minimization.h"

// Constructor for Minimization
template <typename Vec>
Minimization<Vec>::Minimization(const Vec &a_init, const Vec &b_init,
                                const Scheme<RealT> &scheme, const vector<RealT> &W_vec,
                                const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls, const bool &bderivs)
    : a_vec(a_init), b_vec(b_init), scheme(scheme), n(scheme.n), q(scheme.q),
      n_iter(n_iter), eps(eps), Ls(Ls), derivs(nullptr)
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
  default:
    throw invalid_argument("Invalid order");
  }

  if (bderivs)
  {
    derivs = new vector<vector<Tensor<RealT>>>;
    derivs->resize(m, vector<Tensor<RealT>>(nps[0] + nps[1]));
    for (int j{0}; j < nps[0] + nps[1]; ++j)
    {
      // Construct the differentiation vectors according to the index i
      vector<int> a_dev(nps[1]);
      vector<int> b_dev(nps[0]);
      if (j < nps[1])
      {
        a_dev[j] = 1;
      }
      else
      {
        b_dev[j - nps[1]] = 1;
      }
      // Iterate over the polynomials and take their derivatives
      for (int i{0}; i < m; ++i)
      {
        Tensor<RealT> tens;
        if (i == 0)
        {
          tens = scheme.coefs.alpha;
        }
        else if (i == 1)
        {
          tens = scheme.coefs.beta;
        }
        else if (i < 8)
        {
          tens = scheme.coefs.gammas[i - 2];
        }
        else if (i < 26)
        {
          tens = scheme.coefs.deltas[i - 8];
        }
        (*derivs)[i][j] = Tensor<RealT>(tens.derivatives(a_dev, b_dev));
      }
    }
  }
  else
  {
    derivs = nullptr;
  }
}

// Overloaded constructor for Minimization
template <typename Vec>
Minimization<Vec>::Minimization(const Scheme<RealT> &scheme, const vector<RealT> &W_vec,
                                const int &n_iter, const array<RealT, 4> eps, const array<RealT, 2> Ls, const bool &bderivs)
    : scheme(scheme), n(scheme.n), q(scheme.q),
      n_iter(n_iter), eps(eps), Ls(Ls)
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
  default:
    throw invalid_argument("Invalid order");
  }

  if (bderivs)
  {
    derivs = new vector<vector<Tensor<RealT>>>;
    derivs->resize(m, vector<Tensor<RealT>>(nps[0] + nps[1]));
    for (int j{0}; j < nps[0] + nps[1]; ++j)
    {
      // Construct the differentiation vectors according to the index i
      vector<int> a_dev(nps[1]);
      vector<int> b_dev(nps[0]);
      if (j < nps[1])
      {
        a_dev[j] = 1;
      }
      else
      {
        b_dev[j - nps[1]] = 1;
      }
      // Iterate over the polynomials and take their derivatives
      for (int i{0}; i < m; ++i)
      {
        Tensor<RealT> tens;
        if (i == 0)
        {
          tens = scheme.coefs.alpha;
        }
        else if (i == 1)
        {
          tens = scheme.coefs.beta;
        }
        else if (i < 8)
        {
          tens = scheme.coefs.gammas[i - 2];
        }
        else if (i < 26)
        {
          tens = scheme.coefs.deltas[i - 8];
        }
        (*derivs)[i][j] = Tensor<RealT>(tens.derivatives(a_dev, b_dev));
      }
    }
  }
}

// Destructor for Minimization
template <typename Vec>
Minimization<Vec>::~Minimization()
{
}

// Method to construct and evaluate the vector y
template <typename Vec>
Vec Minimization<Vec>::y_eval(const Vec &a_eval, const Vec &b_eval)
{
  Vec y(m);
  y(0) = scheme.coefs.alpha.evaluate(a_eval, b_eval);
  y(1) = scheme.coefs.beta.evaluate(a_eval, b_eval);
  if (n >= 4)
  {
    for (size_t i{2}; i < scheme.coefs.gammas.size() + 2; ++i)
    {
      y(i) = scheme.coefs.gammas[i - 2].evaluate(a_eval, b_eval);
    }
  }
  if (n >= 6)
  {
    for (size_t i{8}; i < scheme.coefs.deltas.size() + 8; ++i)
    {
      y(i) = scheme.coefs.deltas[i - 8].evaluate(a_eval, b_eval);
    }
  }
  return y;
}

// Method to construct and evaluate the Jacobian
template <typename Vec>
typename Minimization<Vec>::Mat Minimization<Vec>::jacobian()
{
  // Initialize the Jacobian matrix of size m x n, n = nA + nB - 2 (because of the constraint on ai/bi)
  Mat J(m, nps[0] + nps[1] - 2);
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < nps[0] + nps[1] - 2; ++j)
    {
      // If derivatives are already computed, use them
      if (derivs)
      {
        int ind{j};
        int chind{nps[1] - 1};
        if (j >= nps[1] - 1)
        {
          ind += 1;
          chind += nps[0];
        }
        // Use chain rule to calculate the polynomial derivative
        J(i, j) = (*derivs)[i][ind].evaluate(a_vec, b_vec);
        J(i, j) -= (*derivs)[i][chind].evaluate(a_vec, b_vec);
      }
      // Otherwise compute them and calculate the Jacobian element
      else
      {
        Tensor<RealT> tens;
        if (i == 0)
        {
          tens = scheme.coefs.alpha;
        }
        else if (i == 1)
        {
          tens = scheme.coefs.beta;
        }
        else if (i < 8)
        {
          tens = scheme.coefs.gammas[i - 2];
        }
        else if (i < 26)
        {
          tens = scheme.coefs.deltas[i - 8];
        }

        // First term in the chain rule
        vector<int> a_dev(nps[1]);
        vector<int> b_dev(nps[0]);
        if (j < nps[1] - 1)
        {
          a_dev[j] = 1;
        }
        else
        {
          b_dev[j - nps[1] + 1] = 1;
        }
        Tensor<RealT> dtens{tens.derivatives(a_dev, b_dev)};
        J(i, j) = dtens.evaluate(a_vec, b_vec);

        // Second term in the chain rule
        vector<int> a_chdev(nps[1]);
        vector<int> b_chdev(nps[0]);
        if (j < nps[1] - 1)
        {
          a_chdev[nps[1] - 1] = 1;
        }
        else
        {
          b_chdev[nps[0] - 1] = 1;
        }
        Tensor<RealT> dchtens{tens.derivatives(a_chdev, b_chdev)};
        J(i, j) -= dchtens.evaluate(a_vec, b_vec);
      }
    }
  }
  return J;
}

// Method to minimize the polynomial manifold
// Input: - a_vec, vector of ai variables
//        - b_vec, vector of bi variables
//        - lambda, initial value for the damping coefficent
//        - hessian, fixed hessian matrix (if provided it is used throughout the minimization)
//        - verbose, if true the method displays the intermediate results of the minimization
template <typename Vec>
MinResult<Vec> Minimization<Vec>::minimize(RealT lambda, const HessMatrixForVector_t<Vec> *hess, const bool &verbose)
{
  // Declare minimization results
  MinResult<Vec> mini(nps);

  // Declare variables and initialize some of them in the first step
  Vec y{y_eval(a_vec, b_vec)}; // Vector of evaluated polynomials
  Vec yn{};                    // Vector of evaluated polynomials for the new vectors a_vecn, b_vecn

  Mat J{jacobian()};          // The Jacobian matrix
  Mat JtW{compute_JtW(J, W)}; // Compute the transpose/adjoint J * weights
  Vec JtWy{-JtW * y};         // The righthand side of the step h equation
  HessMatrixForVector_t<Vec> JtWJ;              // The approximate Hessian matrix
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

    // Accept step if it is goes downhill (with a slight tolerance)
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
    for (int i{2}; i < 8; ++i)
    {
      mini.err2 += pow(abs(y[i]), 2);
    }
  }
  else if (n == 6)
  {
    for (int i{8}; i < 26; ++i)
    {
      mini.err2 += pow(abs(y[i]), 2);
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
// Input: - a_vec, vector of ai variables
//        - b_vec, vector of bi variables
//        - lambda, initial value for the damping coefficent
//        - verbose, if true the method displays the intermediate results of the minimization
template <typename Vec>
MinResult<Vec> Minimization<Vec>::min_twostep(RealT lambda, const array<RealT, 4> *eps2, const bool &verbose, const bool &freeze)
{
  // Minimize normally in the first step
  RealT lambda0{lambda};
  MinResult<Vec> mini1{minimize(lambda0, nullptr, verbose)};
  // Display results after 1. step if desired
  if (verbose)
  {
    mini1.display();
  }

  // Adjust the weights to impose the constraint
  switch (n)
  {
  case 2:
  {
    throw invalid_argument("Double step minimization does not make sense for order 2 schemes");
    break;
  }
  case 4:
  {
    for (int i{2}; i < 8; ++i)
    {
      W.diagonal()(i) = 0.0;
    }
    break;
  }
  case 6:
  {
    RealT w4 = W.diagonal()(2);
    for (int i{2}; i < 8; ++i)
    {
      W.diagonal()(i) = w4;
    }
    for (int i{8}; i < 26; ++i)
    {
      W.diagonal()(i) = 0.0;
    }
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

// Method to find as many minima of the polynomial manifold
// Input: - N, Number of maximal initial values
//        - lambda, initial value for the damping coefficent
//        - steps, number of minimization steps to take (not the number of iterations with a minimization process)
//        - mu, the average of the initial value distribution
//        - sigma, the variance of the initial value distribution
//        - eps2 - Convergence criteria in the second step if needed
//        - bsort, if true the method sorts the minima by err2 values low->high
template <typename Vec>
pair<array<vector<Vec>, 4>, vector<int>> Minimization<Vec>::find(const int &N, RealT &lambda, const int &steps,
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
      // Change to include more steps if needed
      DiagonalMatrix<RealT, Dynamic> W_copy{W};
      mini = min_twostep(lambda, eps2, false, freeze);
      W = W_copy;
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
    // Find the indices of the sorted array and sort the Chi2 values
    vector<size_t> inds(err2s.size());
    iota(inds.begin(), inds.end(), 0);
    sort(inds.begin(), inds.end(), [&err2s](size_t i1, size_t i2)
         { return err2s[i1] < err2s[i2]; });

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

// Explicit instantiation for the Minimization class
template class Minimization<VectorXd>;
template class Minimization<VectorXcd>;
template class Minimization<VectorXld>;
template class Minimization<VectorXcld>;
