#include "Minim_helpers.h"

// Helper function to compare two vectors
template <typename Vec>
bool compare_min(const Vec &v1, const Vec &v2, const typename Eigen::NumTraits<typename Vec::Scalar>::Real &tol)
{
  auto norm = (v1 - v2).norm();
  return norm < tol;
}

// Helper to transfrom from the symmetric basis to the standard basis
template <typename Vec>
pair<Vec, Vec> to_standard(const int &q, Vec const &a_vec, Vec const &b_vec)
{
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;
  Vec a_sym(a_vec.size());
  Vec b_sym(b_vec.size());
  if (q % 2 == 0)
  {
    for (size_t i{0}; i < a_vec.size(); ++i)
    {
      if (i == q / 2)
      {
        a_sym[i] = a_vec[0];
      }
      else
      {
        a_sym[i] = a_vec[a_vec.size() - 1 - i] / RealT(2.0);
      }
    }
    for (size_t i{0}; i < b_vec.size(); ++i)
    {
      b_sym[i] = b_vec[b_vec.size() - 1 - i] / RealT(2.0);
    }
  }
  else
  {
    for (size_t i{0}; i < a_vec.size(); ++i)
    {
      a_sym[i] = a_vec[a_vec.size() - 1 - i] / RealT(2.0);
    }
    for (size_t i{0}; i < b_vec.size(); ++i)
    {
      if (i == (q - 1) / 2)
      {
        b_sym[i] = b_vec[0];
      }
      else
      {
        b_sym[i] = b_vec[b_vec.size() - 1 - i] / RealT(2.0);
      }
    }
  }
  pair<Vec, Vec> ab_sym(a_sym, b_sym);
  return ab_sym;
}

// Helper to transfrom from the symmetric basis to the ramp basis
template <typename Vec>
pair<Vec, Vec> to_ramp(const int &q, Vec const &a_vec, Vec const &b_vec)
{
  pair<Vec, Vec> ab_sym{to_standard(q, a_vec, b_vec)};
  const Vec a_sym = ab_sym.first;
  const Vec b_sym = ab_sym.second;
  Vec c_vec(q);
  Vec d_vec(q);
  int back = 0;
  if (q % 2 == 0)
  {
    for (int i{0}; i < q; ++i)
    {
      if (i == 0)
      {
        c_vec[0] = a_sym[0];
        d_vec[0] = b_sym[0] - c_vec[0];
      }
      else if (i < q / 2)
      {
        c_vec[i] = a_sym[i] - d_vec[i - 1];
        d_vec[i] = b_sym[i] - c_vec[i];
      }
      else
      {
        c_vec[i] = a_sym[q / 2 - back] - d_vec[i - 1];
        d_vec[i] = b_sym[(q - 2) / 2 - back] - c_vec[i];
        back += 1;
      }
    }
  }
  else
  {
    for (int i{0}; i < q; ++i)
    {
      if (i == 0)
      {
        c_vec[0] = a_sym[0];
        d_vec[0] = b_sym[0] - c_vec[0];
      }
      else if (i < (q + 1) / 2)
      {
        c_vec[i] = a_sym[i] - d_vec[i - 1];
        d_vec[i] = b_sym[i] - c_vec[i];
      }
      else
      {
        c_vec[i] = a_sym[(q - 1) / 2 - back] - d_vec[i - 1];
        d_vec[i] = b_sym[(q - 3) / 2 - back] - c_vec[i];
        back += 1;
      }
    }
  }
  pair<Vec, Vec> cd_vec(c_vec, d_vec);
  return cd_vec;
}

// Constructor for MinResult
template <typename Vec>
MinResult<Vec>::MinResult(const array<int, 2> &nps)
    : hess(), err2(0.0), conv(0), iter(0),
      a_vec(nps[1]), b_vec(nps[0]), da_vec(nps[1]), db_vec(nps[0])
{
}

// Default constructor for MinResult
template <typename Vec>
MinResult<Vec>::MinResult()
    : hess(), err2(0.0), conv(0), iter(0),
      a_vec(), b_vec(), da_vec(), db_vec()
{
}

// Deconstructor for MinResult
template <typename Vec>
MinResult<Vec>::~MinResult()
{
}

// Method to display the minimization results
template <typename Vec>
void MinResult<Vec>::display()
{
  if (is_same_v<Vec, VectorXd>)
  {
    cout << string(28, '=') << endl;
    cout << "Convergence: " << conv << endl;
    cout << "Number of iterations: " << iter + 1 << endl;
    cout << fixed << std::setprecision(16);
    cout << "a_vec:" << endl
         << a_vec << endl;
    cout << "da_vec:" << endl
         << da_vec << endl;
    cout << "b_vec:" << endl
         << b_vec << endl;
    cout << "db_vec:" << endl
         << db_vec << endl;
    cout << string(28, '=') << endl;
  }
  else if (is_same_v<Vec, VectorXcd>)
  {
    cout << string(52, '=') << endl;
    cout << "Convergence: " << conv << endl;
    cout << "Number of iterations: " << iter + 1 << endl;
    cout << fixed << setprecision(16);
    cout << "a_vec:" << endl;
    for (int i{0}; i < a_vec.size(); ++i)
    {
      cout << real(a_vec[i]);
      if (imag(a_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(a_vec[i]));
      cout << endl;
    }

    cout << "da_vec:" << endl;
    for (int i{0}; i < da_vec.size(); ++i)
    {
      cout << real(da_vec[i]);
      if (imag(da_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(da_vec[i]));
      cout << endl;
    }

    cout << "b_vec:" << endl;
    for (int i{0}; i < b_vec.size(); ++i)
    {
      cout << real(b_vec[i]);
      if (imag(b_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(b_vec[i]));
      cout << endl;
    }

    cout << "db_vec:" << endl;
    for (int i{0}; i < db_vec.size(); ++i)
    {
      cout << real(db_vec[i]);
      if (imag(db_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(db_vec[i]));
      cout << endl;
    }
    cout << string(52, '=') << endl;
  }
  else if (is_same_v<Vec, VectorXld>)
  {
    cout << string(28, '=') << endl;
    cout << "Convergence: " << conv << endl;
    cout << "Number of iterations: " << iter + 1 << endl;
    cout << fixed << setprecision(20);
    cout << "a_vec:" << endl
         << a_vec << endl;
    cout << "da_vec:" << endl
         << da_vec << endl;
    cout << "b_vec:" << endl
         << b_vec << endl;
    cout << "db_vec:" << endl
         << db_vec << endl;
    cout << string(28, '=') << endl;
  }
  else if (is_same_v<Vec, VectorXcld>)
  {
    cout << string(52, '=') << endl;
    cout << "Convergence: " << conv << endl;
    cout << "Number of iterations: " << iter + 1 << endl;
    cout << fixed << setprecision(20);
    cout << "a_vec:" << endl;
    for (int i{0}; i < a_vec.size(); ++i)
    {
      cout << real(a_vec[i]);
      if (imag(a_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(a_vec[i]));
      cout << endl;
    }

    cout << "da_vec:" << endl;
    for (int i{0}; i < da_vec.size(); ++i)
    {
      cout << real(da_vec[i]);
      if (imag(da_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(da_vec[i]));
      cout << endl;
    }

    cout << "b_vec:" << endl;
    for (int i{0}; i < b_vec.size(); ++i)
    {
      cout << real(b_vec[i]);
      if (imag(b_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(b_vec[i]));
      cout << endl;
    }

    cout << "db_vec:" << endl;
    for (int i{0}; i < db_vec.size(); ++i)
    {
      cout << real(db_vec[i]);
      if (imag(db_vec[i]) > 0)
      {
        cout << " + i";
      }
      else
      {
        cout << " - i";
      }
      cout << abs(imag(db_vec[i]));
      cout << endl;
    }
    cout << string(52, '=') << endl;
  }
}

// Explicit template instantiation for compare_min helper
template bool compare_min<VectorXd>(const VectorXd &v1, const VectorXd &v2, const Eigen::NumTraits<typename VectorXd::Scalar>::Real &tol);
template bool compare_min<VectorXcd>(const VectorXcd &v1, const VectorXcd &v2, const Eigen::NumTraits<typename VectorXcd::Scalar>::Real &tol);
template bool compare_min<VectorXld>(const VectorXld &v1, const VectorXld &v2, const Eigen::NumTraits<typename VectorXld::Scalar>::Real &tol);
template bool compare_min<VectorXcld>(const VectorXcld &v1, const VectorXcld &v2, const Eigen::NumTraits<typename VectorXcld::Scalar>::Real &tol);
template bool compare_min<VectorXQ>(const VectorXQ &v1, const VectorXQ &v2, const Eigen::NumTraits<typename VectorXQ::Scalar>::Real &tol);
template bool compare_min<VectorXcQ>(const VectorXcQ &v1, const VectorXcQ &v2, const Eigen::NumTraits<typename VectorXcQ::Scalar>::Real &tol);

// Explicit template instantiation for to_standard helper
template std::pair<VectorXd, VectorXd> to_standard<VectorXd>(const int&, const VectorXd&, const VectorXd&);
template std::pair<VectorXcd, VectorXcd> to_standard<VectorXcd>(const int&, const VectorXcd&, const VectorXcd&);
template std::pair<VectorXld, VectorXld> to_standard<VectorXld>(const int&, const VectorXld&, const VectorXld&);
template std::pair<VectorXcld, VectorXcld> to_standard<VectorXcld>(const int&, const VectorXcld&, const VectorXcld&);
template std::pair<VectorXQ, VectorXQ> to_standard<VectorXQ>(const int&, const VectorXQ&, const VectorXQ&);
template std::pair<VectorXcQ, VectorXcQ> to_standard<VectorXcQ>(const int&, const VectorXcQ&, const VectorXcQ&);

// Explicit template instantiation for to_ramp helper
template std::pair<VectorXd, VectorXd> to_ramp<VectorXd>(const int&, const VectorXd&, const VectorXd&);
template std::pair<VectorXcd, VectorXcd> to_ramp<VectorXcd>(const int&, const VectorXcd&, const VectorXcd&);
template std::pair<VectorXld, VectorXld> to_ramp<VectorXld>(const int&, const VectorXld&, const VectorXld&);
template std::pair<VectorXcld, VectorXcld> to_ramp<VectorXcld>(const int&, const VectorXcld&, const VectorXcld&);
template std::pair<VectorXQ, VectorXQ> to_ramp<VectorXQ>(const int&, const VectorXQ&, const VectorXQ&);
template std::pair<VectorXcQ, VectorXcQ> to_ramp<VectorXcQ>(const int&, const VectorXcQ&, const VectorXcQ&);

// Explicit instantiation for the MinResult class
template class MinResult<VectorXd>;
template class MinResult<VectorXcd>;
template class MinResult<VectorXld>;
template class MinResult<VectorXcld>;
template class MinResult<VectorXQ>;
template class MinResult<VectorXcQ>;
