#include "Tensor.h"

// Helper function which evaluates a term in the polynomial
template <typename Vec>
auto poly_eval(const Vec &a_vec, const Vec &b_vec,
               const int &ind_A, const int &ind_B, Metadata metadata)
{
  using Scalar = typename Vec::value_type; // Extract the type (RealT or complex<RealT>)
  Scalar eval{Scalar{1}};

  // Get the powers of a from the lexicographical index
  vector<int> pow_A{inv_lex_to_pow(metadata.dims[1] - ind_A - 1, metadata.pows[1], metadata.nps[1])};
  int PA = metadata.pows[1];
  // Iterate over a powers and multiply them
  for (int powi{0}; powi < pow_A.size(); ++powi)
  {
    eval *= pow(a_vec[powi], pow_A[powi]);
  }

  // Get the powers of b from the lexicographical index
  vector<int> pow_B{inv_lex_to_pow(metadata.dims[0] - ind_B - 1, metadata.pows[0], metadata.nps[0])};
  int PB = metadata.pows[0];
  // Iterate over b powers and multiply them
  for (int powi{0}; powi < pow_B.size(); ++powi)
  {
    eval *= pow(b_vec[powi], pow_B[powi]);
  }
  return eval;
}

// Helper function which constructs the polynomial terms
string poly_terms(const int &ind_A, const int &ind_B, Metadata metadata)
{
  stringstream poly_ab;

  // Get the powers of a from the lexicographical index
  vector<int> pow_A{inv_lex_to_pow(metadata.dims[1] - ind_A - 1, metadata.pows[1], metadata.nps[1])};
  int PA = metadata.pows[1];
  // Iterate over a powers to construct the a part of the term
  for (int powi{0}; powi < pow_A.size(); ++powi)
  {
    if (pow_A[powi] > 0)
    {
      if (pow_A[powi] == 1)
      {
        poly_ab << "a" << powi + 1;
      }
      else
      {
        poly_ab << "a" << powi + 1 << "^" << pow_A[powi];
      }
      PA -= pow_A[powi];
      poly_ab << (PA != 0 ? "*" : "");
    }
  }
  poly_ab << "*";

  // Get the powers of b from the lexicographical index
  vector<int> pow_B{inv_lex_to_pow(metadata.dims[0] - ind_B - 1, metadata.pows[0], metadata.nps[0])};
  int PB = metadata.pows[0];
  // Iterate over b powers to construct the b part of the term
  for (int powi{0}; powi < pow_B.size(); ++powi)
  {
    if (pow_B[powi] > 0)
    {
      if (pow_B[powi] == 1)
      {
        poly_ab << "b" << powi + 1;
      }
      else
      {
        poly_ab << "b" << powi + 1 << "^" << pow_B[powi];
      }
      PB -= pow_B[powi];
      poly_ab << (PB != 0 ? "*" : "");
    }
  }
  return poly_ab.str();
}

// Constructor for null Tensor
template <typename RealT>
Tensor<RealT>::Tensor(const int &q, const array<int, 2> &coef)
    : metadata(q, coef)
{
  // Resize T based on the dimensions
  T.resize(metadata.dims[0], vector<RealT>(metadata.dims[1]));
}

// Overloaded constructor for null Tensor
// Used to directly specify powers
template <typename RealT>
Tensor<RealT>::Tensor(const int &q, const array<int, 2> &coef, array<int, 2> pows)
    : metadata(q, coef, pows)
{
  // Resize T based on the dimensions
  T.resize(metadata.dims[0], vector<RealT>(metadata.dims[1]));
}

// Constructor for Tensor with one value
template <typename RealT>
Tensor<RealT>::Tensor(const RealT &Q, Powers &powers, const int &q, const array<int, 2> &coef)
    : metadata(q, coef)
{
  // Resize T based on the dimensions
  T.resize(metadata.dims[0], vector<RealT>(metadata.dims[1]));
  // Find the lexicographical index of the powers
  array<const int, 2> lex_inds{powers.get_lex_inds(metadata.dims, metadata.pows)};
  T[lex_inds[0]][lex_inds[1]] = Q; // Add the RealT Q value to the matrix
}

// Overloaded constructor for Tensor with one value
template <typename RealT>
Tensor<RealT>::Tensor(const RealT &Q, Powers &powers, const int &q, const array<int, 2> &coef, const array<int, 2> &pows)
    : metadata(q, coef, pows)
{
  // Resize T based on the dimensions
  T.resize(metadata.dims[0], vector<RealT>(metadata.dims[1]));
  // Find the lexicographical index of the powers
  array<const int, 2> lex_inds{powers.get_lex_inds(metadata.dims, metadata.pows)};
  T[lex_inds[0]][lex_inds[1]] = Q; // Add the RealT Q value to the matrix
}

// Default constructor for Tensor
template <typename RealT>
Tensor<RealT>::Tensor()
{
}

// Copy constructor for Tensor
template <typename RealT>
Tensor<RealT>::Tensor(const Tensor<RealT> &other)
{
  metadata = other.metadata;
  T.resize(metadata.dims[0], vector<RealT>(metadata.dims[1]));
  T = other.T;
}

// Copy assignment operator for Tensor
template <typename RealT>
Tensor<RealT> &Tensor<RealT>::operator=(const Tensor &other)
{
  if (this != &other)
  {
    metadata = other.metadata;
    T = other.T;
  }
  return *this;
}

// Move constructor
template <typename RealT>
Tensor<RealT>::Tensor(Tensor &&other) noexcept
    : metadata(std::move(other.metadata)), T(std::move(other.T))
{
  other.metadata = {}; // Reset to default state
  other.T = {};        // Reset to default state
}

// Move assignment operator
template <typename RealT>
Tensor<RealT> &Tensor<RealT>::operator=(Tensor &&other) noexcept
{
  if (this != &other)
  {
    metadata = std::move(other.metadata); // Move metadata
    T = std::move(other.T);               // Move tensor data
    other.metadata = {};                  // Reset moved-from object
    other.T = {};
  }
  return *this;
}

// Destructor for Tensor
template <typename RealT>
Tensor<RealT>::~Tensor()
{
}

// Method to transform a Tensor to a higher number of cycles or coefficient
// Input: - meta_other: Metadata of the desired new tensor
//        - powers: Additional Powers to add to the transformed tensor
//        - Q: RealT to multiply the tensor with
template <typename RealT>
Tensor<RealT> Tensor<RealT>::transform(const Metadata &meta_other, Powers &powers, const RealT &Q)
{
  if (meta_other.q < metadata.q)
  {
    throw invalid_argument("Error: Must expand to a higher or equal number of cycles q");
  }
  if (meta_other.coef != metadata.coef && powers.iszero())
  {
    throw invalid_argument("Error: Powers class cannot be zero if the coefficient changes");
  }
  // Create the new Tensor
  Tensor<RealT> transf_tens(meta_other.q, meta_other.coef);
  // Iterate over the original Tensor
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      Powers powersi(powers); // Copy the powers
      // Find the original powers and add the new powers to them for a and b
      auto pow_A{inv_lex_to_pow(metadata.dims[1] - j - 1, metadata.pows[1], metadata.nps[1])};
      for (size_t powi{0}; powi < pow_A.size(); ++powi)
      {
        powersi.a_pow[powi] += pow_A[powi];
      }
      auto pow_B{inv_lex_to_pow(metadata.dims[0] - i - 1, metadata.pows[0], metadata.nps[0])};
      for (size_t powi{0}; powi < pow_B.size(); ++powi)
      {
        powersi.b_pow[powi] += pow_B[powi];
      }

      // Find the indices of the powers in the transformed tensor
      int ni{pow_to_lex(powersi.b_pow, transf_tens.metadata.dims[0], transf_tens.metadata.pows[0])};
      int nj{pow_to_lex(powersi.a_pow, transf_tens.metadata.dims[1], transf_tens.metadata.pows[1])};
      // Assign the factor to the correct indices
      transf_tens.T[ni][nj] += Q * T[i][j];
    }
  }
  return transf_tens;
}

// Overloaded method to transform a Tensor to a higher number of cycles
// Input: - meta_other: Metadata of the desired new tensor
//        - powers: Additional Powers to add to the transformed tensor
//        - Q: RealT to multiply the tensor with
template <typename RealT>
Tensor<RealT> Tensor<RealT>::transform(const int &nq)
{
  if (nq < metadata.q)
  {
    throw invalid_argument("Error: Must expand to a higher or equal number of cycles q");
  }
  // Create the new Tensor
  Tensor<RealT> transf_tens(nq, metadata.coef);
  // Iterate over the original Tensor
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      // Find the original powers
      auto pow_A{inv_lex_to_pow(metadata.dims[1] - j - 1, metadata.pows[1], metadata.nps[1])};
      auto pow_B{inv_lex_to_pow(metadata.dims[0] - i - 1, metadata.pows[0], metadata.nps[0])};
      pow_A.resize(transf_tens.metadata.nps[1]);
      pow_B.resize(transf_tens.metadata.nps[0]);
      // Find the indices of the powers in the transformed tensor
      int ni{pow_to_lex(pow_B, transf_tens.metadata.dims[0], transf_tens.metadata.pows[0])};
      int nj{pow_to_lex(pow_A, transf_tens.metadata.dims[1], transf_tens.metadata.pows[1])};

      // Add the constructed Tensor to the new one
      transf_tens.T[ni][nj] += T[i][j];
    }
  }
  return transf_tens;
}

// Method to multiply the Tensor with another one to gain a higher order Tensor
template <typename RealT>
Tensor<RealT> Tensor<RealT>::multiply(const Tensor &other) &
{
  // Create the new Tensor with higher powers
  array<int, 2> pows{{metadata.pows[0] + other.metadata.pows[0], metadata.pows[1] + other.metadata.pows[1]}};
  Tensor<RealT> multi_tens(metadata.q, metadata.coef, pows);
  // Iterate over this tensor
  for (int iB{0}; iB < metadata.dims[0]; ++iB)
  {
    for (int iA{0}; iA < metadata.dims[1]; ++iA)
    {
      // Find the powers of this element
      Powers powersi(inv_lex_to_pow(metadata.dims[1] - iA - 1, metadata.pows[1], metadata.nps[1]),
                     inv_lex_to_pow(metadata.dims[0] - iB - 1, metadata.pows[0], metadata.nps[0]));
      // Iterate over the other tensor
      for (int jB{0}; jB < other.metadata.dims[0]; ++jB)
      {
        for (int jA{0}; jA < other.metadata.dims[1]; ++jA)
        {
          Powers powersj(powersi); // Copy the powers
          // Find the powers of the other element and add them to powersj
          auto pow_A{inv_lex_to_pow(other.metadata.dims[1] - jA - 1, other.metadata.pows[1], other.metadata.nps[1])};
          for (size_t powi{0}; powi < pow_A.size(); ++powi)
          {
            powersj.a_pow[powi] += pow_A[powi];
          }
          auto pow_B{inv_lex_to_pow(other.metadata.dims[0] - jB - 1, other.metadata.pows[0], other.metadata.nps[0])};
          for (size_t powi{0}; powi < pow_B.size(); ++powi)
          {
            powersj.b_pow[powi] += pow_B[powi];
          }

          // Find the indices of the powers in the transformed tensor
          int ni{pow_to_lex(powersj.b_pow, multi_tens.metadata.dims[0], multi_tens.metadata.pows[0])};
          int nj{pow_to_lex(powersj.a_pow, multi_tens.metadata.dims[1], multi_tens.metadata.pows[1])};

          // Add the constructed Tensor to the new one
          multi_tens.T[ni][nj] += T[iB][iA] * other.T[jB][jA];
        }
      }
    }
  }
  return multi_tens;
}

// Method to evaluate the polynomial of the Tensor
// Input: - a_vec, Eigen vector with the dimension of nA
//        - b_vec, Eigen vector with the dimension of nB
template <typename RealT>
template <typename Vec>
typename Vec::value_type Tensor<RealT>::evaluate(const Vec &a_vec, const Vec &b_vec) const
{
  if (a_vec.size() != metadata.nps[1] || b_vec.size() != metadata.nps[0])
  {
    throw invalid_argument("Error: a_vec and b_vec must match the number of parameters");
  }
  using Scalar = typename Vec::value_type; // Extract the type (RealT or complex<RealT>)
  Scalar eval{Scalar{0}};
  // Iterate over the Tensor and add to eval after evaluating term
  if (!(iszero()))
  {
    for (int i{0}; i < metadata.dims[0]; ++i)
    {
      for (int j{0}; j < metadata.dims[1]; ++j)
      {
        // Multiply the Tensor and the evaluation of the term
        eval += T[i][j] * poly_eval(a_vec, b_vec, j, i, metadata);
      }
    }
  }
  return eval;
}

// Method to take a single derivative of a tensor (Used in the derivatives method)
// Input: - dtensor = Tensor to be derivated
//        - param : parameter to be derivated: true = a, false = b
template <typename RealT>
void Tensor<RealT>::derivative(Tensor<RealT> &dtensor, const bool &param, const int &k)
{
  // Find powers of the differentiation Tensor
  array<int, 2> dpows{dtensor.metadata.pows};
  // Cases for differentiation over a or b parameters
  if (param)
  {
    if (dpows[1] > 0)
    {
      dpows[1] -= 1;
    }
  }
  else
  {
    if (dpows[0] > 0)
    {
      dpows[0] -= 1;
    }
  }
  // Construct a temporary tensor with the new powers
  Tensor<RealT> tmptensor(dtensor.metadata.q, dtensor.metadata.coef, dpows);
  Metadata tmpmeta{tmptensor.metadata};
  // Do not differentiate over a or b if the powers of a or b are 0 -> The result is 0
  if ((dtensor.metadata.pows[0] != 0 && !param) || (dtensor.metadata.pows[1] != 0 && param))
  {
    // Iterate over the temporary tensor, pull values from dtensor and multiply them with the proper power
    for (int i{0}; i < tmpmeta.dims[0]; ++i)
    {
      for (int j{0}; j < tmpmeta.dims[1]; ++j)
      {
        // Find the power from the lexicographical index
        vector<int> dpow_A{inv_lex_to_pow(tmpmeta.dims[1] - j - 1, tmpmeta.pows[1], tmpmeta.nps[1])};
        vector<int> dpow_B{inv_lex_to_pow(tmpmeta.dims[0] - i - 1, tmpmeta.pows[0], tmpmeta.nps[0])};
        // Find index A and B according to parameter differentiation and the correct powers
        int indA;
        int indB;
        if (param)
        {
          dpow_A[k] += 1;
          indA = pow_to_lex(dpow_A, dtensor.metadata.dims[1], dtensor.metadata.pows[1]);
          indB = pow_to_lex(dpow_B, dtensor.metadata.dims[0], dtensor.metadata.pows[0]);
          tmptensor.T[i][j] = dtensor.T[indB][indA] * dpow_A[k];
        }
        else
        {
          dpow_B[k] += 1;
          indA = pow_to_lex(dpow_A, dtensor.metadata.dims[1], dtensor.metadata.pows[1]);
          indB = pow_to_lex(dpow_B, dtensor.metadata.dims[0], dtensor.metadata.pows[0]);
          tmptensor.T[i][j] = dtensor.T[indB][indA] * dpow_B[k];
        }
      }
    }
  }
  // Move the tmptensor to dtensor
  dtensor = move(tmptensor);
}

// Method to take multiple derivatives of a tensor
// Input: - a_dev: vector of which ai to derivate
//        - b_dev: vector of which bi to derivate
template <typename RealT>
Tensor<RealT> Tensor<RealT>::derivatives(vector<int> a_dev, vector<int> b_dev)
{

  if (a_dev.empty() && b_dev.empty())
  {
    throw invalid_argument("Error: a_dev and b_dev must not be empty");
  }

  int a_sum{accumulate(a_dev.begin(), a_dev.end(), 0)};
  int b_sum{accumulate(b_dev.begin(), b_dev.end(), 0)};
  if (a_sum == 0 && b_sum == 0)
  {
    throw invalid_argument("Error: a_dev and b_dev must not be zero");
  }

  // Create the new tensor through copy
  Tensor<RealT> dtensor = *this;
  int ka{0};
  // Iterate until the derivative vector a is zero
  while (a_sum != 0)
  {
    // If the vector at index ka is zero, move on
    if (a_dev[ka] == 0)
    {
      ka += 1;
      // Else differentiate
    }
    else
    {
      derivative(dtensor, true, ka);
      a_sum -= 1;
      a_dev[ka] -= 1;
    }
  }
  int kb{0};
  // Iterate until the derivative vector a is zero
  while (b_sum != 0)
  {
    // If the vector at index kb is zero, move on
    if (b_dev[kb] == 0)
    {
      kb += 1;
      // Else differentiate
    }
    else
    {
      derivative(dtensor, false, kb);
      b_sum -= 1;
      b_dev[kb] -= 1;
    }
  }
  return dtensor;
}

// Method to display the Tensor
template <typename RealT>
void Tensor<RealT>::display_tensor()
{
  cout << scientific << setprecision(5);
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      cout << T[i][j] << ", ";
    }
    cout << endl;
  }
}

// Method to display the Tensor to RealT precision
// Specialization for RealT = double
template <>
void Tensor<double>::display_tensor_precise()
{
  cout << scientific << setprecision(20);
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      cout << T[i][j] << ", ";
    }
    cout << endl;
  }
}

// Specialization for RealT = long double
template <>
void Tensor<long double>::display_tensor_precise()
{
  cout << scientific << setprecision(20);
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      cout << T[i][j] << ", ";
    }
    cout << endl;
  }
}

// Method to display the polynomial of the Tensor
template <typename RealT>
void Tensor<RealT>::display_poly()
{
  if (iszero())
  {
    cout << "0" << endl;
  }
  else
  {
    for (int i{0}; i < metadata.dims[0]; ++i)
    {
      for (int j{0}; j < metadata.dims[1]; ++j)
      {
        if (T[i][j] > 0)
        {
          cout << ((i == 0 && j == 0) ? "" : " + ");
        }
        else
        {
          cout << ((i == 0 && j == 0) ? "-" : " - ");
        }
        cout << scientific << setprecision(8) << abs(T[i][j]);
        cout << poly_terms(j, i, metadata)
             << ((i < metadata.dims[0] && j < metadata.dims[1]) ? "" : ")");
      }
    }
    cout << endl;
  }
}

// Method to save the Tensor
// Specialization for RealT = double
template <>
void Tensor<double>::save(const string &dir, const string &label)
{
  // Open the file for the Tensor
  ofstream T_file(dir + label + ".dat");
  T_file << scientific << setprecision(16);
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      T_file << T[i][j];
      if (j != metadata.dims[1] - 1)
      {
        T_file << "\t";
      }
      else
      {
        T_file << endl;
      }
    }
  }
}

// Specialization for RealT = long double
template <>
void Tensor<long double>::save(const string &dir, const string &label)
{
  // Open the file for the Tensor
  ofstream T_file(dir + label + ".dat");
  T_file << scientific << setprecision(20);
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      T_file << T[i][j];
      if (j != metadata.dims[1] - 1)
      {
        T_file << "\t";
      }
      else
      {
        T_file << endl;
      }
    }
  }
}

// Method to load the Tensor
// Specialization for RealT = double
template <>
void Tensor<double>::load(const string &dir, const string &label)
{
  // Open the file for the Tensor
  ifstream T_file(dir + label + ".dat");

  // Check if the file opened correctly
  if (!T_file.is_open())
  {
    throw runtime_error("Error opening file!");
  }

  string line;
  int i{0};
  while (getline(T_file, line))
  {
    stringstream lstream(line);
    string token;
    // Split line by \t
    int j{0};
    while (getline(lstream, token, '\t'))
    {
      // Split token by '/'
      T[i][j] = stod(token);
      j += 1;
    }
    i += 1;
  }
  // Close the Tensor file
  T_file.close();
}

// Specialization for RealT = long double
template <>
void Tensor<long double>::load(const string &dir, const string &label)
{
  // Open the file for the Tensor
  ifstream T_file(dir + label + ".dat");

  // Check if the file opened correctly
  if (!T_file.is_open())
  {
    throw runtime_error("Error opening file!");
  }

  string line;
  int i{0};
  while (getline(T_file, line))
  {
    stringstream lstream(line);
    string token;
    // Split line by \t
    int j{0};
    while (getline(lstream, token, '\t'))
    {
      // Split token by '/'
      T[i][j] = stold(token);
      j += 1;
    }
    i += 1;
  }
  // Close the Tensor file
  T_file.close();
}

// Method to check whether the Tensor is null
template <typename RealT>
bool Tensor<RealT>::iszero() const
{
  for (int i{0}; i < metadata.dims[0]; ++i)
  {
    for (int j{0}; j < metadata.dims[1]; ++j)
    {
      if (T[i][j] != 0.0)
      {
        return false;
      };
    }
  }
  return true;
}

// Overload the += operator for Tensor addition
template <typename RealT>
Tensor<RealT> &Tensor<RealT>::operator+=(const Tensor<RealT> &other)
{
  // Check if dimensions match
  if (metadata.dims != other.metadata.dims)
  {
    throw invalid_argument("Tensor dimensions do not match for addition.");
  }
  // Perform element-wise addition
  for (size_t i = 0; i < T.size(); ++i)
  {
    for (size_t j = 0; j < T[i].size(); ++j)
    {
      T[i][j] = T[i][j] + other.T[i][j];
    }
  }
  return *this;
}

// Overload the -= operator for Tensor addition
template <typename RealT>
Tensor<RealT> &Tensor<RealT>::operator-=(const Tensor<RealT> &other)
{
  // Check if dimensions match
  if (metadata.dims != other.metadata.dims)
  {
    throw invalid_argument("Tensor dimensions do not match for addition.");
  }
  // Perform element-wise subtraction
  for (size_t i = 0; i < T.size(); ++i)
  {
    for (size_t j = 0; j < T[i].size(); ++j)
    {
      T[i][j] = T[i][j] - other.T[i][j];
    }
  }
  return *this;
}

// Overload the * operator for multiplication with RealT
template <typename RealT>
Tensor<RealT> Tensor<RealT>::operator*(const RealT &val)
{
  // Perform element-wise multiplication
  Tensor tens(metadata.q, metadata.coef);
  for (size_t i = 0; i < tens.T.size(); ++i)
  {
    for (size_t j = 0; j < tens.T[i].size(); ++j)
    {
      tens.T[i][j] = val * T[i][j];
    }
  }
  return tens;
}

// Overload the *= operator for multiplication of Tensors and integers
template <typename RealT>
Tensor<RealT> &Tensor<RealT>::operator*=(const int &val)
{
  // Perform element-wise multiplication
  for (size_t i = 0; i < T.size(); ++i)
  {
    for (size_t j = 0; j < T[i].size(); ++j)
    {
      T[i][j] = T[i][j] * val;
    }
  }
  return *this;
}

// Explicit instantiation for the template class
template class Tensor<double>;
template class Tensor<long double>;

// Explicit instantiation for the evaluate method
template double Tensor<double>::evaluate<VectorXd>(const VectorXd &, const VectorXd &) const;
template complex<double> Tensor<double>::evaluate<VectorXcd>(const VectorXcd &, const VectorXcd &) const;

template long double Tensor<long double>::evaluate<VectorXld>(const VectorXld &, const VectorXld &) const;
template complex<long double> Tensor<long double>::evaluate<VectorXcld>(const VectorXcld &, const VectorXcld &) const;
