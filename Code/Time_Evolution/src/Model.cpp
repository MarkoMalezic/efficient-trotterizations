#include "Model.h"

// Model constructor
template <typename RealT>
Model<RealT>::Model(const int &L_, const string &label_, const bool &bound_)
    : L(L_), bound(bound_), label(label_)
{
  termsX.resize(L);
  termsY.resize(L);
  termsZ.resize(L);
}

// Model deconstructor
template <typename RealT>
Model<RealT>::~Model()
{
}

// Method to add a term to the Hamiltonian
template <typename RealT>
void Model<RealT>::add_term(const int &ind, const int &dir, const vector<array<Operators::MatrixC2<RealT>, 2>> &ops,
                            std::vector<RealT> coefs, const std::array<int, 2> &sites)
{
  Hterm term(ops, coefs, sites);
  switch (dir)
  {
  case 0:
    termsX[ind] = term;
    break;
  case 1:
    termsY[ind] = term;
    break;
  case 2:
    termsZ[ind] = term;
    break;
  default:
    throw invalid_argument("Invalid direction");
    break;
  }
}

// Default method for building the Hamiltonian (to be overridden by subclasses)
template <typename RealT>
void Model<RealT>::build_model()
{
}

// Method to build the XZ Heisenberg model
template <typename RealT>
void ModelXZ<RealT>::build_model()
{
  array<int, 2> sites{{0, 0}};

  // J * X_{i}.X_{i+1} term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsX{{Operators::PauliX<RealT>(), Operators::PauliX<RealT>()}};
  vector<RealT> coefsX;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0
    if (this->bound && i == this->L - 1)
    {
      coefsX = {0.0};
    }
    // Periodic boundary: J = J
    else
    {
      coefsX = {J};
    }

    this->add_term(i, 0, opsX, coefsX, sites);
  }

  // J * Z_{i}.Z_{i+1} + hmag_{i} * Z term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsZ{{Operators::PauliZ<RealT>(), Operators::PauliZ<RealT>()},
                                                    {Operators::PauliZ<RealT>(), Operators::Id<RealT>()}};
  // Add terms with different magnetic fields on each site
  vector<RealT> coefsZ;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0, hmag = hmag[L-1]
    if (this->bound && i == this->L - 1)
    {
      coefsZ = {0.0, hmag[i]};
    }
    // Periodic boundary: J = J, hmag = hmag[L-1]
    else
    {
      coefsZ = {J, hmag[i]};
    }
    this->add_term(i, 2, opsZ, coefsZ, sites);
  }
}

// Method to build the XXZ Heisenberg model
template <typename RealT>
void ModelXXZ<RealT>::build_model()
{
  array<int, 2> sites{{0, 0}};

  // J * X_{i}.X_{i+1} term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsX{{Operators::PauliX<RealT>(), Operators::PauliX<RealT>()}};
  vector<RealT> coefsX;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0
    if (this->bound && i == this->L - 1)
    {
      coefsX = {0.0};
    }
    // Periodic boundary: J = J
    else
    {
      coefsX = {J};
    }
    this->add_term(i, 0, opsX, coefsX, sites);
  }

  // J * Y_{i}.Y_{i+1} term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsY{{Operators::PauliY<RealT>(), Operators::PauliY<RealT>()}};
  vector<RealT> coefsY;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0
    if (this->bound && i == this->L - 1)
    {
      coefsY = {0.0};
    }
    // Periodic boundary: J = J
    else
    {
      coefsY = {J};
    }
    this->add_term(i, 1, opsY, coefsY, sites);
  }

  // J * Z_{i}.Z_{i+1} + hmag_{i} * Z term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsZ{{Operators::PauliZ<RealT>(), Operators::PauliZ<RealT>()},
                                                    {Operators::PauliZ<RealT>(), Operators::Id<RealT>()}};
  // Add terms with different magnetic fields on each site
  vector<RealT> coefsZ;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0, hmag = hmag[L-1]
    if (this->bound && i == this->L - 1)
    {
      coefsZ = {0.0, hmag[i]};
    }
    // Periodic boundary: J = J, hmag = hmag[L-1]
    else
    {
      coefsZ = {J, hmag[i]};
    }
    this->add_term(i, 2, opsZ, coefsZ, sites);
  }
}

// Method to build the general Heisenberg model
template <typename RealT>
void ModelHeisenberg<RealT>::build_model()
{
  array<int, 2> sites{{0, 0}};

  // Jx * X_{i}.X_{i+1} term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsX{{Operators::PauliX<RealT>(), Operators::PauliX<RealT>()}};
  vector<RealT> coefsX;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0
    if (this->bound && i == this->L - 1)
    {
      coefsX = {0.0};
    }
    // Periodic boundary: J = Jx
    else
    {
      coefsX = {J[0]};
    }
    this->add_term(i, 0, opsX, coefsX, sites);
  }

  // Jy * Y_{i}.Y_{i+1} term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsY{{Operators::PauliY<RealT>(), Operators::PauliY<RealT>()}};
  vector<RealT> coefsY;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0
    if (this->bound && i == this->L - 1)
    {
      coefsY = {0.0};
    }
    // Periodic boundary: J = Jy
    else
    {
      coefsY = {J[1]};
    }
    this->add_term(i, 1, opsY, coefsY, sites);
  }

  // Jz * Z_{i}.Z_{i+1} + hmag_{i} * Z term
  vector<array<Operators::MatrixC2<RealT>, 2>> opsZ{{Operators::PauliZ<RealT>(), Operators::PauliZ<RealT>()},
                                                    {Operators::PauliZ<RealT>(), Operators::Id<RealT>()}};
  // Add terms with different magnetic fields on each site
  vector<RealT> coefsZ;
  for (int i{0}; i < this->L; ++i)
  {
    sites = {i, (i + 1) % this->L}; // Using modulo for periodic boundary
    // Open boundary: J = 0, hmag = hmag[L-1]
    if (this->bound && i == this->L - 1)
    {
      coefsZ = {0.0, hmag[i]};
    }
    // Periodic boundary: J = Jz, hmag = hmag[L-1]
    else
    {
      coefsZ = {J[2], hmag[i]};
    }
    this->add_term(i, 2, opsZ, coefsZ, sites);
  }
}

// Explicit instantiation of Models
template class Model<double>;
template class Model<long double>;
template class ModelXZ<double>;
template class ModelXZ<long double>;
template class ModelXXZ<double>;
template class ModelXXZ<long double>;
template class ModelHeisenberg<double>;
template class ModelHeisenberg<long double>;
