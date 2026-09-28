#include "IO.h"

// Helper function to trim a single line in the input file
string trim(const string &str)
{
  size_t first = str.find_first_not_of(" \t");
  if (first == string::npos)
    return "";
  size_t last = str.find_last_not_of(" \t");
  return str.substr(first, (last - first + 1));
}

// Helper function to parse a boolean value
bool parse_bool(const string &str)
{
  string lower_str = str;
  transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
  if (lower_str == "true" || lower_str == "1")
  {
    return true;
  }
  else if (lower_str == "false" || lower_str == "0")
  {
    return false;
  }
  else
  {
    throw invalid_argument("Invalid boolean value: " + str);
  }
}

// Helper function to parse a real or complex value
// Specialization for double
template <>
double parse_value(const string &str)
{
  istringstream iss(str);
  double value;
  iss >> value;
  return value;
}

// Specialization for complex<double>
template <>
complex<double> parse_value(const string &str)
{
  istringstream iss(str);
  double real, imag;
  char sign, i;
  iss >> real >> sign >> i >> imag;
  ;
  if (sign == '-')
  {
    imag = -imag;
  }
  return complex<double>(real, imag);
}

// Specialization for long double
template <>
long double parse_value(const string &str)
{
  istringstream iss(str);
  long double value;
  iss >> value;
  return value;
}

// Specialization for complex<long double>
template <>
complex<long double> parse_value(const string &str)
{
  istringstream iss(str);
  long double real, imag;
  char sign, i;
  iss >> real >> sign >> i >> imag;
  ;
  if (sign == '-')
  {
    imag = -imag;
  }
  return complex<long double>(real, imag);
}

// Specialization for quad
template <>
quad parse_value(const string &str)
{
  istringstream iss(str);
  quad value;
  iss >> value;
  return value;
}

// Specialization for complex<quad>
template <>
complex<quad> parse_value(const string &str)
{
  istringstream iss(str);
  quad real, imag;
  char sign, i;
  iss >> real >> sign >> i >> imag;
  ;
  if (sign == '-')
  {
    imag = -imag;
  }
  return complex<quad>(real, imag);
}

// Helper function to parse a list to extract vectors and arrays
// Specialization for array<double, 2>
template <>
array<double, 2> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<double, 2> list;
  for (int i = 0; i < 2; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<double, 2>");
    }
    list[i] = stod(trim(token));
  }
  return list;
}

// Specialization for array<long double, 2>
template <>
array<long double, 2> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<long double, 2> list;
  for (int i = 0; i < 2; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<long double, 2>");
    }
    list[i] = stold(trim(token));
  }
  return list;
}

// Specialization for array<quad, 2>
template <>
array<quad, 2> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<quad, 2> list;
  for (int i = 0; i < 2; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<quad, 2>");
    }
    list[i] = stold(trim(token));
  }
  return list;
}

// Specialization for array<double, 4>
template <>
array<double, 4> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<double, 4> list;
  for (int i = 0; i < 4; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<double, 2>");
    }
    list[i] = stod(trim(token));
  }
  return list;
}

// Specialization for array<long double, 4>
template <>
array<long double, 4> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<long double, 4> list;
  for (int i = 0; i < 4; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<long double, 4>");
    }
    list[i] = stold(trim(token));
  }
  return list;
}

// Specialization for array<quad, 4>
template <>
array<quad, 4> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<quad, 4> list;
  for (int i = 0; i < 4; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<quad, 4>");
    }
    list[i] = stold(trim(token));
  }
  return list;
}

// Specialization for vector<double>
template <>
vector<double> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  vector<double> list;
  while (getline(iss, token, ','))
  {
    list.push_back(stod(trim(token)));
  }
  return list;
}

// Specialization for vector<complex<double>>
template <>
vector<complex<double>> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  vector<complex<double>> list;
  while (getline(iss, token, ','))
  {
    list.push_back(parse_value<complex<double>>(trim(token)));
  }
  return list;
}

// Specialization for vector<long double>
template <>
vector<long double> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  vector<long double> list;
  while (getline(iss, token, ','))
  {
    list.push_back(stold(trim(token)));
  }
  return list;
}

// Specialization for vector<complex<long double>>
template <>
vector<complex<long double>> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  vector<complex<long double>> list;
  while (getline(iss, token, ','))
  {
    list.push_back(parse_value<complex<long double>>(trim(token)));
  }
  return list;
}

// Specialization for vector<quad>
template <>
vector<quad> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  vector<quad> list;
  while (getline(iss, token, ','))
  {
    quad value(trim(token));
    list.push_back(value);
  }
  return list;
}

// Specialization for vector<complex<quad>>
template <>
vector<complex<quad>> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  vector<complex<quad>> list;
  while (getline(iss, token, ','))
  {
    list.push_back(parse_value<complex<quad>>(trim(token)));
  }
  return list;
}

// Specialization for VectorXd
template <>
VectorXd parse_list(const string &str)
{
  // Use a temporary list to extract the values
  vector<double> temp_list = parse_list<vector<double>>(str);

  // Fill the VectorXd with the extracted values
  VectorXd list(temp_list.size());
  for (size_t i = 0; i < temp_list.size(); ++i)
  {
    list(i) = temp_list[i];
  }
  return list;
}

// Specialization for VectorXcd
template <>
VectorXcd parse_list(const string &str)
{
  // Use a temporary list to extract the values
  vector<complex<double>> temp_list = parse_list<vector<complex<double>>>(str);

  // Fill the VectorXd with the extracted values
  VectorXcd list(temp_list.size());
  for (size_t i = 0; i < temp_list.size(); ++i)
  {
    list(i) = temp_list[i];
  }
  return list;
}

// Specialization for VectorXld
template <>
VectorXld parse_list(const string &str)
{
  // Use a temporary list to extract the values
  vector<long double> temp_list = parse_list<vector<long double>>(str);

  // Fill the VectorXd with the extracted values
  VectorXld list(temp_list.size());
  for (size_t i = 0; i < temp_list.size(); ++i)
  {
    list(i) = temp_list[i];
  }
  return list;
}

// Specialization for VectorXcld
template <>
VectorXcld parse_list(const string &str)
{
  // Use a temporary list to extract the values
  vector<complex<long double>> temp_list = parse_list<vector<complex<long double>>>(str);

  // Fill the VectorXd with the extracted values
  VectorXcld list(temp_list.size());
  for (size_t i = 0; i < temp_list.size(); ++i)
  {
    list(i) = temp_list[i];
  }
  return list;
}

// Specialization for VectorXQ
template <>
VectorXQ parse_list(const string &str)
{
  // Use a temporary list to extract the values
  vector<quad> temp_list = parse_list<vector<quad>>(str);

  // Fill the VectorXd with the extracted values
  VectorXQ list(temp_list.size());
  for (size_t i = 0; i < temp_list.size(); ++i)
  {
    list(i) = temp_list[i];
  }
  return list;
}

// Specialization for VectorXcQ
template <>
VectorXcQ parse_list(const string &str)
{
  // Use a temporary list to extract the values
  vector<complex<quad>> temp_list = parse_list<vector<complex<quad>>>(str);

  // Fill the VectorXd with the extracted values
  VectorXcQ list(temp_list.size());
  for (size_t i = 0; i < temp_list.size(); ++i)
  {
    list(i) = temp_list[i];
  }
  return list;
}


// Helper functions for printing in the minimization routines
template <typename RealT>
void print_num_min(ostream &stream, const int &order, const int &no_cycles,
                   const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                   const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int steps,
                   const RealT *geod_eps)
{
  stream << "Minimize (" << steps << " step) results" << endl;
  stream << "Order: " << order << ", No. cycles: " << no_cycles << endl;
  stream << "Convergence criteria: - 1. step: " << scientific << setprecision(1) << "eps1 = " << eps1[0] << ", eps2 = " << eps1[1] << ", eps3 = " << eps1[2] << ", eps4 = " << eps1[3] << endl;

  if (steps == 2)
  {
    stream << string(22, ' ') << "- 2. step: " << scientific << setprecision(1) << "eps2 = " << eps2[0] << ", eps2 = " << eps2[1] << ", eps3 = " << eps2[2] << ", eps4 = " << eps2[3] << endl;
  }

  stream << "Weigths per order: " << scientific << setprecision(1) << "w2 = " << wi[0];

  if (order >= 4)
  {
    stream << ", w4 = " << wi[1];
  }
  if (order >= 6)
  {
    stream << ", w6 = " << wi[2];
  }
  if (order >= 8)
  {
    stream << ", w8 = " << wi[3];
  }

  stream << endl;
  stream << "Derivative step: " << scientific << setprecision(1) << step << endl;
  stream << "Number of maximal iterations: " << fixed << setprecision(1) << n_iter << endl;
  stream << "Ls = (" << fixed << setprecision(2) << Ls[0] << ", " << Ls[1] << "), lambda = " << lambda << endl;
  if (geod_eps)
  {
    stream << "geod_eps = " << scientific << setprecision(2) << *geod_eps << endl;
  }
}


// Helper functions for printing in the minimization routines (with origin term)
template <typename RealT>
void print_num_min_origin(ostream &stream, const int &order, const int &no_cycles,
                          const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                          const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const RealT ratio, const int steps)
{
  stream << "Minimize (" << steps << " step) results" << endl;
  stream << "Order: " << order << ", No. cycles: " << no_cycles << endl;
  stream << "Convergence criteria: - 1. step: " << scientific << setprecision(0) << "eps1 = " << eps1[0] << ", eps2 = " << eps1[1] << ", eps3 = " << eps1[2] << ", eps4 = " << eps1[3] << endl;

  if (steps == 2)
  {
    stream << string(22, ' ') << "- 2. step: " << scientific << setprecision(0) << "eps2 = " << eps2[0] << ", eps2 = " << eps2[1] << ", eps3 = " << eps2[2] << ", eps4 = " << eps2[3] << endl;
  }

  stream << "Weigths per order: " << scientific << setprecision(1) << "w2 = " << wi[0];

  if (order >= 4)
  {
    stream << ", w4 = " << wi[1];
  }
  if (order >= 6)
  {
    stream << ", w6 = " << wi[2];
  }
  if (order >= 8)
  {
    stream << ", w8 = " << wi[3];
  }

  stream << endl;
  stream << "Derivative step: " << scientific << setprecision(1) << step << endl;
  stream << "Number of maximal iterations: " << fixed << setprecision(1) << n_iter << endl;
  stream << "Ls = (" << fixed << setprecision(2) << Ls[0] << ", " << Ls[1] << "), lambda = " << lambda << endl;
  stream << "ratio = " << scientific << setprecision(5) << ratio << endl;
}


// Helper functions for printing in the minimization routines
template <typename RealT>
void print_sym_min(ostream &stream, const int &order, const int &no_cycles,
                   const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                   const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int steps)
{
  stream << "Minimize (" << steps << " steps) results" << endl;
  stream << "Order: " << order << ", No. cycles: " << no_cycles << endl;
  stream << "Convergence criteria: - 1. step: " << scientific << setprecision(0) << "eps1 = " << eps1[0] << ", eps2 = " << eps1[1] << ", eps3 = " << eps1[2] << ", eps4 = " << eps1[3] << endl;

  if (steps == 2)
  {
    stream << string(22, ' ') << "- 2. step: " << scientific << setprecision(0) << "eps2 = " << eps2[0] << ", eps2 = " << eps2[1] << ", eps3 = " << eps2[2] << ", eps4 = " << eps2[3] << endl;
  }

  stream << "Weigths per order: " << scientific << setprecision(1) << "w2 = " << wi[0];

  if (order >= 4)
  {
    stream << ", w4 = " << wi[1];
  }
  if (order >= 6)
  {
    stream << ", w6 = " << wi[2];
  }

  stream << endl;
  stream << "Number of maximal iterations: " << fixed << setprecision(1) << n_iter << endl;
  stream << "Ls = (" << fixed << setprecision(2) << Ls[0] << ", " << Ls[1] << "), lambda = " << lambda << endl;
}


// Specialization for NScheme
// Specialization for VectorXd
template <>
void print_num_min1(ostream &stream, NScheme<VectorXd> &scheme, const VectorXd &a_init, const VectorXd &b_init)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i];
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i];
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i];
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i];
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl;
}

// Specialization for VectorXcd
template <>
void print_num_min1(ostream &stream, NScheme<VectorXcd> &scheme, const VectorXcd &a_init, const VectorXcd &b_init)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i].real();
    if (a_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_init[i].imag());
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i].real();
    if (b_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_init[i].imag());
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i].real();
    if (scheme.a_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.a_eval[i].imag());
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i].real();
    if (scheme.b_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.b_eval[i].imag());
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl;
}

// Specialization for VectorXld
template <>
void print_num_min1(ostream &stream, NScheme<VectorXld> &scheme, const VectorXld &a_init, const VectorXld &b_init)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i];
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i];
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i];
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i];
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl;
}

// Specialization for VectorXcld
template <>
void print_num_min1(ostream &stream, NScheme<VectorXcld> &scheme, const VectorXcld &a_init, const VectorXcld &b_init)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i].real();
    if (a_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_init[i].imag());
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i].real();
    if (b_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_init[i].imag());
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i].real();
    if (scheme.a_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.a_eval[i].imag());
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i].real();
    if (scheme.b_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.b_eval[i].imag());
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl;
}

// Specialization for VectorXQ
template <>
void print_num_min1(ostream &stream, NScheme<VectorXQ> &scheme, const VectorXQ &a_init, const VectorXQ &b_init)
{
  stream << string(99, '-') << endl;
  stream << fixed << setprecision(36);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i];
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i];
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(99, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i];
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i];
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(99, '-') << endl;
  stream << scientific << setprecision(36);
  stream << "| Order |                  Error^2                   |                 Efficiency                 |" << endl;
  stream << string(99, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |                    0.0                     |                  Infinity                  |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(99, '-') << endl;
}

// Specialization for VectorXcQ
template <>
void print_num_min1(ostream &stream, NScheme<VectorXcQ> &scheme, const VectorXcQ &a_init, const VectorXcQ &b_init)
{
  stream << string(99, '-') << endl;
  stream << fixed << setprecision(36);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i].real();
    if (a_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_init[i].imag());
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i].real();
    if (b_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_init[i].imag());
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(99, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i].real();
    if (scheme.a_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.a_eval[i].imag());
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i].real();
    if (scheme.b_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.b_eval[i].imag());
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(99, '-') << endl;
  stream << scientific << setprecision(36);
  stream << "| Order |                  Error^2                   |                 Efficiency                 |" << endl;
  stream << string(99, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |                    0.0                     |                  Infinity                  |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(99, '-') << endl;
}


// Specialization for Scheme
// Specialization for VectorXd
template <>
void print_sym_min1(ostream &stream, Scheme<double> &scheme,
                    const VectorXd &a_init, const VectorXd &b_init, const VectorXd &a_vec, const VectorXd &b_vec)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i];
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i];
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i];
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i];
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl;
}

// Specialization for VectorXcd
template <>
void print_sym_min1(ostream &stream, Scheme<double> &scheme,
                    const VectorXcd &a_init, const VectorXcd &b_init, const VectorXcd &a_vec, const VectorXcd &b_vec)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i].real();
    if (a_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_init[i].imag());
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i].real();
    if (b_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_init[i].imag());
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i].real();
    if (a_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_vec[i].imag());
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i].real();
    if (b_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_vec[i].imag());
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl;
}

// Specialization for VectorXld
template <>
void print_sym_min1(ostream &stream, Scheme<long double> &scheme,
                    const VectorXld &a_init, const VectorXld &b_init, const VectorXld &a_vec, const VectorXld &b_vec)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i];
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i];
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i];
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i];
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl;
}

// Specialization for VectorXcld
template <>
void print_sym_min1(ostream &stream, Scheme<long double> &scheme,
                    const VectorXcld &a_init, const VectorXcld &b_init, const VectorXcld &a_vec, const VectorXcld &b_vec)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Initial vector a: ";
  for (int i = 0; i < a_init.size(); ++i)
  {
    stream << a_init[i].real();
    if (a_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_init[i].imag());
    if (i != a_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Initial vector b: ";
  for (int i = 0; i < b_init.size(); ++i)
  {
    stream << b_init[i].real();
    if (b_init[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_init[i].imag());
    if (i != b_init.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << "Final vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i].real();
    if (a_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_vec[i].imag());
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Final vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i].real();
    if (b_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_vec[i].imag());
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl;
}


// Helper functions for printing in the finder routine
template <typename Scalar, typename RealT>
void print_num_find(ostream &stream, const int &order, const int &no_cycles,
                    const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                    const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int N, const int steps,
                    Scalar mu, Scalar sigma, int sum_convs, const RealT tol,
                    const RealT *geod_eps)
{
  stream << "Minima finder results" << endl;
  stream << "Order: " << order << ", No. cycles: " << no_cycles << endl;
  stream << "Convergence criteria: - 1. step: " << scientific << setprecision(0) << "eps1 = " << eps1[0] << ", eps2 = " << eps1[1] << ", eps3 = " << eps1[2] << ", eps4 = " << eps1[3] << endl;

  if (steps == 2)
  {
    stream << string(22, ' ') << "- 2. step: " << scientific << setprecision(0) << "eps2 = " << eps2[0] << ", eps2 = " << eps2[1] << ", eps3 = " << eps2[2] << ", eps4 = " << eps2[3] << endl;
  }

  stream << "Weigths per order: " << scientific << setprecision(1) << "w2 = " << wi[0];

  if (order >= 4)
  {
    stream << ", w4 = " << wi[1];
  }
  if (order >= 6)
  {
    stream << ", w6 = " << wi[2];
  }
  if (order >= 8)
  {
    stream << ", w8 = " << wi[3];
  }

  stream << endl;
  stream << "Derivative step: " << scientific << setprecision(1) << step << endl;
  stream << "Number of maximal iterations: " << fixed << setprecision(1) << n_iter << endl;
  stream << "Ls = (" << fixed << setprecision(2) << Ls[0] << ", " << Ls[1] << "), lambda = " << lambda << endl;
  if (geod_eps)
  {
    stream << "geod_eps = " << scientific << setprecision(2) << *geod_eps << endl;
  }
  stream << "Number of initial conditions: " << fixed << N << endl;
  stream << "Number of steps taken: " << fixed << steps << endl;
  stream << "Comparison tolerance: " << scientific << setprecision(2) << tol << endl;

  stream << fixed << setprecision(4);
  if (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
  {
    stream << "Normal distribution: - mu = " << mu << endl;
    stream << string(21, ' ') << "- sigma = " << sigma << endl;
  }
  else if (is_same_v<Scalar, complex<double>> || is_same_v<Scalar, complex<long double>>)
  {
    stream << "Normal distribution: - mu = " << real(mu);
    if (imag(mu) > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(imag(mu)) << endl;

    stream << string(21, ' ') << "- sigma = " << real(sigma);
    if (imag(sigma) > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(imag(sigma)) << endl;
  }

  stream << string(38, '-') << endl;
  stream << "Ratio of unconverged samples: " << fixed << setprecision(2) << (static_cast<double>(N - sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;

  if constexpr (is_same_v<Scalar, double>)
  {
    stream << string(106, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |          a          |          b          |\n";
    stream << string(106, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, complex<double>>)
  {
    stream << string(150, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                     a                     |                     b                     |\n";
    stream << string(150, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, long double>)
  {
    stream << string(114, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |            a            |            b            |\n";
    stream << string(114, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, complex<long double>>)
  {
    stream << string(166, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                         a                         |                         b                         |\n";
    stream << string(166, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, quad>)
  {
    stream << string(146, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                    a                    |                    b                    |\n";
    stream << string(146, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, complex<quad>>)
  {
    stream << string(230, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                                         a                                         |                                         b                                         |\n";
    stream << string(230, '=') << endl;
  }
}


// Helper functions for printing in the finder routine (with origin term)
template <typename Scalar, typename RealT>
void print_num_find_origin(ostream &stream, const int &order, const int &no_cycles,
                           const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                           const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const RealT &ratio, const int N, const int steps,
                           Scalar mu, Scalar sigma, int sum_convs, const RealT tol)
{
  stream << "Minima finder results" << endl;
  stream << "Order: " << order << ", No. cycles: " << no_cycles << endl;
  stream << "Convergence criteria: - 1. step: " << scientific << setprecision(0) << "eps1 = " << eps1[0] << ", eps2 = " << eps1[1] << ", eps3 = " << eps1[2] << ", eps4 = " << eps1[3] << endl;

  if (steps == 2)
  {
    stream << string(22, ' ') << "- 2. step: " << scientific << setprecision(0) << "eps2 = " << eps2[0] << ", eps2 = " << eps2[1] << ", eps3 = " << eps2[2] << ", eps4 = " << eps2[3] << endl;
  }

  stream << "Weigths per order: " << scientific << setprecision(1) << "w2 = " << wi[0];

  if (order >= 4)
  {
    stream << ", w4 = " << wi[1];
  }
  if (order >= 6)
  {
    stream << ", w6 = " << wi[2];
  }
  if (order >= 8)
  {
    stream << ", w8 = " << wi[3];
  }

  stream << endl;
  stream << "Derivative step: " << scientific << setprecision(1) << step << endl;
  stream << "Number of maximal iterations: " << fixed << setprecision(1) << n_iter << endl;
  stream << "Ls = (" << fixed << setprecision(2) << Ls[0] << ", " << Ls[1] << "), lambda = " << lambda << endl;
  stream << "ratio a: " << scientific << setprecision(5) << ratio << endl;
  stream << "Number of initial conditions: " << fixed << N << endl;
  stream << "Number of steps taken: " << fixed << steps << endl;
  stream << "Comparison tolerance: " << scientific << setprecision(2) << tol << endl;

  stream << fixed << setprecision(4);
  if (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
  {
    stream << "Normal distribution: - mu = " << mu << endl;
    stream << string(21, ' ') << "- sigma = " << sigma << endl;
  }
  else if (is_same_v<Scalar, complex<double>> || is_same_v<Scalar, complex<long double>>)
  {
    stream << "Normal distribution: - mu = " << real(mu);
    if (imag(mu) > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(imag(mu)) << endl;

    stream << string(21, ' ') << "- sigma = " << real(sigma);
    if (imag(sigma) > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(imag(sigma)) << endl;
  }

  stream << string(38, '-') << endl;
  stream << "Ratio of unconverged samples: " << fixed << setprecision(2) << (static_cast<double>(N - sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;

  if constexpr (is_same_v<Scalar, double>)
  {
    stream << string(106, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |          a          |          b          |\n";
    stream << string(106, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, complex<double>>)
  {
    stream << string(150, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                     a                     |                     b                     |\n";
    stream << string(150, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, long double>)
  {
    stream << string(114, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |            a            |            b            |\n";
    stream << string(114, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, complex<long double>>)
  {
    stream << string(166, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                         a                         |                         b                         |\n";
    stream << string(166, '=') << endl;
  }
}


template <typename Scalar, typename RealT>
void print_sym_find(ostream &stream, const int &order, const int &no_cycles,
                    const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                    const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int N, const int steps,
                    Scalar mu, Scalar sigma, int sum_convs, const RealT tol)
{
  stream << "Minima finder results" << endl;
  stream << "Order: " << order << ", No. cycles: " << no_cycles << endl;
  stream << "Convergence criteria: - 1. step: " << scientific << setprecision(0) << "eps1 = " << eps1[0] << ", eps2 = " << eps1[1] << ", eps3 = " << eps1[2] << ", eps4 = " << eps1[3] << endl;

  if (steps == 2)
  {
    stream << string(22, ' ') << "- 2. step: " << scientific << setprecision(0) << "eps2 = " << eps2[0] << ", eps2 = " << eps2[1] << ", eps3 = " << eps2[2] << ", eps4 = " << eps2[3] << endl;
  }

  stream << "Weigths per order: " << scientific << setprecision(1) << "w2 = " << wi[0];

  if (order >= 4)
  {
    stream << ", w4 = " << wi[1];
  }
  if (order >= 6)
  {
    stream << ", w6 = " << wi[2];
  }

  stream << endl;
  stream << "Number of maximal iterations: " << fixed << setprecision(1) << n_iter << endl;
  stream << "Ls = (" << fixed << setprecision(2) << Ls[0] << ", " << Ls[1] << "), lambda = " << lambda << endl;
  stream << "Number of initial conditions: " << fixed << N << endl;
  stream << "Number of steps taken: " << fixed << steps << endl;
  stream << "Comparison tolerance: " << scientific << setprecision(2) << tol << endl;

  stream << fixed << setprecision(4);
  if (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
  {
    stream << "Normal distribution: - mu = " << mu << endl;
    stream << string(21, ' ') << "- sigma = " << sigma << endl;
  }
  else if (is_same_v<Scalar, complex<double>> || is_same_v<Scalar, complex<long double>>)
  {
    stream << "Normal distribution: - mu = " << real(mu);
    if (imag(mu) > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(imag(mu)) << endl;

    stream << string(21, ' ') << "- sigma = " << real(sigma);
    if (imag(sigma) > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(imag(sigma)) << endl;
  }

  stream << string(38, '-') << endl;
  stream << "Ratio of unconverged samples: " << fixed << setprecision(2) << (static_cast<double>(N - sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;

  if constexpr (is_same_v<Scalar, double>)
  {
    stream << string(106, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |          a          |          b          |\n";
    stream << string(106, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, complex<double>>)
  {
    stream << string(150, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                     a                     |                     b                     |\n";
    stream << string(150, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, long double>)
  {
    stream << string(114, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |            a            |            b            |\n";
    stream << string(114, '=') << endl;
  }
  else if constexpr (is_same_v<Scalar, complex<long double>>)
  {
    stream << string(166, '-') << endl;
    stream << "| i   ||  Ratio  |   Eff2   |   Eff4   |   Eff6   |   Eff8   |                         a                         |                         b                         |\n";
    stream << string(166, '=') << endl;
  }
}


template <typename RealT>
void print_find1(ostream &stream, const int &tracker, const double &ratio,
                 const RealT &eff2, const RealT &eff4, const RealT &eff6, const RealT &eff8)
{
  stream << "| " << tracker;

  if (tracker < 10)
  {
    stream << "   || ";
  }
  else if (tracker < 100)
  {
    stream << "  || ";
  }
  else
  {
    stream << " || ";
  }

  stream << fixed;
  if (ratio == 100.0)
  {
    stream << setprecision(1);
  }
  else if (ratio >= 10.0)
  {
    stream << setprecision(2);
  }
  else
  {
    stream << setprecision(3);
  }

  stream << ratio << " % | ";

  if (eff2 == INFINITY)
  {
    stream << "infinity | ";
  }
  else
  {
    if (eff2 > 1e3 || eff2 <= 1e-4)
    {
      stream << scientific << setprecision(2);
    }
    else if (eff2 > 1e2)
    {
      stream << fixed << setprecision(4);
    }
    else if (eff2 > 1e1)
    {
      stream << fixed << setprecision(5);
    }
    else
    {
      stream << fixed << setprecision(6);
    }
    stream << eff2 << " | ";
  }

  if (eff4 == INFINITY)
  {
    stream << "infinity | ";
  }
  else
  {
    if (eff4 > 1e3 || eff4 <= 1e-4)
    {
      stream << scientific << setprecision(2);
    }
    else if (eff4 > 1e2)
    {
      stream << fixed << setprecision(4);
    }
    else if (eff4 > 1e1)
    {
      stream << fixed << setprecision(5);
    }
    else
    {
      stream << fixed << setprecision(6);
    }
    stream << eff4 << " | ";
  }

  if (eff6 == INFINITY)
  {
    stream << "infinity | ";
  }
  else
  {
    if (eff6 > 1e3 || eff6 <= 1e-4)
    {
      stream << scientific << setprecision(2);
    }
    else if (eff6 > 1e2)
    {
      stream << fixed << setprecision(4);
    }
    else if (eff6 > 1e1)
    {
      stream << fixed << setprecision(5);
    }
    else
    {
      stream << fixed << setprecision(6);
    }
    stream << eff6 << " | ";
  }

  if (eff8 == INFINITY)
  {
    stream << "infinity | ";
  }
  else
  {
    if (eff8 > 1e3 || eff8 <= 1e-4)
    {
      stream << scientific << setprecision(2);
    }
    else if (eff8 > 1e2)
    {
      stream << fixed << setprecision(4);
    }
    else if (eff8 > 1e1)
    {
      stream << fixed << setprecision(5);
    }
    else
    {
      stream << fixed << setprecision(6);
    }
    stream << eff8 << " | ";
  }
}


// Specialization for VectorXd
template <>
void print_find2(ostream &stream, const VectorXd &a_vec, const VectorXd &b_vec)
{
  stream << fixed << setprecision(16);

  if (a_vec[0] >= 0)
  {
    stream << " ";
  }
  stream << a_vec[0] << " | ";

  if (b_vec[0] >= 0)
  {
    stream << " ";
  }
  stream << b_vec[0] << " |\n";

  if (a_vec.size() > b_vec.size())
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (a_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << a_vec[j] << " | ";
      if (j != a_vec.size() - 1)
      {
        if (b_vec[j] >= 0)
        {
          stream << " ";
        }
        stream << b_vec[j] << " |\n";
      }
      else
      {
        stream << string(20, ' ') << "|" << endl;
      }
    }
  }
  else if (b_vec.size() > a_vec.size())
  {
    for (int j{1}; j < b_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (j != b_vec.size() - 1)
      {
        if (a_vec[j] >= 0)
        {
          stream << " ";
        }
        stream << a_vec[j] << " | ";
      }
      else
      {
        stream << string(20, ' ') << "| ";
      }
      if (b_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << b_vec[j] << " |\n";
    }
  }
  else
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (a_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << a_vec[j] << " | ";
      if (b_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << b_vec[j] << " |\n";
    }
  }
  stream << string(106, '-') << endl;
}

// Specialization for VectorXcd
template <>
void print_find2(ostream &stream, const VectorXcd &a_vec, const VectorXcd &b_vec)
{
  stream << fixed << setprecision(16);

  if (real(a_vec[0]) >= 0)
  {
    stream << " ";
  }
  stream << real(a_vec[0]);

  if (imag(a_vec[0]) >= 0)
  {
    stream << " + i";
  }
  else
  {
    stream << " - i";
  }
  stream << abs(imag(a_vec[0])) << " | ";

  if (real(b_vec[0]) >= 0)
  {
    stream << " ";
  }
  stream << real(b_vec[0]);

  if (imag(b_vec[0]) >= 0)
  {
    stream << " + i";
  }
  else
  {
    stream << " - i";
  }
  stream << abs(imag(b_vec[0])) << " |\n";

  if (a_vec.size() > b_vec.size())
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (real(a_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(a_vec[j]);
      if (imag(a_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(a_vec[j])) << " | ";

      if (j != a_vec.size() - 1)
      {
        if (real(b_vec[j]) >= 0)
        {
          stream << " ";
        }
        stream << real(b_vec[j]);
        if (imag(b_vec[j]) >= 0)
        {
          stream << " + i";
        }
        else
        {
          stream << " - i";
        }
        stream << abs(imag(b_vec[j])) << " |\n";
      }
      else
      {
        stream << string(42, ' ') << "|" << endl;
      }
    }
  }
  else if (b_vec.size() > a_vec.size())
  {
    for (int j{1}; j < b_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (j != b_vec.size() - 1)
      {
        if (real(a_vec[j]) >= 0)
        {
          stream << " ";
        }
        stream << real(a_vec[j]);
        if (imag(a_vec[j]) >= 0)
        {
          stream << " + i";
        }
        else
        {
          stream << " - i";
        }
        stream << abs(imag(a_vec[j])) << " | ";
      }
      else
      {
        stream << string(42, ' ') << "| ";
      }
      if (real(b_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(b_vec[j]);
      if (imag(b_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(b_vec[j])) << " |\n";
    }
  }
  else
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (real(a_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(a_vec[j]);
      if (imag(a_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(a_vec[j])) << " | ";
      if (real(b_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(b_vec[j]);
      if (imag(b_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(b_vec[j])) << " |\n";
    }
  }
  stream << string(150, '-') << endl;
}

// Specialization for VectorXld
template <>
void print_find2(ostream &stream, const VectorXld &a_vec, const VectorXld &b_vec)
{
  stream << fixed << setprecision(20);

  if (a_vec[0] >= 0)
  {
    stream << " ";
  }
  stream << a_vec[0] << " | ";

  if (b_vec[0] >= 0)
  {
    stream << " ";
  }
  stream << b_vec[0] << " |\n";

  if (a_vec.size() > b_vec.size())
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (a_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << a_vec[j] << " | ";
      if (j != a_vec.size() - 1)
      {
        if (b_vec[j] >= 0)
        {
          stream << " ";
        }
        stream << b_vec[j] << " |\n";
      }
      else
      {
        stream << string(24, ' ') << "|" << endl;
      }
    }
  }
  else if (b_vec.size() > a_vec.size())
  {
    for (int j{1}; j < b_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (j != b_vec.size() - 1)
      {
        if (a_vec[j] >= 0)
        {
          stream << " ";
        }
        stream << a_vec[j] << " | ";
      }
      else
      {
        stream << string(24, ' ') << "| ";
      }
      if (b_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << b_vec[j] << " |\n";
    }
  }
  else
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (a_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << a_vec[j] << " | ";
      if (b_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << b_vec[j] << " |\n";
    }
  }
  stream << string(114, '-') << endl;
}

// Specialization for VectorXcld
template <>
void print_find2(ostream &stream, const VectorXcld &a_vec, const VectorXcld &b_vec)
{
  stream << fixed << setprecision(20);

  if (real(a_vec[0]) >= 0)
  {
    stream << " ";
  }
  stream << real(a_vec[0]);

  if (imag(a_vec[0]) >= 0)
  {
    stream << " + i";
  }
  else
  {
    stream << " - i";
  }
  stream << abs(imag(a_vec[0])) << " | ";

  if (real(b_vec[0]) >= 0)
  {
    stream << " ";
  }
  stream << real(b_vec[0]);

  if (imag(b_vec[0]) >= 0)
  {
    stream << " + i";
  }
  else
  {
    stream << " - i";
  }
  stream << abs(imag(b_vec[0])) << " |\n";

  if (a_vec.size() > b_vec.size())
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (real(a_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(a_vec[j]);
      if (imag(a_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(a_vec[j])) << " | ";

      if (j != a_vec.size() - 1)
      {
        if (real(b_vec[j]) >= 0)
        {
          stream << " ";
        }
        stream << real(b_vec[j]);
        if (imag(b_vec[j]) >= 0)
        {
          stream << " + i";
        }
        else
        {
          stream << " - i";
        }
        stream << abs(imag(b_vec[j])) << " |\n";
      }
      else
      {
        stream << string(50, ' ') << "|" << endl;
      }
    }
  }
  else if (b_vec.size() > a_vec.size())
  {
    for (int j{1}; j < b_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (j != b_vec.size() - 1)
      {
        if (real(a_vec[j]) >= 0)
        {
          stream << " ";
        }
        stream << real(a_vec[j]);
        if (imag(a_vec[j]) >= 0)
        {
          stream << " + i";
        }
        else
        {
          stream << " - i";
        }
        stream << abs(imag(a_vec[j])) << " | ";
      }
      else
      {
        stream << string(50, ' ') << "| ";
      }
      if (real(b_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(b_vec[j]);
      if (imag(b_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(b_vec[j])) << " |\n";
    }
  }
  else
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (real(a_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(a_vec[j]);
      if (imag(a_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(a_vec[j])) << " | ";
      if (real(b_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(b_vec[j]);
      if (imag(b_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(b_vec[j])) << " |\n";
    }
  }
  stream << string(166, '-') << endl;
}

// Specialization for VectorXQ
template <>
void print_find2(ostream &stream, const VectorXQ &a_vec, const VectorXQ &b_vec)
{
  stream << fixed << setprecision(36);

  if (a_vec[0] >= 0)
  {
    stream << " ";
  }
  stream << a_vec[0] << " | ";

  if (b_vec[0] >= 0)
  {
    stream << " ";
  }
  stream << b_vec[0] << " |\n";

  if (a_vec.size() > b_vec.size())
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (a_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << a_vec[j] << " | ";
      if (j != a_vec.size() - 1)
      {
        if (b_vec[j] >= 0)
        {
          stream << " ";
        }
        stream << b_vec[j] << " |\n";
      }
      else
      {
        stream << string(24, ' ') << "|" << endl;
      }
    }
  }
  else if (b_vec.size() > a_vec.size())
  {
    for (int j{1}; j < b_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (j != b_vec.size() - 1)
      {
        if (a_vec[j] >= 0)
        {
          stream << " ";
        }
        stream << a_vec[j] << " | ";
      }
      else
      {
        stream << string(24, ' ') << "| ";
      }
      if (b_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << b_vec[j] << " |\n";
    }
  }
  else
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (a_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << a_vec[j] << " | ";
      if (b_vec[j] >= 0)
      {
        stream << " ";
      }
      stream << b_vec[j] << " |\n";
    }
  }
  stream << string(146, '-') << endl;
}

// Specialization for VectorXcQ
template <>
void print_find2(ostream &stream, const VectorXcQ &a_vec, const VectorXcQ &b_vec)
{
  stream << fixed << setprecision(36);

  if (real(a_vec[0]) >= 0)
  {
    stream << " ";
  }
  stream << real(a_vec[0]);

  if (imag(a_vec[0]) >= 0)
  {
    stream << " + i";
  }
  else
  {
    stream << " - i";
  }
  stream << abs(imag(a_vec[0])) << " | ";

  if (real(b_vec[0]) >= 0)
  {
    stream << " ";
  }
  stream << real(b_vec[0]);

  if (imag(b_vec[0]) >= 0)
  {
    stream << " + i";
  }
  else
  {
    stream << " - i";
  }
  stream << abs(imag(b_vec[0])) << " |\n";

  if (a_vec.size() > b_vec.size())
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (real(a_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(a_vec[j]);
      if (imag(a_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(a_vec[j])) << " | ";

      if (j != a_vec.size() - 1)
      {
        if (real(b_vec[j]) >= 0)
        {
          stream << " ";
        }
        stream << real(b_vec[j]);
        if (imag(b_vec[j]) >= 0)
        {
          stream << " + i";
        }
        else
        {
          stream << " - i";
        }
        stream << abs(imag(b_vec[j])) << " |\n";
      }
      else
      {
        stream << string(50, ' ') << "|" << endl;
      }
    }
  }
  else if (b_vec.size() > a_vec.size())
  {
    for (int j{1}; j < b_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (j != b_vec.size() - 1)
      {
        if (real(a_vec[j]) >= 0)
        {
          stream << " ";
        }
        stream << real(a_vec[j]);
        if (imag(a_vec[j]) >= 0)
        {
          stream << " + i";
        }
        else
        {
          stream << " - i";
        }
        stream << abs(imag(a_vec[j])) << " | ";
      }
      else
      {
        stream << string(50, ' ') << "| ";
      }
      if (real(b_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(b_vec[j]);
      if (imag(b_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(b_vec[j])) << " |\n";
    }
  }
  else
  {
    for (int j{1}; j < a_vec.size(); ++j)
    {
      stream << "|" << string(5, ' ') << "||" << string(9, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "|" << string(10, ' ') << "| ";
      if (real(a_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(a_vec[j]);
      if (imag(a_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(a_vec[j])) << " | ";
      if (real(b_vec[j]) >= 0)
      {
        stream << " ";
      }
      stream << real(b_vec[j]);
      if (imag(b_vec[j]) >= 0)
      {
        stream << " + i";
      }
      else
      {
        stream << " - i";
      }
      stream << abs(imag(b_vec[j])) << " |\n";
    }
  }
  stream << string(230, '-') << endl;
}


// Specialization for VectorXd
template <>
void print_sym_scheme(ostream &stream, Scheme<double> &scheme, const VectorXd &a_vec, const VectorXd &b_vec)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i];
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i];
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl
         << endl;
}

// Specialization for VectorXcd
template <>
void print_sym_scheme(ostream &stream, Scheme<double> &scheme, const VectorXcd &a_vec, const VectorXcd &b_vec)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i].real();
    if (a_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_vec[i].imag());
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i].real();
    if (b_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_vec[i].imag());
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl
         << endl;
}

// Specialization for VectorXld
template <>
void print_sym_scheme(ostream &stream, Scheme<long double> &scheme, const VectorXld &a_vec, const VectorXld &b_vec)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i];
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i];
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl
         << endl;
}

// Specialization for VectorXcld
template <>
void print_sym_scheme(ostream &stream, Scheme<long double> &scheme, const VectorXcld &a_vec, const VectorXcld &b_vec)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Vector a: ";
  for (int i = 0; i < a_vec.size(); ++i)
  {
    stream << a_vec[i].real();
    if (a_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(a_vec[i].imag());
    if (i != a_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < b_vec.size(); ++i)
  {
    stream << b_vec[i].real();
    if (b_vec[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(b_vec[i].imag());
    if (i != b_vec.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2, a_vec, b_vec) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2, a_vec, b_vec) << " | " << scheme.eff(2, a_vec, b_vec) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4, a_vec, b_vec) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4, a_vec, b_vec) << " | " << scheme.eff(4, a_vec, b_vec) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6, a_vec, b_vec) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6, a_vec, b_vec) << " | " << scheme.eff(6, a_vec, b_vec) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl
         << endl;
}


// Specialization for VectorXd
template <>
void print_num_scheme(ostream &stream, NScheme<VectorXd> &scheme)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i];
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i];
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl
         << endl;
}

// Specialization for VectorXcd
template <>
void print_num_scheme(ostream &stream, NScheme<VectorXcd> &scheme)
{
  stream << string(59, '-') << endl;
  stream << fixed << setprecision(16);
  stream << "Vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i].real();
    if (scheme.a_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.a_eval[i].imag());
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i].real();
    if (scheme.b_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.b_eval[i].imag());
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(59, '-') << endl;
  stream << scientific << setprecision(16);
  stream << "| Order |        Error^2         |       Efficiency       |" << endl;
  stream << string(59, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |          0.0           |        Infinity        |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |          0.0           |        Infinity        |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(59, '-') << endl
         << endl;
}

// Specialization for VectorXld
template <>
void print_num_scheme(ostream &stream, NScheme<VectorXld> &scheme)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i];
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i];
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl
         << endl;
}

// Specialization for VectorXcld
template <>
void print_num_scheme(ostream &stream, NScheme<VectorXcld> &scheme)
{
  stream << string(67, '-') << endl;
  stream << fixed << setprecision(20);
  stream << "Vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i].real();
    if (scheme.a_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.a_eval[i].imag());
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i].real();
    if (scheme.b_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.b_eval[i].imag());
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(67, '-') << endl;
  stream << scientific << setprecision(20);
  stream << "| Order |          Error^2           |         Efficiency         |" << endl;
  stream << string(67, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |            0.0             |          Infinity          |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |            0.0             |          Infinity          |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(67, '-') << endl
         << endl;
}

// Specialization for VectorXQ
template <>
void print_num_scheme(ostream &stream, NScheme<VectorXQ> &scheme)
{
  stream << string(99, '-') << endl;
  stream << fixed << setprecision(36);
  stream << "Vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i];
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i];
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(99, '-') << endl;
  stream << scientific << setprecision(36);
  stream << "| Order |                  Error^2                   |                 Efficiency                 |" << endl;
  stream << string(99, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |                    0.0                     |                  Infinity                  |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(99, '-') << endl
         << endl;
}

// Specialization for VectorXcQ
template <>
void print_num_scheme(ostream &stream, NScheme<VectorXcQ> &scheme)
{
  stream << string(99, '-') << endl;
  stream << fixed << setprecision(36);
  stream << "Vector a: ";
  for (int i = 0; i < scheme.a_eval.size(); ++i)
  {
    stream << scheme.a_eval[i].real();
    if (scheme.a_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.a_eval[i].imag());
    if (i != scheme.a_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << "Vector b: ";
  for (int i = 0; i < scheme.b_eval.size(); ++i)
  {
    stream << scheme.b_eval[i].real();
    if (scheme.b_eval[i].imag() > 0)
    {
      stream << " + i";
    }
    else
    {
      stream << " - i";
    }
    stream << abs(scheme.b_eval[i].imag());
    if (i != scheme.b_eval.size() - 1)
    {
      stream << ", ";
    }
  }
  stream << endl;
  stream << string(99, '-') << endl;
  stream << scientific << setprecision(36);
  stream << "| Order |                  Error^2                   |                 Efficiency                 |" << endl;
  stream << string(99, '-') << endl;
  if (scheme.eff(2) == INFINITY)
  {
    stream << "|   2   |                    0.0                     |                  Infinity                  |" << endl;
  }
  else
  {
    stream << "|   2   | " << scheme.err2(2) << " | " << scheme.eff(2) << " |" << endl;
  }
  if (scheme.n >= 4)
  {
    if (scheme.eff(4) == INFINITY)
    {
      stream << "|   4   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   4   | " << scheme.err2(4) << " | " << scheme.eff(4) << " |" << endl;
    }
  }
  if (scheme.n >= 6)
  {
    if (scheme.eff(6) == INFINITY)
    {
      stream << "|   6   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   6   | " << scheme.err2(6) << " | " << scheme.eff(6) << " |" << endl;
    }
  }
  if (scheme.n >= 8)
  {
    if (scheme.eff(8) == INFINITY)
    {
      stream << "|   8   |                    0.0                     |                  Infinity                  |" << endl;
    }
    else
    {
      stream << "|   8   | " << scheme.err2(8) << " | " << scheme.eff(8) << " |" << endl;
    }
  }
  stream << string(99, '-') << endl
         << endl;
}

// Constructor for IO
IO::IO(const string read_file)
{
  rfile.open(read_file);
  if (!rfile)
  {
    runtime_error("Error opening input file.\n");
  }

  string line;
  while (getline(rfile, line) && line != "# Specific parameters")
  {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue; // Skip comments and empty lines

    size_t pos = line.find('=');
    if (pos == string::npos)
      continue; // Skip malformed lines

    string key = trim(line.substr(0, pos));
    string value = trim(line.substr(pos + 1));

    // Parse known keys
    if (key == "routine")
    {
      routine = stoi(value);
    }
    else if (key == "mode")
    {
      mode = stoi(value);
    }
    else if (key == "scalar_type")
    {
      scalar_type = value;
    }
    else if (key == "save_dir")
    {
      save_dir = value;
    }
    else if (key == "order")
    {
      order = stoi(value);
    }
    else if (key == "no_cycles")
    {
      no_cycles = stoi(value);
    }
  }
}

// Destructor for IO
IO::~IO()
{
}

// *** Routines ***

// ** Scheme computations **

// Symbolic Scheme computation
// Methods: - scratch: compute scheme from scratch
//          - load: load scheme from file
//          - load_iterate: load scheme from file and iterate for 2 cycles
// Optional methods (Execute if provided): - Compute errors and efficiencies (needs a_vec and b_vec)
//                                         - Save the scheme (needs a save folder)
template <typename Scalar>
void IO::sym_scheme()
{
  // Determine the Vec type from Scalar
  using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

  // Specific parameters
  string method; // Method to execute (scratch, load, load_iterate)

  // ai and bi evaluation vectors
  Vec a_eval;
  Vec b_eval;

  // Load directory
  string load_dir;

  // Iteration steps
  int steps;

  // Verbose option (0: nothing, 1: cycles, 2: values)
  char verbose;

  // Read the required parameters for the desired method
  string line;
  while (getline(rfile, line))
  {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;

    size_t pos = line.find('=');
    if (pos == string::npos)
      continue;

    string key = trim(line.substr(0, pos));
    string value = trim(line.substr(pos + 1));

    if (key == "method")
    {
      method = value;
    }
    // Evaluation vectors for ai and bi
    else if (key == "a_eval")
    {
      a_eval = parse_list<Vec>(value);
    }
    else if (key == "b_eval")
    {
      b_eval = parse_list<Vec>(value);
    }
    // Load directory
    else if (key == "load_dir")
    {
      load_dir = value;
    }
    // Iteration steps
    else if (key == "steps")
    {
      steps = stoi(value);
    }
    // Verbose option
    else if (key == "verbose")
    {
      verbose = stoi(value);
    }
  }

  // Close the read_file
  rfile.close();

  // Determine the scheme order based on the number of cycles (Only for printing)
  int scheme_order;
  if (no_cycles <= 2)
  {
    scheme_order = 2;
  }
  else if (no_cycles <= 6)
  {
    scheme_order = 4;
  }
  else if (no_cycles <= 14)
  {
    scheme_order = 6;
  }
  else if (no_cycles <= 30)
  {
    scheme_order = 8;
  }
  else
  {
    scheme_order = 10;
  }

  if (method == "scratch")
  {
    cout << "Running Symbolic Scheme up to order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(64, '=') << endl;
    cout << "Computing the scheme from scratch:" << endl;

    if (is_same_v<Scalar, double> || is_same_v<Scalar, complex<double>>)
    {
      Scheme<double> scheme(order, no_cycles, verbose);
      scheme.iterate();

      scheme.coefs.alpha.display_poly();
      scheme.coefs.beta.display_poly();

      if (no_cycles > 2)
      {
        for (size_t c{0}; c < scheme.coefs.gammas.size(); ++c)
        {
          scheme.coefs.gammas[c].display_poly();
        }
      }
      if (no_cycles > 6)
      {
        for (size_t c{0}; c < scheme.coefs.deltas.size(); ++c)
        {
          //scheme.coefs.deltas[c].display_poly();
        }
      }

      // Save the scheme if a save directory is provided
      if (!save_dir.empty())
      {
        cout << "Saving Scheme to " << save_dir << "saved_schemes/" << endl;
        scheme.save(save_dir + "saved_schemes/");
      }

      // Compute errors and efficiencies if evaluation vectors are provided
      if (a_eval.size() > 0 && b_eval.size() > 0)
      {
        if (!save_dir.empty())
        {
          string write_file = save_dir + "sym_scheme_q" + to_string(no_cycles) + ".out";
          ofstream wfile(write_file, ios::app);
          if (!wfile.is_open())
          {
            runtime_error("Error opening output file: " + write_file);
          }
          wfile << "Scheme evaluation results:" << endl;
          wfile << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          cout << "Scheme evaluation results:" << endl;
          cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          print_sym_scheme(wfile, scheme, a_eval, b_eval);
          print_sym_scheme(cout, scheme, a_eval, b_eval);
        }
        else
        {
          cout << "Scheme evaluation results:" << endl;
          cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;
          
          print_sym_scheme(cout, scheme, a_eval, b_eval);
        }
      }
    }
    else if (is_same_v<Scalar, long double> || is_same_v<Scalar, complex<long double>>)
    {
      Scheme<long double> scheme(order, no_cycles, verbose);
      scheme.iterate();

      // Save the scheme if a save directory is provided
      if (!save_dir.empty())
      {
        cout << "Saving Scheme to " << save_dir << "saved_schemes/" << endl;
        scheme.save(save_dir + "saved_schemes/");
      }

      // Compute errors and efficiencies if evaluation vectors are provided
      if (a_eval.size() > 0 && b_eval.size() > 0)
      {
        if (!save_dir.empty())
        {
          string write_file = save_dir + "sym_scheme_q" + to_string(no_cycles) + ".out";
          ofstream wfile(write_file, ios::app);
          if (!wfile.is_open())
          {
            runtime_error("Error opening output file: " + write_file);
          }
          wfile << "Scheme evaluation results:" << endl;
          wfile << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          cout << "Scheme evaluation results:" << endl;
          cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          print_sym_scheme(wfile, scheme, a_eval, b_eval);
          print_sym_scheme(cout, scheme, a_eval, b_eval);
        }
        else
        {
          cout << "Scheme evaluation results:" << endl;
          cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          print_sym_scheme(cout, scheme, a_eval, b_eval);
        }
      }
    }

    cout << string(64, '=') << endl;
  }
  else if (method == "load")
  {
    cout << "Running Symbolic Scheme up to order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(64, '=') << endl;
    cout << "Loading the scheme from " << load_dir << endl;

    if (is_same_v<Scalar, double> || is_same_v<Scalar, complex<double>>)
    {
      Scheme<double> scheme(load_dir, order, no_cycles, verbose);

      // Compute errors and efficiencies if evaluation vectors are provided
      if (a_eval.size() > 0 && b_eval.size() > 0)
      {
        scheme.display_err2(a_eval, b_eval);
        cout << string(41, '-') << endl;
        scheme.display_eff(a_eval, b_eval);
      }
    }
    else if (is_same_v<Scalar, long double> || is_same_v<Scalar, complex<long double>>)
    {
      Scheme<long double> scheme(load_dir, order, no_cycles, verbose);

      // Compute errors and efficiencies if evaluation vectors are provided
      if (a_eval.size() > 0 && b_eval.size() > 0)
      {
        if (!save_dir.empty())
        {
          string write_file = save_dir + "sym_scheme_q" + to_string(no_cycles) + ".out";
          ofstream wfile(write_file, ios::app);
          if (!wfile.is_open())
          {
            runtime_error("Error opening output file: " + write_file);
          }
          wfile << "Scheme evaluation results:" << endl;
          wfile << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          cout << "Scheme evaluation results:" << endl;
          cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          print_sym_scheme(wfile, scheme, a_eval, b_eval);
          print_sym_scheme(cout, scheme, a_eval, b_eval);
        }
        else
        {
          cout << "Scheme evaluation results:" << endl;
          cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

          print_sym_scheme(cout, scheme, a_eval, b_eval);
        }
      }
    }

    cout << string(64, '=') << endl;
  }
  else if (method == "load_iterate")
  {
    cout << "Running Symbolic Scheme up to order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(64, '=') << endl;

    cout << "Loading the scheme from " << load_dir << endl;
    if (is_same_v<Scalar, double> || is_same_v<Scalar, complex<double>>)
    {
      Scheme<double> scheme(load_dir, order, no_cycles, verbose);

      cout << "Iterating " << 2 * steps << " cycles:" << endl;
      scheme.q += 2 * steps;
      no_cycles += 2 * steps;
      scheme.transform(scheme.q);
      scheme.iterate(steps);

      // Save the scheme if a save directory is provided
      if (!save_dir.empty())
      {
        cout << "Saving Scheme to " << save_dir << "saved_schemes/" << endl;
        scheme.save(save_dir + "saved_schemes/");
      }
    }
    else if (is_same_v<Scalar, long double> || is_same_v<Scalar, complex<long double>>)
    {
      Scheme<long double> scheme(load_dir, order, no_cycles, verbose);

      cout << "Iterating " << 2 * steps << " cycles:" << endl;
      scheme.q += 2 * steps;
      no_cycles += 2 * steps;
      scheme.transform(scheme.q);
      scheme.iterate(steps);

      // Save the scheme if a save directory is provided
      if (!save_dir.empty())
      {
        cout << "Saving Scheme to " << save_dir << "saved_schemes/" << endl;
        scheme.save(save_dir + "saved_schemes/");
      }
    }

    cout << string(64, '=') << endl;
  }
  else
  {
    runtime_error("Method not supported!");
  }
}


// Numeric Scheme
// Computes the errors and efficiencies of the scheme for the provided a_vec and b_vec
template <typename Scalar>
void IO::num_scheme()
{
  // Determine the Vec type from Scalar
  using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

  // Specific paramters
  // Evaluation vectors for ai and bi
  Vec a_eval;
  Vec b_eval;
  // Verbose option (0-nothing, 1-cycles, 2-values)
  char verbose;

  // Read the required parameters for the desired method
  string line;
  while (getline(rfile, line))
  {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;

    size_t pos = line.find('=');
    if (pos == string::npos)
      continue;

    string key = trim(line.substr(0, pos));
    string value = trim(line.substr(pos + 1));

    // Evaluation vectors for ai and bi
    if (key == "a_eval")
    {
      a_eval = parse_list<Vec>(value);
    }
    else if (key == "b_eval")
    {
      b_eval = parse_list<Vec>(value);
    }
    // Verbose option
    else if (key == "verbose")
    {
      verbose = stoi(value);
    }
  }

  // Close the read_file
  rfile.close();

  // Determine the scheme order based on the number of cycles (Only for printing)
  int scheme_order;
  if (no_cycles <= 2)
  {
    scheme_order = 2;
  }
  else if (no_cycles <= 6)
  {
    scheme_order = 4;
  }
  else if (no_cycles <= 14)
  {
    scheme_order = 6;
  }
  else if (no_cycles <= 30)
  {
    scheme_order = 8;
  }
  else
  {
    scheme_order = 10;
  }

  cout << "Running Numeric Scheme for order " << order << " and number of cycles q = " << no_cycles << endl;
  cout << string(64, '=') << endl;

  NScheme<Vec> scheme(order, no_cycles, a_eval, b_eval, verbose);
  scheme.iterate();

  if (!save_dir.empty())
  {
    string write_file = save_dir + "num_scheme_q" + to_string(no_cycles) + ".out";
    ofstream wfile(write_file, ios::app);
    if (!wfile.is_open())
    {
      runtime_error("Error opening output file: " + write_file);
    }
    wfile << "Scheme evaluation results:" << endl;
    wfile << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

    cout << "Scheme evaluation results:" << endl;
    cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

    print_num_scheme(wfile, scheme);
    print_num_scheme(cout, scheme);
  }
  else
  {
    cout << "Scheme evaluation results:" << endl;
    cout << "Order: " << scheme_order << ", No. cycles: " << no_cycles << endl;

    print_num_scheme(cout, scheme);
  }

  cout << string(64, '=') << endl;
}

// ** Minimizations **


// Numeric Minimization
// Methods: - minimize: minimize the computed scheme manifold once using initial vector a_vec, b_vec
//          - min_twostep: minimize the computed scheme manifold twice (constraint minimization in the 2. step)
//          - find: find as many minima for the computed scheme manifold
template <typename Scalar>
void IO::num_minim()
{
  // Determine the Vec type from Scalar
  using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;

  // Specific parameters
  string method; // Method to use for minimization (minimize, min_twostep, find)

  // Initial ai and bi vectors
  Vec a_init;
  Vec b_init;

  // Minimization hyperparameters
  array<RealT, 4> eps1;            // Convergence criteria for the 1. step (eps1, eps2, eps3, eps4)
  array<RealT, 4> eps2;            // Convergence criteria for the 2. step (eps1, eps2, eps3, eps4)
  vector<RealT> wi;                // Weigths per order (w2, w4, w6, ...)
  RealT step{1e-5};                // Numeric derivative step
  int n_iter{500};                 // Number of iterations before quitting if convergence isn't achieved
  array<RealT, 2> Ls{{9.0, 11.0}}; // L_up, L_down
  RealT lambda{0.25};              // Damping parameter

  // Geodesic acceleration parameters
  RealT geod_eps{0.1}; // Convergence parameter for the geodesic acceleration (2|a|/|v| <= geod_eps)

  // Minima finder parameters
  int N;                  // Number of initial conditions
  int steps;              // Number of steps (1: minimize, 2: min_twostep)
  RealT threshold{1e-18}; // Threshold to pass the minimal error
  // Normal distribution parameters
  Scalar mu;
  Scalar sigma;

  // Sort the found minima accoring to Chi2
  bool bsort{false};

  // Verbose option (false: nothing, true: Change of ai, bi during iterations)
  bool verbose;

  // Minima comparison parameters
  RealT tol{1e-15};

  // Hessian freeze option
  bool freeze{false};

  // Ratio for the origin methods
  RealT ratio{0.0};

  // Read the required parameters for the desired method
  string line;
  while (getline(rfile, line))
  {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;

    size_t pos = line.find('=');
    if (pos == string::npos)
      continue;

    string key = trim(line.substr(0, pos));
    string value = trim(line.substr(pos + 1));

    if (key == "method")
    {
      method = value;
    }
    // Initial ai and bi vectors
    else if (key == "a_init")
    {
      a_init = parse_list<Vec>(value);
    }
    else if (key == "b_init")
    {
      b_init = parse_list<Vec>(value);
    }
    // Minimization hyperparameters
    else if (key == "eps1")
    {
      eps1 = parse_list<array<RealT, 4>>(value);
    }
    else if (key == "eps2")
    {
      eps2 = parse_list<array<RealT, 4>>(value);
    }
    else if (key == "wi")
    {
      wi = parse_list<vector<RealT>>(value);
    }
    else if (key == "step")
    {
      step = parse_value<RealT>(value);
    }
    else if (key == "n_iter")
    {
      n_iter = stoi(value);
    }
    else if (key == "Ls")
    {
      Ls = parse_list<array<RealT, 2>>(value);
    }
    else if (key == "lambda")
    {
      lambda = parse_value<RealT>(value);
    }
    else if (key == "geod_eps")
    {
      geod_eps = parse_value<RealT>(value);
    }
    else if (key == "verbose")
    {
      verbose = stoi(value);
    }
    // Minima finder parameters
    else if (key == "N")
    {
      N = stoi(value);
    }
    else if (key == "steps")
    {
      steps = stoi(value);
    }
    else if (key == "threshold")
    {
      threshold = parse_value<RealT>(value);
    }
    else if (key == "mu")
    {
      mu = parse_value<Scalar>(value);
    }
    else if (key == "sigma")
    {
      sigma = parse_value<Scalar>(value);
    }
    else if (key == "bsort")
    {
      bsort = stoi(value);
    }
    // Minima comparison parameters
    else if (key == "tol")
    {
      tol = parse_value<RealT>(value);
    }
    // Hessian freeze option
    else if (key == "freeze")
    {
      freeze = stoi(value);
    }
    else if (key == "ratio")
    {
      ratio = parse_value<RealT>(value);
    }
  }

  // Close the read_file
  rfile.close();

  // Fill the weights vector according to order
  std::vector<RealT> W_vec;
  W_vec.push_back(wi[0]);
  W_vec.push_back(wi[0]);
  if (order == 4)
  {
    for (size_t i{0}; i < 6; ++i)
    {
      W_vec.push_back(wi[1]);
    }
  }
  if (order == 6)
  {
    for (size_t i{0}; i < 6; ++i)
    {
      W_vec.push_back(wi[1]);
    }
    for (size_t i{0}; i < 18; ++i)
    {
      W_vec.push_back(wi[2]);
    }
  }
  if (order == 8)
  {
    for (size_t i{0}; i < 6; ++i)
    {
      W_vec.push_back(wi[1]);
    }
    for (size_t i{0}; i < 18; ++i)
    {
      W_vec.push_back(wi[2]);
    }
    for (size_t i{0}; i < 56; ++i)
    {
      W_vec.push_back(wi[3]);
    }
  }

  if (method == "minimize")
  {
    cout << "Running Numeric Minimization (minimize) for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(81, '=') << endl;

    NMinimization<Vec> minim(a_init, b_init, step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    MinResult<Vec> mini{minim.minimize(lambda, nullptr, verbose)};

    if (verbose)
    {
      mini.display();
      cout << string(81, '-') << endl;
    }

    NScheme<Vec> scheme(8, no_cycles, mini.a_vec, mini.b_vec);
    scheme.iterate();

    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_minimize_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_min(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 1);
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 1);

      print_num_min1(wfile, scheme, a_init, b_init);
      print_num_min1(cout, scheme, a_init, b_init);

      wfile << "\n\n";
    }
    else
    {
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 1);
      print_num_min1(cout, scheme, a_init, b_init);
    }

    cout << string(81, '=') << endl;
  }
  else if (method == "minimize_geodesic")
  {
    cout << "Running Numeric Minimization (minimize_geodesic) for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(81, '=') << endl;

    NMinimization<Vec> minim(a_init, b_init, step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    MinResult<Vec> mini{minim.minimize_geodesic(lambda, nullptr, verbose, geod_eps)};

    if (verbose)
    {
      mini.display();
      cout << string(81, '-') << endl;
    }

    NScheme<Vec> scheme(8, no_cycles, mini.a_vec, mini.b_vec);
    scheme.iterate();

    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_minimize_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_min(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 1, &geod_eps);
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 1, &geod_eps);

      print_num_min1(wfile, scheme, a_init, b_init);
      print_num_min1(cout, scheme, a_init, b_init);

      wfile << "\n\n";
    }
    else
    {
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 1, &geod_eps);
      print_num_min1(cout, scheme, a_init, b_init);
    }

    cout << string(81, '=') << endl;
  }
  else if (method == "minimize_origin")
  {
    cout << "Running Numeric Minimization (minimize_origin) for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(81, '=') << endl;

    // Add the additional weight
    if (order == 4)
    {
      W_vec.push_back(wi[1]);
    }
    if (order == 6)
    {
      W_vec.push_back(wi[2]);
    }
    if (order == 8)
    {
      W_vec.push_back(wi[3]);
    }

    NMinimization<Vec> minim(a_init, b_init, step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    MinResult<Vec> mini{minim.minimize_origin(lambda, ratio, nullptr, verbose)};

    if (verbose)
    {
      mini.display();
      cout << string(81, '-') << endl;
    }

    NScheme<Vec> scheme(8, no_cycles, mini.a_vec, mini.b_vec);
    scheme.iterate();

    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_minimize_origin_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_min_origin(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, ratio, 1);
      print_num_min_origin(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, ratio, 1);

      print_num_min1(wfile, scheme, a_init, b_init);
      print_num_min1(cout, scheme, a_init, b_init);

      wfile << "\n\n";
    }
    else
    {
      print_num_min_origin(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, ratio, 1);
      print_num_min1(cout, scheme, a_init, b_init);
    }

    cout << string(81, '=') << endl;
  }
  else if (method == "min_twostep")
  {
    cout << "Running Numeric Minimization (min_twostep) for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(81, '=') << endl;

    NMinimization<Vec> minim(a_init, b_init, step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    MinResult<Vec> mini{minim.min_twostep(lambda, &eps2, verbose, freeze)};

    if (verbose)
    {
      mini.display();
      cout << string(81, '-') << endl;
    }

    NScheme<Vec> scheme(8, no_cycles, mini.a_vec, mini.b_vec);
    scheme.iterate();

    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_min_twostep_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_min(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 2);
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 2);

      print_num_min1(wfile, scheme, a_init, b_init);
      print_num_min1(cout, scheme, a_init, b_init);

      wfile << "\n\n";
    }
    else
    {
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 2);

      print_num_min1(cout, scheme, a_init, b_init);
    }

    cout << string(81, '=') << endl;
  }
  else if (method == "min_twostep_geodesic")
  {
    cout << "Running Numeric Minimization (min_twostep_geodesic) for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(81, '=') << endl;

    NMinimization<Vec> minim(a_init, b_init, step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    MinResult<Vec> mini{minim.min_twostep_geodesic(lambda, &eps2, verbose, freeze, geod_eps)};

    if (verbose)
    {
      mini.display();
      cout << string(81, '-') << endl;
    }

    NScheme<Vec> scheme(8, no_cycles, mini.a_vec, mini.b_vec);
    scheme.iterate();

    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_min_twostep_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_min(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 2, &geod_eps);
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 2, &geod_eps);

      print_num_min1(wfile, scheme, a_init, b_init);
      print_num_min1(cout, scheme, a_init, b_init);

      wfile << "\n\n";
    }
    else
    {
      print_num_min(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, 2, &geod_eps);

      print_num_min1(cout, scheme, a_init, b_init);
    }

    cout << string(81, '=') << endl;
  }
  else if (method == "min_twostep_origin")
  {
    cout << "Running Numeric Minimization (min_twostep_origin) for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(81, '=') << endl;

    // Add the additional weight
    if (order == 4)
    {
      W_vec.push_back(wi[1]);
    }
    if (order == 6)
    {
      W_vec.push_back(wi[2]);
    }
    if (order == 8)
    {
      W_vec.push_back(wi[3]);
    }

    NMinimization<Vec> minim(a_init, b_init, step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    MinResult<Vec> mini{minim.min_twostep_origin(lambda, ratio, &eps2, verbose, freeze)};

    if (verbose)
    {
      mini.display();
      cout << string(81, '-') << endl;
    }

    NScheme<Vec> scheme(8, no_cycles, mini.a_vec, mini.b_vec);
    scheme.iterate();

    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_min_twostep_origin_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_min_origin(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, ratio, 2);
      print_num_min_origin(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, ratio, 2);

      print_num_min1(wfile, scheme, a_init, b_init);
      print_num_min1(cout, scheme, a_init, b_init);

      wfile << "\n\n";
    }
    else
    {
      print_num_min_origin(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, ratio, 2);

      print_num_min1(cout, scheme, a_init, b_init);
    }

    cout << string(81, '=') << endl;
  }
  else if (method == "find")
  {
    cout << "Running Numeric Minimization (find, " << steps << ") for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(88, '=') << endl;

    NMinimization<Vec> minim(step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    pair<array<vector<Vec>, 4>, vector<int>> minima;

    if (steps == 1)
    {
      minima = minim.find(N, lambda, steps, mu, sigma, nullptr, true, verbose, tol, freeze);
    }
    else if (steps == 2)
    {
      minima = minim.find(N, lambda, steps, mu, sigma, &eps2, true, verbose, tol, freeze);
    }

    array<vector<Vec>, 4> ab_vecs{minima.first};

    int sum_convs{accumulate(minima.second.begin(), minima.second.end(), 0)};
    // Write to file if specified otherwise just output to terminal
    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_minim_find_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_find(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, N, steps, mu, sigma, sum_convs, tol);
      print_num_find(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, N, steps, mu, sigma, sum_convs, tol);

      double ratioN;
      RealT eff2;
      RealT eff4;
      RealT eff6;
      RealT eff8;
      int tracker{0};
      for (size_t i{0}; i < ab_vecs[0].size(); ++i)
      {
        NScheme<Vec> scheme(8, no_cycles, ab_vecs[0][i], ab_vecs[2][i]);
        scheme.iterate();

        ratioN = (static_cast<double>(minima.second[i]) / static_cast<double>(N)) * 100.0;
        eff2 = 1 / (pow(no_cycles, 2) * sqrt(scheme.err2(2)));
        eff4 = 1 / (pow(no_cycles, 4) * sqrt(scheme.err2(4)));
        eff6 = 1 / (pow(no_cycles, 6) * sqrt(scheme.err2(6)));
        eff8 = 1 / (pow(no_cycles, 8) * sqrt(scheme.err2(8)));

        bool pass{true};
        if (order == 4 && scheme.err2(2) > threshold)
        {
          pass = false;
        }
        else if (order == 6 && scheme.err2(2) > threshold && scheme.err2(4) > threshold)
        {
          pass = false;
        }
        else if (order == 8 && scheme.err2(2) > threshold && scheme.err2(4) > threshold && scheme.err2(6) > threshold)
        {
          pass = false;
        }

        if (pass)
        {
          tracker += 1;
          sum_convs -= minima.second[i];

          print_find1(cout, tracker, ratioN, eff2, eff4, eff6, eff8);
          print_find1(wfile, tracker, ratioN, eff2, eff4, eff6, eff8);

          print_find2(cout, ab_vecs[0][i], ab_vecs[2][i]);
          print_find2(wfile, ab_vecs[0][i], ab_vecs[2][i]);
        }
      }

      cout << "Ratio of discarded samples: " << fixed << setprecision(2) << (static_cast<double>(sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;
      cout << string(38, '-') << endl;

      wfile << "Ratio of discarded samples: " << fixed << setprecision(2) << (static_cast<double>(sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;
      wfile << string(38, '-') << endl;

      wfile << "\n\n";

      wfile.close();
    }
    else
    {

      print_num_find(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, N, steps, mu, sigma, sum_convs, tol);

      double ratioN;
      RealT eff2;
      RealT eff4;
      RealT eff6;
      RealT eff8;
      int tracker{0};
      for (size_t i{0}; i < ab_vecs[0].size(); ++i)
      {
        NScheme<Vec> scheme(8, no_cycles, ab_vecs[0][i], ab_vecs[2][i]);
        scheme.iterate();

        ratioN = (static_cast<double>(minima.second[i]) / static_cast<double>(N)) * 100.0;
        eff2 = 1 / (pow(no_cycles, 2) * sqrt(scheme.err2(2)));
        eff4 = 1 / (pow(no_cycles, 4) * sqrt(scheme.err2(4)));
        eff6 = 1 / (pow(no_cycles, 6) * sqrt(scheme.err2(6)));
        eff8 = 1 / (pow(no_cycles, 8) * sqrt(scheme.err2(8)));

        bool pass{true};
        if (order == 4 && scheme.err2(2) > threshold)
        {
          pass = false;
        }
        else if (order == 6 && scheme.err2(2) > threshold && scheme.err2(4) > threshold)
        {
          pass = false;
        }
        else if (order == 8 && scheme.err2(2) > threshold && scheme.err2(4) > threshold && scheme.err2(6) > threshold)
        {
          pass = false;
        }

        if (pass)
        {
          tracker += 1;
          sum_convs -= minima.second[i];

          print_find1(cout, tracker, ratioN, eff2, eff4, eff6, eff8);
          print_find2(cout, ab_vecs[0][i], ab_vecs[2][i]);
        }
      }

      cout << "Ratio of discarded samples: " << fixed << setprecision(2) << (static_cast<double>(sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;
      cout << string(38, '-') << endl;
    }
    cout << string(88, '=') << endl;
  }
  else if (method == "find_geodesic")
  {
    cout << "Running Numeric Minimization (find_geodesic, " << steps << ") for order " << order << " and number of cycles q = " << no_cycles << endl;
    cout << string(88, '=') << endl;

    NMinimization<Vec> minim(step, order, no_cycles, W_vec, n_iter, eps1, Ls);
    pair<array<vector<Vec>, 4>, vector<int>> minima;

    if (steps == 1)
    {
      minima = minim.find_geodesic(N, lambda, steps, mu, sigma, nullptr, true, verbose, tol, freeze, geod_eps);
    }
    else if (steps == 2)
    {
      minima = minim.find_geodesic(N, lambda, steps, mu, sigma, &eps2, true, verbose, tol, freeze, geod_eps);
    }

    array<vector<Vec>, 4> ab_vecs{minima.first};

    int sum_convs{accumulate(minima.second.begin(), minima.second.end(), 0)};
    // Write to file if specified otherwise just output to terminal
    if (!save_dir.empty())
    {
      string write_file = save_dir + "num_minim_find_q" + to_string(no_cycles) + ".out";
      ofstream wfile(write_file, ios::app);
      if (!wfile.is_open())
      {
        runtime_error("Error opening output file: " + write_file);
      }

      print_num_find(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, N, steps, mu, sigma, sum_convs, tol, &geod_eps);
      print_num_find(wfile, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, N, steps, mu, sigma, sum_convs, tol, &geod_eps);

      double ratioN;
      RealT eff2;
      RealT eff4;
      RealT eff6;
      RealT eff8;
      int tracker{0};
      for (size_t i{0}; i < ab_vecs[0].size(); ++i)
      {
        NScheme<Vec> scheme(8, no_cycles, ab_vecs[0][i], ab_vecs[2][i]);
        scheme.iterate();

        ratioN = (static_cast<double>(minima.second[i]) / static_cast<double>(N)) * 100.0;
        eff2 = 1 / (pow(no_cycles, 2) * sqrt(scheme.err2(2)));
        eff4 = 1 / (pow(no_cycles, 4) * sqrt(scheme.err2(4)));
        eff6 = 1 / (pow(no_cycles, 6) * sqrt(scheme.err2(6)));
        eff8 = 1 / (pow(no_cycles, 8) * sqrt(scheme.err2(8)));

        bool pass{true};
        if (order == 4 && scheme.err2(2) > threshold)
        {
          pass = false;
        }
        else if (order == 6 && scheme.err2(2) > threshold && scheme.err2(4) > threshold)
        {
          pass = false;
        }
        else if (order == 8 && scheme.err2(2) > threshold && scheme.err2(4) > threshold && scheme.err2(6) > threshold)
        {
          pass = false;
        }

        if (pass)
        {
          tracker += 1;
          sum_convs -= minima.second[i];

          print_find1(cout, tracker, ratioN, eff2, eff4, eff6, eff8);
          print_find1(wfile, tracker, ratioN, eff2, eff4, eff6, eff8);

          print_find2(cout, ab_vecs[0][i], ab_vecs[2][i]);
          print_find2(wfile, ab_vecs[0][i], ab_vecs[2][i]);
        }
      }

      cout << "Ratio of discarded samples: " << fixed << setprecision(2) << (static_cast<double>(sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;
      cout << string(38, '-') << endl;

      wfile << "Ratio of discarded samples: " << fixed << setprecision(2) << (static_cast<double>(sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;
      wfile << string(38, '-') << endl;

      wfile << "\n\n";

      wfile.close();
    }
    else
    {

      print_num_find(cout, order, no_cycles, eps1, eps2, wi, step, n_iter, Ls, lambda, N, steps, mu, sigma, sum_convs, tol, &geod_eps);

      double ratioN;
      RealT eff2;
      RealT eff4;
      RealT eff6;
      RealT eff8;
      int tracker{0};
      for (size_t i{0}; i < ab_vecs[0].size(); ++i)
      {
        NScheme<Vec> scheme(8, no_cycles, ab_vecs[0][i], ab_vecs[2][i]);
        scheme.iterate();

        ratioN = (static_cast<double>(minima.second[i]) / static_cast<double>(N)) * 100.0;
        eff2 = 1 / (pow(no_cycles, 2) * sqrt(scheme.err2(2)));
        eff4 = 1 / (pow(no_cycles, 4) * sqrt(scheme.err2(4)));
        eff6 = 1 / (pow(no_cycles, 6) * sqrt(scheme.err2(6)));
        eff8 = 1 / (pow(no_cycles, 8) * sqrt(scheme.err2(8)));

        bool pass{true};
        if (order == 4 && scheme.err2(2) > threshold)
        {
          pass = false;
        }
        else if (order == 6 && scheme.err2(2) > threshold && scheme.err2(4) > threshold)
        {
          pass = false;
        }
        else if (order == 8 && scheme.err2(2) > threshold && scheme.err2(4) > threshold && scheme.err2(6) > threshold)
        {
          pass = false;
        }

        if (pass)
        {
          tracker += 1;
          sum_convs -= minima.second[i];

          print_find1(cout, tracker, ratioN, eff2, eff4, eff6, eff8);
          print_find2(cout, ab_vecs[0][i], ab_vecs[2][i]);
        }
      }

      cout << "Ratio of discarded samples: " << fixed << setprecision(2) << (static_cast<double>(sum_convs) / static_cast<double>(N)) * 100.0 << "%" << endl;
      cout << string(38, '-') << endl;
    }
    cout << string(88, '=') << endl;
  }
  else
  {
    runtime_error("Method not supported!");
  }
}


// Explicit routine instantiations for double
template void IO::sym_scheme<double>();
template void IO::num_scheme<double>();
template void IO::num_minim<double>();

// Explicit routine instantiations for complex<double>
template void IO::sym_scheme<complex<double>>();
template void IO::num_scheme<complex<double>>();
template void IO::num_minim<complex<double>>();


// Explicit routine instantiations for long double
template void IO::sym_scheme<long double>();
template void IO::num_scheme<long double>();
template void IO::num_minim<long double>();

// Explicit routine instantiations for complex<long double>
template void IO::sym_scheme<complex<long double>>();
template void IO::num_scheme<complex<long double>>();
template void IO::num_minim<complex<long double>>();

// Explicit routine instantiations for quad
//template void IO::sym_scheme<quad>();
template void IO::num_scheme<quad>();
template void IO::num_minim<quad>();

// Explicit routine instantiations for complex<quad>
//template void IO::sym_scheme<complex<quad>>();
template void IO::num_scheme<complex<quad>>();
template void IO::num_minim<complex<quad>>();
