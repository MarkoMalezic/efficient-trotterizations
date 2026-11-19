#include "InOut.h"

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

// Helper function to parse a list to extract vectors and arrays
// Specialization for array<double, 3>
template <>
array<double, 3> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<double, 3> list;
  for (int i = 0; i < 3; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<double, 2>");
    }
    list[i] = stod(trim(token));
  }
  return list;
}

// Specialization for array<long double, 3>
template <>
array<long double, 3> parse_list(const string &str)
{
  istringstream iss(str);
  string token;

  array<long double, 3> list;
  for (int i = 0; i < 3; ++i)
  {
    if (!getline(iss, token, ','))
    {
      throw runtime_error("Insufficient elements in list to parse array<long double, 2>");
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

  // Fill the VectorXcd with the extracted values
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

  // Fill the VectorXld with the extracted values
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

  // Fill the VectorXcld with the extracted values
  VectorXcld list(temp_list.size());
  for (size_t i = 0; i < temp_list.size(); ++i)
  {
    list(i) = temp_list[i];
  }
  return list;
}

// Helper function to generate a random vector
template <typename RealT>
Eigen::Matrix<RealT, Eigen::Dynamic, 1> generate_random_vector(int size, const RealT &min, const RealT &max, const int &seed)
{
    std::mt19937 gen(seed);  // Fixed seed for reproducibility
    std::uniform_real_distribution<RealT> dist(min, max);

    Eigen::Matrix<RealT, Eigen::Dynamic, 1> vec(size);
    for (int i = 0; i < size; ++i) {
        vec[i] = dist(gen);
    }
    return vec;
}

template <typename RealT>
void write_complex_ab_vec(const Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>& ab_vec, H5::H5File& file,
                  const std::string& name) {
    // Cast to double and split into real/imag
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ab_vec_double(ab_vec.size(), 2);
    for (int i = 0; i < ab_vec.size(); ++i) {
        ab_vec_double(i, 0) = static_cast<double>(real(ab_vec(i)));
        ab_vec_double(i, 1) = static_cast<double>(imag(ab_vec(i)));
    }
    // Create dataspace
    hsize_t dims[2] = {static_cast<hsize_t>(ab_vec.rows()), 2};
    H5::DataSpace dataspace(2, dims);

    // Create and write attribute
    H5::Attribute attr = file.createAttribute(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    attr.write(H5::PredType::NATIVE_DOUBLE, ab_vec_double.data());
}

// Helper function to write model metadata
template <typename RealT>
void write_model_metadata(H5::H5File &file, const Model<RealT> &model)
{
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE);

  file.createAttribute("model", strdatatype, H5::DataSpace()).write(strdatatype, &model.label);
  file.createAttribute("L", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &model.L);
  file.createAttribute("bound", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &model.bound);
  if (model.label == "ModelXZ")
  {
    auto xz_model = dynamic_cast<const ModelXZ<RealT> *>(&model);
    // Create the J attribute
    double J_double = static_cast<double>(xz_model->J);
    file.createAttribute("J", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &J_double);
    // Create the hmag attribute
    Eigen::Matrix<double, Eigen::Dynamic, 1> hmag_double = xz_model->hmag.template cast<double>();
    hsize_t hmag_dims[1] = {static_cast<hsize_t>(xz_model->hmag.size())}; // Store the size in a variable
    H5::Attribute hmag_attr = file.createAttribute("hmag", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, hmag_dims));
    hmag_attr.write(H5::PredType::NATIVE_DOUBLE, hmag_double.data());
  }
  else if (model.label == "ModelXXZ")
  {
    auto xxz_model = dynamic_cast<const ModelXXZ<RealT> *>(&model);
    // Create the J attribute
    double J_double = static_cast<double>(xxz_model->J);
    file.createAttribute("J", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &J_double);
    // Create the hmag attribute
    Eigen::Matrix<double, Eigen::Dynamic, 1> hmag_double = xxz_model->hmag.template cast<double>();
    hsize_t hmag_dims[1] = {static_cast<hsize_t>(xxz_model->hmag.size())}; // Store the size in a variable
    H5::Attribute hmag_attr = file.createAttribute("hmag", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, hmag_dims));
    hmag_attr.write(H5::PredType::NATIVE_DOUBLE, hmag_double.data());
  }
  else if (model.label == "ModelHeisenberg")
  {
    auto heisenberg_model = dynamic_cast<const ModelHeisenberg<RealT> *>(&model);
    // Create the J attribute
    array<double, 3> J_double = {static_cast<double>(heisenberg_model->J[0]), static_cast<double>(heisenberg_model->J[1]), static_cast<double>(heisenberg_model->J[2])};
    hsize_t dims[1] = {3}; // Define the dimensions of the dataspace
    H5::Attribute Jarr_attr = file.createAttribute("J", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dims));
    Jarr_attr.write(H5::PredType::NATIVE_DOUBLE, J_double.data());
    // Create the hmag attribute
    Eigen::Matrix<double, Eigen::Dynamic, 1> hmag_double = heisenberg_model->hmag.template cast<double>();
    hsize_t hmag_dims[1] = {static_cast<hsize_t>(heisenberg_model->hmag.size())}; // Store the size in a variable
    H5::Attribute hmag_attr = file.createAttribute("hmag", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, hmag_dims));
    hmag_attr.write(H5::PredType::NATIVE_DOUBLE, hmag_double.data());
  }
}

// Helper function to write a state to an hdf5 file
template <typename RealT>
void write_psi(const Operators::VectorCX<RealT> &psi, H5::H5File &file, const string &dataset_base, int step)
{
  const size_t N = psi.size();

  // Create flat buffers for real/imag
  vector<double> real_part(N);
  vector<double> imag_part(N);

  for (size_t i = 0; i < N; ++i)
  {
    // Cast to double if RealT is long double
    real_part[i] = static_cast<double>(psi[i].real());
    imag_part[i] = static_cast<double>(psi[i].imag());
  }

  // Define dataspace: 1 row (this time step), N columns
  hsize_t dims[2] = {1, static_cast<hsize_t>(N)};
  H5::DataSpace dataspace(2, dims);

  // Create dataset names
  string dset_real = dataset_base + "_real";
  string dset_imag = dataset_base + "_imag";

  // If first step, create datasets
  if (step == 0)
  {
    // Max dims for extendable dataset
    hsize_t maxdims[2] = {H5S_UNLIMITED, N};
    H5::DSetCreatPropList plist;
    hsize_t chunk_dims[2] = {1, N};
    plist.setChunk(2, chunk_dims);

    file.createDataSet(dset_real, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dims, maxdims), plist);
    file.createDataSet(dset_imag, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dims, maxdims), plist);
  }

  // Open datasets
  H5::DataSet dataset_real = file.openDataSet(dset_real);
  H5::DataSet dataset_imag = file.openDataSet(dset_imag);

  // Extend datasets to fit this step
  hsize_t new_dims[2];
  dataset_real.getSpace().getSimpleExtentDims(new_dims);
  new_dims[0] = step + 1;
  dataset_real.extend(new_dims);
  dataset_imag.extend(new_dims);

  // Select hyperslab to write this row
  H5::DataSpace filespace_real = dataset_real.getSpace();
  H5::DataSpace filespace_imag = dataset_imag.getSpace();
  hsize_t offset[2] = {static_cast<hsize_t>(step), 0};
  filespace_real.selectHyperslab(H5S_SELECT_SET, dims, offset);
  filespace_imag.selectHyperslab(H5S_SELECT_SET, dims, offset);

  // Memory space for this row
  H5::DataSpace memspace(2, dims);

  // Write the data
  dataset_real.write(real_part.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace_real);
  dataset_imag.write(imag_part.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace_imag);
}

// Constructor for InOut
template <typename Scalar>
InOut<Scalar>::InOut(const string read_file)
    : model(make_unique<Model<RealT>>(0, "", false))
{

  rfile.open(read_file);
  if (!rfile)
  {
    throw runtime_error("Error opening input file: " + read_file);
  }

  cout << "Reading input file: " << read_file << endl
       << string(50, '-') << endl;

  string line;
  while (getline(rfile, line) && line != "# General parameters")
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
    else if (key == "save_file")
    {
      save_file = value;
    }
    else if (key == "model")
    {
      cout << "Building model: " << value << endl;
      build_model(value);
      cout << "Model built" << endl
           << string(50, '-') << endl;
    }
  }
}

// Destructor for InOut
template <typename Scalar>
InOut<Scalar>::~InOut()
{
}

// Method to build the model object
template <typename Scalar>
void InOut<Scalar>::build_model(const string &model_str)
{
  using RealT = typename Eigen::NumTraits<Scalar>::Real;

  // Parse model parameters
  int L{0};
  bool bound{false};
  RealT J{0.0};
  array<RealT, 3> J_array{{0.0, 0.0, 0.0}};
  Eigen::Matrix<RealT, Eigen::Dynamic, 1> hmag;

  string line;
  while (getline(rfile, line) && line != "# Time evolution parameters")
  {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue; // Skip comments and empty lines

    size_t pos = line.find('=');
    if (pos == string::npos)
      continue; // Skip malformed lines

    string key = trim(line.substr(0, pos));
    string value = trim(line.substr(pos + 1));

    if (key == "L")
    {
      L = stoi(value);
    }
    else if (key == "bound")
    {
      bound = parse_bool(value);
    }
    else if (key == "J")
    {
      if (value.find(',') != string::npos)
      {
        // Parse J as an array
        J_array = parse_list<array<RealT, 3>>(value);
      }
      else
      {
        // Parse J as a single value
        J = parse_value<RealT>(value);
      }
    }
    else if (key == "hmag")
    {
      if (value == "random")
      {
        // Generate random magnetic field
        hmag = generate_random_vector<RealT>(L);
      }
      else
      {
        // Parse hmag as a vector
        hmag = parse_list<Eigen::Matrix<RealT, Eigen::Dynamic, 1>>(value);
      }
    }
  }

  // Dynamically allocate the appropriate model
  if (model_str == "ModelXZ")
  {
    cout << "J = " << J << endl;
    cout << "hmag = " << hmag.transpose() << endl;
    model = make_unique<ModelXZ<RealT>>(L, J, hmag, model_str, bound);
  }
  else if (model_str == "ModelXXZ")
  {
    cout << "J = " << J << endl;
    cout << "hmag = " << hmag.transpose() << endl;
    model = make_unique<ModelXXZ<RealT>>(L, J, hmag, model_str, bound);
  }
  else if (model_str == "ModelHeisenberg")
  {
    cout << "Jx = " << J_array[0] << "Jy = " << J_array[1] << "Jz = " << J_array[2] << endl;
    cout << "hmag =" << hmag.transpose() << endl;
    model = make_unique<ModelHeisenberg<RealT>>(L, J_array, hmag, model_str, bound);
  }
  else
  {
    throw invalid_argument("Unknown model type: " + model_str);
  }

  if (bound)
  {
    cout << "Open boundary conditions" << endl;
  }
  else
  {
    cout << "Periodic boundary conditions" << endl;
  }

  // Call the build_model() method
  model->build_model();
}

// *** Routines ***

// Time evolution
// Methods: - exact: time evolve a state using exact diagonalization
//          - trotter: time evolve a state using Trotter decomposition
template <typename Scalar>
void InOut<Scalar>::time_evolve()
{
  // Determine the Vec type from Scalar
  using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;

  // Specific parameters
  // Method to use for time evolution (exact, trotter)
  string method;

  // Initial state
  Operators::VectorCX<RealT> psi; // Initial state vector

  // Time evolution parameters
  RealT t{0.0}; // Time to evolve
  int N_steps{1}; // Time step

  // Trotter decomposition parameters
  int q{0};               // Number of cycles
  bool group{false};      // Grouping of stages (local = 0, global = 1)
  bool imag_time{false};  // Time evolution (imaginary time = 1, real time = 0)
  Vec a_vec;              // The a vector
  Vec b_vec;              // The b vector

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
    // Initial state
    else if (key == "psi")
    {
      if (value == "random")
      {
        // Generate random psi field
        psi = Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>::Random(pow(2, model->L));
        psi /= psi.norm(); // Normalize the state
      }
      else
      {
        // Parse psi as a vector
        psi = parse_list<Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>>(value);
      }
    }
    // Time evolution parameters
    else if (key == "t")
    {
      t = parse_value<RealT>(value);
    }
    else if (key == "N_steps")
    {
      N_steps = stoi(value);
    }
    // Trotter decomposition parameters
    else if (key == "q")
    {
      q = stoi(value);
    }
    else if (key == "group")
    {
      group = parse_bool(value);
    }
    else if (key == "imag_time")
    {
      imag_time = parse_bool(value);
    }
    else if (key == "a_vec")
    {
      a_vec = parse_list<Vec>(value); // Vector a in the symmetric basis
    }
    else if (key == "b_vec")
    {
      b_vec = parse_list<Vec>(value); // Vector b in the symmetric basis
    }
  }

  // Close the read_file
  rfile.close();

  RealT h = static_cast<RealT>(t / N_steps);

  cout << "Evolving state using method: " << method << endl;
  cout << "Time to evolve: " << t << endl;
  cout << "Number of time steps: " << N_steps << endl;
  cout << "Step size: " << h << endl;

  if (method == "exact")
  {
    // Construct the Exact diagonalization object
    Exact<RealT> diag(*model);

    // Prepare the hdf5 file for writing
    H5::H5File file(save_file, H5F_ACC_TRUNC);

    // Write the model metadata
    write_model_metadata(file, *model);

    H5::DataSpace scalar_space(H5S_SCALAR);
    // Write the Time evolution parameters
    double t_double = static_cast<double>(t);
    file.createAttribute("t", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &t_double);
    double h_double = static_cast<double>(h);
    file.createAttribute("h", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &h_double);
    file.createAttribute("imag_time", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &imag_time);

    write_psi(psi, file, "psi", 0);

    // Loop over the time steps and write the
    for (int i = 0; i < N_steps; ++i)
    {
      cout << "\rStep: " << i + 1 << " / " << N_steps << flush;
      // Evolve the state using the exact method
      if (imag_time)
      {
        psi = diag.imag_evolve(h, psi);
      }
      else
      {
        psi = diag.real_evolve(h, psi);
      }
      write_psi(psi, file, "psi", i + 1);
    }
    cout << endl
         << string(50, '-') << endl;
  }
  else if (method == "trotter")
  {
    // Construct the Trotter decomposition object
    Trotter<Vec> trotter(*model, q, group, a_vec, b_vec);

    // Prepare the hdf5 file for writing
    H5::H5File file(save_file, H5F_ACC_TRUNC);

    // Write the model metadata
    write_model_metadata(file, *model);

    H5::DataSpace scalar_space(H5S_SCALAR);
    // Write the Time evolution parameters
    double t_double = static_cast<double>(t);
    file.createAttribute("t", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &t_double);
    file.createAttribute("N_steps", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &N_steps);
    double h_double = static_cast<double>(h);
    file.createAttribute("h", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &h_double);
    file.createAttribute("imag_time", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &imag_time);

    // Write the Trotter parameters
    file.createAttribute("q", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &q);
    file.createAttribute("group", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &group);

    if constexpr (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
    {
      Eigen::Matrix<double, Eigen::Dynamic, 1> a_vec_double = a_vec.template cast<double>();
      hsize_t a_dim = a_vec_double.size();
      H5::DataSpace a_space(1, &a_dim);
      file.createAttribute("a_vec", H5::PredType::NATIVE_DOUBLE, a_space).write(H5::PredType::NATIVE_DOUBLE, a_vec_double.data());

      Eigen::Matrix<double, Eigen::Dynamic, 1> b_vec_double = b_vec.template cast<double>();
      hsize_t b_dim = b_vec_double.size();
      H5::DataSpace b_space(1, &b_dim);
      file.createAttribute("b_vec", H5::PredType::NATIVE_DOUBLE, b_space).write(H5::PredType::NATIVE_DOUBLE, b_vec_double.data());
    }
    else
    {
      write_complex_ab_vec(a_vec, file, "a_vec");
      write_complex_ab_vec(b_vec, file, "b_vec");
    }
    
    // Write the initial state
    write_psi(psi, file, "psi", 0);

    // Loop over the time steps and write the
    for (int i = 0; i < N_steps; ++i)
    {
      cout << "\rStep: " << i + 1 << " / " << N_steps << flush;
      // Evolve the state using the exact method
      psi = trotter.take_step(h, psi, imag_time);
      write_psi(psi, file, "psi", i + 1);
    }
    cout << endl
         << string(50, '-') << endl;
  }
  else
  {
    throw invalid_argument("Unknown method: " + method);
  }
}

// Error estimation
// Methods: - operator: compare the time evolution operators from exact diagonalization and Trotter decomposition
//          - state: compare the time evolved states from exact diagonalization and Trotter decomposition
template <typename Scalar>
void InOut<Scalar>::error_estimate()
{
  // Determine the Vec type from Scalar
  using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;

  // Specific parameters
  // Method to use for error estimation (operator, state)
  string method;

  // Time evolution parameters
  RealT t{0.0}; // Time to evolve
  int N_steps {1}; // Number of time steps

  // Trotter decomposition parameters
  int q{0};               // Number of cycles
  bool group{false};      // Grouping of stages (local = 0, global = 1)
  bool imag_time{false};  // Time evolution (imaginary time = 1, real time = 0)
  Vec a_vec;              // The a vector
  Vec b_vec;              // The b vector

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
    // Time evolution parameters
    else if (key == "t")
    {
      t = parse_value<RealT>(value);
    }
    else if (key == "N_steps")
    {
      N_steps = stoi(value);
    }
    // Trotter decomposition parameters
    else if (key == "q")
    {
      q = stoi(value);
    }
    else if (key == "group")
    {
      group = parse_bool(value);
    }
    else if (key == "imag_time")
    {
      imag_time = parse_bool(value);
    }
    else if (key == "a_vec")
    {
      a_vec = parse_list<Vec>(value); // Vector a in the symmetric basis
    }
    else if (key == "b_vec")
    {
      b_vec = parse_list<Vec>(value); // Vector b in the symmetric basis
    }
  }

  // Close the read_file
  rfile.close();

  RealT h = static_cast<RealT>(t / N_steps);
  cout << "Estimating error using method: " << method << endl;
  cout << "Time to evolve: " << t << endl;
  cout << "Number of time steps: " << N_steps << endl;
  cout << "Time step: " << h << endl;

  if (method == "operator")
  {
    // Construct the Exact diagonalization object
    Exact<RealT> diag(*model);

    // Construct the Trotter decomposition object
    Trotter<Vec> trotter(*model, q, group, a_vec, b_vec);

    // Prepare the hdf5 file for writing
    H5::H5File file(save_file, H5F_ACC_TRUNC);

    // Write the model metadata
    write_model_metadata(file, *model);

    H5::DataSpace scalar_space(H5S_SCALAR);
    // Write the Time evolution parameters
    double t_double = static_cast<double>(t);
    file.createAttribute("t", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &t_double);
    file.createAttribute("N_steps", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &N_steps);
    double h_double = static_cast<double>(h);
    file.createAttribute("h", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &h_double);
    file.createAttribute("imag_time", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &imag_time);

    // Write the Trotter parameters
    file.createAttribute("q", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &q);
    file.createAttribute("group", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &group);

    if constexpr (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
    {
      Eigen::Matrix<double, Eigen::Dynamic, 1> a_vec_double = a_vec.template cast<double>();
      hsize_t a_dim = a_vec_double.size();
      H5::DataSpace a_space(1, &a_dim);
      file.createAttribute("a_vec", H5::PredType::NATIVE_DOUBLE, a_space).write(H5::PredType::NATIVE_DOUBLE, a_vec_double.data());

      Eigen::Matrix<double, Eigen::Dynamic, 1> b_vec_double = b_vec.template cast<double>();
      hsize_t b_dim = b_vec_double.size();
      H5::DataSpace b_space(1, &b_dim);
      file.createAttribute("b_vec", H5::PredType::NATIVE_DOUBLE, b_space).write(H5::PredType::NATIVE_DOUBLE, b_vec_double.data());
    }
    else
    {
      write_complex_ab_vec(a_vec, file, "a_vec");
      write_complex_ab_vec(b_vec, file, "b_vec");
    }

    hsize_t dims[1] = {0}; // Start with an empty dataset
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hsize_t chunk_dims[1] = {1}; // Chunk size of 1 for appending one value at a time

    H5::DSetCreatPropList plist;
    plist.setChunk(1, chunk_dims);

    H5::DataSet error_dataset = file.createDataSet("error", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dims, maxdims), plist);

    const int N = static_cast<int>(pow(2, model->L));
    Operators::MatrixCX<RealT> U_exact;
    Operators::MatrixCX<RealT> U_trotter{Operators::MatrixCX<RealT>::Identity(N, N)};
    auto S_trotter = trotter.build_step_op(h, imag_time);

    double runtime{0.0};
    // Loop over the time steps and write the
    for (int i = 0; i < N_steps; ++i)
    {
      cout << "\rStep: " << i + 1 << " / " << N_steps << flush;

      auto start = chrono::high_resolution_clock::now();
      // Construct the exact time evolution operator
      if (imag_time)
      {
        U_exact = diag.imag_evolve_op(h * (i + 1));
      }
      else
      {
        U_exact = diag.real_evolve_op(h * (i + 1));
      }
      // Construct the Trotter time evolution operator
      U_trotter *= S_trotter;

      // Compare the two operators
      auto error = (U_exact - U_trotter).norm() / sqrt(N);
      double error_double = static_cast<double>(error);
      auto end = chrono::high_resolution_clock::now();
      runtime += chrono::duration<double>(end - start).count();

      cout << "\rError at step " << i + 1 << ": " << error << endl;

      // Extend the dataset to add one more row
      dims[0] = i + 1;
      error_dataset.extend(dims);

      // Select the hyperslab for the new row
      H5::DataSpace filespace = error_dataset.getSpace();
      hsize_t offset[1] = {static_cast<hsize_t>(i)};
      filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);

      // Write the error value
      H5::DataSpace memspace(1, chunk_dims);
      error_dataset.write(&error_double, H5::PredType::NATIVE_DOUBLE, memspace, filespace);
    }

    file.createAttribute("runtime", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &runtime);

    cout << endl
         << string(50, '-') << endl;
  }
  else if (method == "state")
  {
    // Construct the Exact diagonalization object
    Exact<RealT> diag(*model);

    // Construct the Trotter decomposition object
    Trotter<Vec> trotter(*model, q, group, a_vec, b_vec);

    // Prepare the hdf5 file for writing
    H5::H5File file(save_file, H5F_ACC_TRUNC);

    // Write the model metadata
    write_model_metadata(file, *model);

    H5::DataSpace scalar_space(H5S_SCALAR);
    // Write the Time evolution parameters
    double t_double = static_cast<double>(t);
    file.createAttribute("t", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &t_double);
    file.createAttribute("N_steps", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &N_steps);
    double h_double = static_cast<double>(h);
    file.createAttribute("h", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &h_double);
    file.createAttribute("imag_time", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &imag_time);

    // Write the Trotter parameters
    file.createAttribute("q", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &q);
    file.createAttribute("group", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &group);

    if constexpr (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
    {
      Eigen::Matrix<double, Eigen::Dynamic, 1> a_vec_double = a_vec.template cast<double>();
      hsize_t a_dim = a_vec_double.size();
      H5::DataSpace a_space(1, &a_dim);
      file.createAttribute("a_vec", H5::PredType::NATIVE_DOUBLE, a_space).write(H5::PredType::NATIVE_DOUBLE, a_vec_double.data());

      Eigen::Matrix<double, Eigen::Dynamic, 1> b_vec_double = b_vec.template cast<double>();
      hsize_t b_dim = b_vec_double.size();
      H5::DataSpace b_space(1, &b_dim);
      file.createAttribute("b_vec", H5::PredType::NATIVE_DOUBLE, b_space).write(H5::PredType::NATIVE_DOUBLE, b_vec_double.data());
    }
    else
    {
      write_complex_ab_vec(a_vec, file, "a_vec");
      write_complex_ab_vec(b_vec, file, "b_vec");
    }

    const int N = static_cast<int>(pow(2, model->L));

    hsize_t dims[1] = {static_cast<hsize_t>(N_steps)};
    H5::DataSet error_dataset = file.createDataSet("error", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dims));

    Operators::VectorCX<RealT> psi_exact;
    Operators::VectorCX<RealT> psi_trotter;

    // Vector to accumulate errors for each time step
    VectorXd errors(N_steps);
    errors.setZero();
    
    double runtime{0.0};

    if (imag_time)
    {
      // Initialize the states as random normalized vectors
      psi_exact = Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>::Random(pow(2, model->L));
      psi_exact /= psi_exact.norm();
      psi_trotter = Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>::Random(pow(2, model->L));
      psi_trotter /= psi_trotter.norm();

      // Start the clock
      auto start = chrono::high_resolution_clock::now();
      // Loop over the time steps
      for (int i = 0; i < N_steps; ++i)
      {
        cout << "\rStep: " << i + 1 << " / " << N_steps << flush;
        psi_exact = diag.imag_evolve(h, psi_exact);
        psi_trotter = trotter.take_step(h, psi_trotter, imag_time);

        // Add to the error at time step i
        errors(i) = static_cast<double>((psi_exact - psi_trotter).squaredNorm());
      }
      auto end = chrono::high_resolution_clock::now();
      runtime = chrono::duration<double>(end - start).count();
    }
    else
    {
      // Start the clock
      auto start = chrono::high_resolution_clock::now();
      // Loop over the basis states
      for (int v{0}; v < N; ++v)
      {
        cout << "\rState v = " << v + 1 << " / " << N << flush;
        psi_exact = Operators::VectorCX<RealT>::Zero(N);
        psi_exact(v) = 1.0;
        psi_trotter = Operators::VectorCX<RealT>::Zero(N);
        psi_trotter(v) = 1.0;

        // Loop over the time steps
        for (int i = 0; i < N_steps; ++i)
        {
          //cout << "\rStep: " << i + 1 << " / " << N_steps << flush;
          psi_exact = diag.real_evolve(h, psi_exact);
          psi_trotter = trotter.take_step(h, psi_trotter, imag_time);

          // Add to the error at time step i
          errors(i) += static_cast<double>((psi_exact - psi_trotter).squaredNorm());
        }
      }
      // Normalize the accumulated errors
      errors = (errors / N).array().sqrt();
      // Calculate the runtime
      auto end = chrono::high_resolution_clock::now();
      runtime = chrono::duration<double>(end - start).count();
    }
    
    // Write the errors to the dataset
    error_dataset.write(errors.data(), H5::PredType::NATIVE_DOUBLE);
    // Write the runtime
    file.createAttribute("runtime", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &runtime);

    cout << endl
         << string(50, '-') << endl;
  }
  else if (method == "state_total")
  {
    // Construct the Exact diagonalization object
    Exact<RealT> diag(*model);

    // Construct the Trotter decomposition object
    Trotter<Vec> trotter(*model, q, group, a_vec, b_vec);

    // Prepare the hdf5 file for writing
    H5::H5File file(save_file, H5F_ACC_TRUNC);

    // Write the model metadata
    write_model_metadata(file, *model);

    H5::DataSpace scalar_space(H5S_SCALAR);
    // Write the Time evolution parameters
    double t_double = static_cast<double>(t);
    file.createAttribute("t", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &t_double);
    file.createAttribute("N_steps", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &N_steps);
    double h_double = static_cast<double>(h);
    file.createAttribute("h", H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &h_double);
    file.createAttribute("imag_time", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &imag_time);

    // Write the Trotter parameters
    file.createAttribute("q", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &q);
    file.createAttribute("group", H5::PredType::NATIVE_HBOOL, H5::DataSpace()).write(H5::PredType::NATIVE_HBOOL, &group);

    if constexpr (is_same_v<Scalar, double> || is_same_v<Scalar, long double>)
    {
      Eigen::Matrix<double, Eigen::Dynamic, 1> a_vec_double = a_vec.template cast<double>();
      hsize_t a_dim = a_vec_double.size();
      H5::DataSpace a_space(1, &a_dim);
      file.createAttribute("a_vec", H5::PredType::NATIVE_DOUBLE, a_space).write(H5::PredType::NATIVE_DOUBLE, a_vec_double.data());

      Eigen::Matrix<double, Eigen::Dynamic, 1> b_vec_double = b_vec.template cast<double>();
      hsize_t b_dim = b_vec_double.size();
      H5::DataSpace b_space(1, &b_dim);
      file.createAttribute("b_vec", H5::PredType::NATIVE_DOUBLE, b_space).write(H5::PredType::NATIVE_DOUBLE, b_vec_double.data());
    }
    else
    {
      write_complex_ab_vec(a_vec, file, "a_vec");
      write_complex_ab_vec(b_vec, file, "b_vec");
    }

    const int N = static_cast<int>(pow(2, model->L));

    hsize_t dims[1] = {1};
    H5::DataSet error_dataset = file.createDataSet("error", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dims));

    Operators::VectorCX<RealT> psi_exact;
    Operators::VectorCX<RealT> psi_trotter;

    // The final Trotter error after the time evolution
    double error{0.0};
    double runtime{0.0};

    if (imag_time)
    {
      // Initialize the states as random normalized vectors
      psi_exact = Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>::Random(pow(2, model->L));
      psi_exact /= psi_exact.norm();
      psi_trotter = Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>::Random(pow(2, model->L));
      psi_trotter /= psi_trotter.norm();

      // Start the clock
      auto start = chrono::high_resolution_clock::now();

      // Not implemented

      auto end = chrono::high_resolution_clock::now();
      runtime = chrono::duration<double>(end - start).count();
    }
    else
    {
      // Start the clock
      auto start = chrono::high_resolution_clock::now();
      // Loop over the basis states
      for (int v{0}; v < N; ++v)
      {
        cout << "\rState v = " << v + 1 << " / " << N << flush;
        psi_exact = Operators::VectorCX<RealT>::Zero(N);
        psi_exact(v) = 1.0;
        psi_trotter = Operators::VectorCX<RealT>::Zero(N);
        psi_trotter(v) = 1.0;
        // Evolve the states
        psi_exact = diag.real_evolve(t, psi_exact);
        psi_trotter = trotter.evolve(t, N_steps, psi_trotter, imag_time);
        // Add to the total error
        error += static_cast<double>((psi_exact - psi_trotter).squaredNorm());
      }
      // Normalize the accumulated errors
      error = sqrt(error / N);
      cout << "\rTotal error over all states: " << error << endl;
      // Calculate the runtime
      auto end = chrono::high_resolution_clock::now();
      runtime = chrono::duration<double>(end - start).count();
    }

    // Write the error to the dataset
    error_dataset.write(&error, H5::PredType::NATIVE_DOUBLE);
    // Write the runtime
    file.createAttribute("runtime", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &runtime);

    cout << endl
         << string(50, '-') << endl;
  }
  else
  {
    throw invalid_argument("Unknown method: " + method);
  }
}

// Explicit routine instantiations for double
template class InOut<double>;
template class InOut<complex<double>>;
template class InOut<long double>;
template class InOut<complex<long double>>;
