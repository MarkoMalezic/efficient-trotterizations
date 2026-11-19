#include "InOut.h"

using Eigen::VectorXcd;
using VectorXcld = Eigen::Matrix<complex<long double>, Eigen::Dynamic, 1>;
using Eigen::VectorXd;
using VectorXld = Eigen::Matrix<long double, Eigen::Dynamic, 1>;

using Eigen::MatrixXcd;
using MatrixXcld = Eigen::Matrix<complex<long double>, Eigen::Dynamic, Eigen::Dynamic>;
using Eigen::MatrixXd;
using MatrixXld = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;

// Runs the routine given by the input file (see main below for the error handling)
static int run(int argc, char *argv[])
{
  // Check if input file is provided
  if (argc != 3)
  {
    cerr << "Usage: " << argv[0] << " <input_file> <scalar_type>" << std::endl;
    return 1; // Exit with error if argument is missing
  }

  auto start = chrono::high_resolution_clock::now();

  if (string(argv[2]) == "double")
  {
    InOut<double> input(argv[1]);
    switch (input.routine)
    {
    case 0:
    {
      input.time_evolve();
      break;
    }
    case 1:
    {
      input.error_estimate();
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
  else if (string(argv[2]) == "long_double")
  {
    InOut<long double> input(argv[1]);
    switch (input.routine)
    {
    case 0:
    {
      input.time_evolve();
      break;
    }
    case 1:
    {
      input.error_estimate();
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
  else if (string(argv[2]) == "complex_double")
  {
    InOut<complex<double>> input(argv[1]);
    switch (input.routine)
    {
    case 0:
    {
      input.time_evolve();
      break;
    }
    case 1:
    {
      input.error_estimate();
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
  else if (string(argv[2]) == "complex_long_double")
  {
    InOut<complex<long double>> input(argv[1]);
    switch (input.routine)
    {
    case 0:
    {
      input.time_evolve();
      break;
    }
    case 1:
    {
      input.error_estimate();
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
  else
  {
    throw runtime_error("Invalid scalar type!");
  }

  auto end = chrono::high_resolution_clock::now();
  chrono::duration<double> duration = end - start;
  cout << "Time taken: " << fixed << setprecision(6) << duration.count() << " seconds" << std::endl;

  return 0;
}

int main(int argc, char *argv[])
{
  // Print errors (e.g. an invalid input file) as a message instead of terminating with an uncaught exception
  try
  {
    return run(argc, argv);
  }
  catch (const exception &e)
  {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }
}
