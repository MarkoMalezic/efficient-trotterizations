#include <chrono>
#include "IO.h"

// Runs the routine given by the input file (see main below for the error handling)
static int run(int argc, char *argv[])
{
  // Check if input file is provided
  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1; // Exit with error if argument is missing
  }

  auto start = chrono::high_resolution_clock::now();

  IO input(argv[1]);

  // Precisions that exist but were not compiled into this build (see CMake options USE_QUAD, USE_COMPLEX)
#if !defined(USE_COMPLEX)
  if (input.scalar_type.rfind("complex_", 0) == 0)
  {
    cerr << "Scalar type '" << input.scalar_type << "' is not available in this build: "
         << "reconfigure with  cmake -DUSE_COMPLEX=ON .." << endl;
    return 1;
  }
#endif
#if !defined(USE_QUAD)
  if (input.scalar_type == "quad" || input.scalar_type == "complex_quad")
  {
    cerr << "Scalar type '" << input.scalar_type << "' is not available in this build: "
         << "reconfigure with  cmake -DUSE_QUAD=ON .." << endl;
    return 1;
  }
#endif

  if (input.scalar_type == "double")
  {
    switch (input.routine)
    {
    case 0:
    {
      if (input.mode == 0)
      {
        input.num_minim<double>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic minimization is no longer supported." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    case 1:
    {
      if (input.mode == 0)
      {
        input.num_scheme<double>();
      }
      else if (input.mode == 1)
      {
        input.sym_scheme<double>();
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
#if defined(USE_COMPLEX)
  else if (input.scalar_type == "complex_double")
  {
    switch (input.routine)
    {
    case 0:
    {
      if (input.mode == 0)
      {
        input.num_minim<complex<double>>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic minimization is no longer supported." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    case 1:
    {
      if (input.mode == 0)
      {
        input.num_scheme<complex<double>>();
      }
      else if (input.mode == 1)
      {
        input.sym_scheme<complex<double>>();
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
#endif // defined(USE_COMPLEX)
  else if (input.scalar_type == "long_double")
  {
    switch (input.routine)
    {
    case 0:
    {
      if (input.mode == 0)
      {
        input.num_minim<long double>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic minimization is no longer supported." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    case 1:
    {
      if (input.mode == 0)
      {
        input.num_scheme<long double>();
      }
      else if (input.mode == 1)
      {
        input.sym_scheme<long double>();
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
#if defined(USE_COMPLEX)
  else if (input.scalar_type == "complex_long_double")
  {
    switch (input.routine)
    {
    case 0:
    {
      if (input.mode == 0)
      {
        input.num_minim<complex<long double>>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic minimization is no longer supported." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    case 1:
    {
      if (input.mode == 0)
      {
        input.num_scheme<complex<long double>>();
      }
      else if (input.mode == 1)
      {
        input.sym_scheme<complex<long double>>();
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
#endif // defined(USE_COMPLEX)
#if defined(USE_QUAD)
  else if (input.scalar_type == "quad")
  {
    switch (input.routine)
    {
    case 0:
    {
      if (input.mode == 0)
      {
        input.num_minim<quad>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic minimization for quad is not implemented yet." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    case 1:
    {
      if (input.mode == 0)
      {
        input.num_scheme<quad>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic scheme for quad is not implemented yet." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
#endif // defined(USE_QUAD)
#if defined(USE_QUAD) && defined(USE_COMPLEX)
  else if (input.scalar_type == "complex_quad")
  {
    switch (input.routine)
    {
    case 0:
    {
      if (input.mode == 0)
      {
        input.num_minim<complex<quad>>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic minimization for complex<quad> is not implemented yet." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    case 1:
    {
      if (input.mode == 0)
      {
        input.num_scheme<complex<quad>>();
      }
      else if (input.mode == 1)
      {
        cout << "Symbolic scheme for complex<quad> is not implemented yet." << endl;
      }
      else
      {
        throw runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      throw runtime_error("Invalid routine!");
    }
    }
  }
#endif // defined(USE_QUAD) && defined(USE_COMPLEX)
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
