#include <chrono>
#include "IO.h"

int main(int argc, char *argv[])
{
  // Check if input file is provided
  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1; // Exit with error if argument is missing
  }

  auto start = chrono::high_resolution_clock::now();

  IO input(argv[1]);

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
        input.sym_minim<double>();
      }
      else
      {
        runtime_error("Invalid mode!");
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
        runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      runtime_error("Invalid routine!");
    }
    }
  }
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
        input.sym_minim<complex<double>>();
      }
      else
      {
        runtime_error("Invalid mode!");
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
        runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      runtime_error("Invalid routine!");
    }
    }
  }
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
        input.sym_minim<long double>();
      }
      else
      {
        runtime_error("Invalid mode!");
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
        runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      runtime_error("Invalid routine!");
    }
    }
  }
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
        input.sym_minim<complex<long double>>();
      }
      else
      {
        runtime_error("Invalid mode!");
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
        runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      runtime_error("Invalid routine!");
    }
    }
  }
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
        runtime_error("Invalid mode!");
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
        runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      runtime_error("Invalid routine!");
    }
    }
  }
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
        runtime_error("Invalid mode!");
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
        runtime_error("Invalid mode!");
      }
      break;
    }
    default:
    {
      runtime_error("Invalid routine!");
    }
    }
  }
  else
  {
    runtime_error("Invalid scalar type!");
  }

  auto end = chrono::high_resolution_clock::now();
  chrono::duration<double> duration = end - start;
  cout << "Time taken: " << fixed << setprecision(6) << duration.count() << " seconds" << std::endl;

  return 0;
}
