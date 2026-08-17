# Efficient High Order Trotterizations

A derivation, implementation and application of Trotter-Suzuki decomposition schemes at high orders.

This repository contains the code required to reproduce the results presented in "Efficient Trotter-Suzuki Schemes for Long-time Quantum Dynamics".
It also stores novel schemes at order n = 2, 4, 6, 8 in different formats.
Data for time evolution is not included in the repository, but it will be gladly provided upon request.

For questions concerning the code contact: malezic@hiskp.uni-bonn.de

## Derivation

The main part of the derivations for Trotterizations are the recursive formulae found in a Mathematica notebook `Efficient_Trotterizations.nb` within [Recursive_formulae](Code/Recursive_formulae).
In the same folder we gather the results of the recursive formulae in a slightly more readable form `Recursive_formulae.pdf`.

Applying the formulae it is possible to construct polynomial manifolds, which are to be minimized in order to find efficient schemes.
Computation of such manifolds for 2nd and 4th order and their visualizations can be found in [Manifolds](Code/Manifolds) directory, both as Mathematica notebooks and as standalone PDF images.

## Implementation

Polynomial manifolds constructed from the recursive formulae need to be minimized in order to find efficient Trotterizations.
The code for this is written in C++ (see folder [Trotter-Suzuki_schemes](Code/Trotter-Suzuki_schemes)), which can be built using CMake (look for `build-run_instructions.txt`).
The main dependency is on a linear algebra library called Eigen (found at https://eigen.tuxfamily.org/), while other modules can be found within the standard C++ library.

There are multiple routines, which are first split between the symbolic and numeric implementation of the formulae.
Furthermore, the routines split between scheme evaluation at specific parameters and their minimization.
The program can be run by `/path_to_build/main input_file.in`, where the `input_file.in` specifies the routine and its parameters.
More details on specific subroutines and their usage can be found in `README-Trotter.md` within [Trotter-Suzuki_schemes](Code/Trotter-Suzuki_schemes).

In our minimization efforts we were able to find many Trotter-Suzuki schemes, that are gathered inside [Schemes](/Code/Trotter-Suzuki_schemes/Schemes) in a few formats/notations, which are explained in `README-Schemes.md`

## Application

The schemes found by the minimizer in the implementation can then be used for time evolution of physical models.
Inside [Time_Evolution](Code/Time_Evolution) we model two distinct classes of quantum mechanical systems.

Within `Symplectic_evolution.py` one finds a Python implementation for time evolution of symplectic models.
The code currently includes a class Harmonic, which models the quantum harmonic oscillator.
Time evolution can be done exactly and using arbitrary Trotter-Suzuki schemes, which allows comparison.
To run this code one needs to run: `python Symplectic_evolution.py input_file.in`, where an input file needs to be set up according to the examples.
More details inside `README-Symplectic.md`.

A more complicated implementation of time evolution for the Heisenberg model and its derivatives can also be found in this directory. 
The code for this is again written in C++ with a dependency on Eigen and the building/running works similarly to before (look for `build-run_instructions.txt`).
The time evolution can be evaluated by using exact diagonalization or Trotterizations, which is again useful for comparison and evaluation of the Trotter error.
There also exist a few routines, and which ones are run can be specified in the input file.
Further details can be found in `README-Heisenberg.md`.

### Error evolution animation

Combining the theoretical error computed from schemes in the folder [Schemes](/Code/Trotter-Suzuki_schemes/Data/Schemes) with the distance from the origin term as described in the paper, and plotting it against the experimental error, we can animate how the correlation improves if the extra term is added.
These animations can be found in [Error_comparison_animations](/Error_comparison_animations/) for each cycle and below for schemes at 14 cycles:

https://github.com/user-attachments/assets/2b385df9-96ee-408e-b8ad-27ea129d2447
