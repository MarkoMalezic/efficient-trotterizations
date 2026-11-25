[![DOI](https://zenodo.org/badge/1135085602.svg)](https://doi.org/10.5281/zenodo.18347429)

# Efficient High Order Trotterizations

A derivation, implementation and application of Trotter-Suzuki decomposition schemes at high orders.

This repository contains the code required to reproduce the results presented in our papers:

- M. Maležič and J. Ostmeyer, "Efficient Trotter-Suzuki Schemes for Long-time Quantum Dynamics", J. Phys. A, [10.1088/1751-8121/ae919f](https://doi.org/10.1088/1751-8121/ae919f) ([arXiv:2601.18756](https://arxiv.org/abs/2601.18756))
- M. Maležič and J. Ostmeyer, "Optimised 8th-Order Trotter-Suzuki Schemes for Quantum Simulation" ([arXiv:2611.27896](https://arxiv.org/abs/2611.27896))

The algorithm used to implement such Omelyan-type schemes as quantum circuits (see [Quantum circuits](#quantum-circuits)) is described in our proceedings:

- M. Maležič and J. Ostmeyer, "Reducing the Gate Count with Efficient Trotter-Suzuki Schemes", Proceedings of the 42nd International Symposium on Lattice Field Theory (Lattice 2025) ([arXiv:2602.21145](https://arxiv.org/abs/2602.21145))

It also stores novel schemes at order $n = 2, 4, 6, 8$ in different formats, and a selection of recommended schemes.
Data for time evolution and the quantum circuit simulations is not included in the repository, but it will be gladly provided upon request.

For questions concerning the code contact: malezic@hiskp.uni-bonn.de

## Derivation

The main part of the derivations for Trotterizations are the recursive formulae found in a Mathematica notebook `Efficient_Trotterizations.nb` within [Recursive_formulae](Code/Recursive_formulae).
In the same folder we gather the results of the recursive formulae in a slightly more readable form `Recursive_formulae.pdf`.
For those without access to Mathematica, the same derivation (up to order $n = 8$) is available as an open-source Python package `pybch` in the same folder, see `README-pybch.md`.

Applying the formulae it is possible to construct polynomial manifolds, which are to be minimised in order to find efficient schemes.
Computation of such manifolds for 2nd and 4th order and their visualisations can be found in [Manifolds](Code/Manifolds) directory, both as Mathematica notebooks and as standalone PDF images.
Furthermore, properties of the Jacobian on the manifolds at orders $n = 6, 8$ give insight into the structure of the over-constrained regions, see the notebook `Jacobian_on_manifold.nb` in the same directory.

## Implementation

Polynomial manifolds constructed from the recursive formulae need to be minimised in order to find efficient Trotterizations.
The code for this is written in C++ (see folder [Trotter-Suzuki_schemes](Code/Trotter-Suzuki_schemes)), which can be built using CMake (look for `build-run_instructions.txt`).
The main dependency is on a linear algebra library called Eigen (found at https://eigen.tuxfamily.org/), while other modules can be found within the standard C++ library.
It also needs the header-only Boost.Multiprecision library.
By default only real double and long double parameters are compiled; quad precision and complex parameters are switched on with the CMake options `USE_QUAD` and `USE_COMPLEX`, where quad precision additionally requires GCC and `libquadmath`.

There are multiple routines, which are first split between the symbolic and numeric implementation of the formulae.
Furthermore, the routines split between scheme evaluation at specific parameters and their minimisation.
The minimiser is a Levenberg-Marquardt algorithm, optionally with geodesic acceleration and restricted to compositions of second-order schemes, and it supports schemes up to order $n = 8$ in double, long double and quad precision.
The program can be run by `/path_to_build/main input_file.in`, where the `input_file.in` specifies the routine and its parameters.
More details on specific subroutines and their usage can be found in `README-Trotter.md` within [Trotter-Suzuki_schemes](Code/Trotter-Suzuki_schemes).

In our minimisation efforts we were able to find many Trotter-Suzuki schemes, that are gathered inside [Schemes](/Code/Trotter-Suzuki_schemes/Schemes) in a few formats/notations, which are explained in `README-Schemes.md`.
Besides the archive of all schemes found, we also provide a selection of recommended schemes at each number of cycles, chosen for their theoretical efficiency and uniformity.

## Application

The schemes found by the minimiser in the implementation can then be used for time evolution of physical models.
Inside [Time_Evolution](Code/Time_Evolution) we model two distinct classes of quantum mechanical systems.

Within `Symplectic_evolution.py` one finds a Python implementation for time evolution of symplectic models.
The code currently includes a class Harmonic, which models the quantum harmonic oscillator.
Time evolution can be done exactly and using arbitrary Trotter-Suzuki schemes, which allows comparison.
To run this code one needs to run: `python Symplectic_evolution.py input_file.in`, where an input file needs to be set up according to the examples.
More details inside `README-Symplectic.md`.

A more complicated implementation of time evolution for the Heisenberg model and its derivatives can also be found in this directory. 
The code for this is again written in C++ with a dependency on Eigen and the building/running works similarly to before (look for `build-run_instructions.txt`).
The time evolution can be evaluated by using exact diagonalisation or Trotterizations, which is again useful for comparison and evaluation of the Trotter error.
There also exist a few routines, and which ones are run can be specified in the input file.
Further details can be found in `README-Heisenberg.md`.

### Quantum circuits

Inside [Quantum_circuits](Code/Quantum_circuits) we implement the schemes as quantum circuits for the Heisenberg model and simulate them with and without gate noise.
The circuits are built with our Qiskit plugin [qiskit-omelyan](https://github.com/MarkoMalezic/qiskit-omelyan) (version 1.0.0 or later, [10.5281/zenodo.24371684](https://doi.org/10.5281/zenodo.24371684)), which implements the algorithm described in [arXiv:2602.21145](https://arxiv.org/abs/2602.21145), and simulated with Qiskit Aer.
To run this code one needs to run: `python noisy_circuits.py`, where the model, scheme and noise are set inside the script.
More details inside `README-Circuits.md`.

### Error evolution animation

Combining the theoretical error computed from schemes in the folder [Schemes](/Code/Trotter-Suzuki_schemes/Schemes) with the distance from the origin term as described in the paper, and plotting it against the experimental error, we can animate how the correlation improves if the extra term is added.
These animations can be found in [Error_comparison_animations](/Error_comparison_animations/) for each cycle and below for schemes at 14 cycles:

https://github.com/user-attachments/assets/2b385df9-96ee-408e-b8ad-27ea129d2447
