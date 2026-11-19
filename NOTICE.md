Third-Party Software and Notices
================================

This software makes use of the following third-party components:

Eigen
-----
This project uses the Eigen C++ template library.
Eigen is available at https://eigen.tuxfamily.org/

Eigen is licensed under the Mozilla Public License version 2.0 (MPL-2.0).
A copy of the license is available at https://www.mozilla.org/MPL/2.0/.

Wolfram Mathematica
-------------------
This repository includes scripts written in the Wolfram Language (Mathematica).
Mathematica is available at https://www.wolfram.com/mathematica/

These scripts are licensed under the same license as the rest of the
repository (MIT License). Running them requires a licensed copy of Wolfram Mathematica.
No rights to Mathematica itself are granted by this repository.

SymPy
-----
The Python package in Code/Recursive_formulae/pybch uses SymPy for symbolic output.
SymPy is available at https://www.sympy.org/ and is licensed under the BSD 3-Clause License.

Boost.Multiprecision
--------------------
The C++ program in Code/Trotter-Suzuki_schemes uses the Boost.Multiprecision library (https://www.boost.org/)
for quad precision. Boost is licensed under the Boost Software License, Version 1.0.

GCC libquadmath
---------------
When compiled with quad precision (CMake option USE_QUAD), the C++ program in Code/Trotter-Suzuki_schemes
is linked against libquadmath.
libquadmath is part of GCC (https://gcc.gnu.org/) and is licensed under the
GNU Lesser General Public License (LGPL) version 2.1 or later.

Qiskit, Qiskit Aer and qiskit-omelyan
-------------------------------------
The scripts in Code/Quantum_circuits use Qiskit (https://github.com/Qiskit/qiskit),
Qiskit Aer (https://github.com/Qiskit/qiskit-aer) and the qiskit-omelyan plugin
(https://github.com/MarkoMalezic/qiskit-omelyan), all licensed under the Apache License 2.0.
Optionally, they use Mitiq (https://github.com/unitaryfund/mitiq), licensed under the GNU General Public License v3.0.

NumPy and SciPy
---------------
The Python code uses NumPy (https://numpy.org/) and SciPy (https://scipy.org/),
both licensed under the BSD 3-Clause License.
