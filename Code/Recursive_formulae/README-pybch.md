# pybch: the recursive formulae in Python

`pybch` is an open-source Python port of the derivation in `Efficient_Trotterizations.nb`, so you can reproduce the recursive formulae without Mathematica.
It covers schemes up to order $n = 8$ and gives exactly the same results as the notebook and as `NCoefficients.cpp` / `Prefactors.cpp` in [Trotter-Suzuki_schemes](../Trotter-Suzuki_schemes).

The only dependency is [SymPy](https://www.sympy.org) (`pip install sympy`), which is used only to format the results.
All other arithmetic is done with exact rational numbers from Python's standard library.

## Usage

Run these from this directory (`Code/Recursive_formulae`):

```bash
python -m pybch.prefactors                       # BCH prefactors of all orders (alpha, beta, gamma_k, delta_k, epsilon_k)
python -m pybch.prefactors --order 6             # only order 6
python -m pybch.prefactors --write               # also write pybch/data/prefactors_half.{json,m}

python -m pybch.recursion --order 4              # recursive formulae up to order 4 for an A step
python -m pybch.recursion --order 8 --step B     # B step, up to order 8
python -m pybch.recursion --order 6 --numeric    # with the prefactor values inserted

python -m unittest discover -s pybch/tests -t . -v   # tests (~10 s)
```

From Python:

```python
from pybch.recursion import recursive_formulas, to_sympy
formulas = to_sympy(recursive_formulas('A', order=8))   # {'gammas[0]': <sympy expression>, ...}
```

Variables use the names from the C++ code: `nu`, `sigma`, `alpha`, `beta`, `gammas[k]`, `deltas[k]`, `epsilons[k]` for the scheme coefficients (0-based indices), and `alpha_Q`, `beta_Q`, `gammas_Q[k]`, ... for the BCH prefactors.

## How it works

| Module | Content | Notebook counterpart |
|---|---|---|
| `freealg.py` | Exact arithmetic with non-commuting words in $A$, $B$; commutators; truncated `exp`, `log` and BCH products | FeynCalc `DotSimplify` / `Commutator` |
| `basis_defs.py` | The symmetric compact basis $C_k, D_k, E_k, F_k$ as right-nested commutators | Methods > Symmetric compact basis |
| `basis.py` | Coordinates of any Lie polynomial (odd degree $\le 9$) in that basis | Replacement rules, Hall word transformations |
| `bch.py`, `prefactors.py` | Prefactors of $\log(e^{A/2} e^{B} e^{A/2})$ in the basis | Symmetric compact basis coefficients |
| `recursion.py` | Recursive formulae for one A (or B) step | Recursive formula derivation |

1. **Basis coordinates.** A commutator never changes how many $A$'s and $B$'s a word contains.
   So a Lie polynomial splits into blocks of fixed multidegree $(n_A, n_B)$, and each block is expanded into words.
   For each block, Gaussian elimination picks as many "pivot" words as there are basis elements, such that the square matrix of their coefficients is invertible.
   The coordinates are then $M^{-1}$ applied to the coefficients of the pivot words.
   This also proves that the basis is linearly independent, and its size matches the dimension of the free Lie algebra (2, 6, 18, 56 for degrees 3, 5, 7, 9).
2. **Prefactors.** $\log(e^{A/2} e^{B} e^{A/2})$ is expanded directly as a series in words up to degree 9 and then projected onto the basis.
3. **Recursion.** One step maps $\Phi \to \log(e^{aA/2} e^{\Phi} e^{aA/2}) = aA + \Phi + \sum_K Q_K\, K(aA, \Phi)$.
   Each basis element $K$ is a right-nested commutator, so it is linear in each of its slots.
   Every $A$ slot becomes $aA$, and every $B$ slot is expanded over all components of $\Phi$ (up to the given order).
   Each nested commutator produced this way is expressed in the basis again.
   The prefactors $Q_K$ are kept symbolic, so each result is a rational linear combination of them, e.g.
   `gammas[0] += (gammas_Q[1] + gammas_Q[2]) a^3 nu sigma + ...`.

The tests compare the prefactors with the notebook values.
They also run a random multi-stage scheme in two ways (stage by stage with the recursive formulae, and through the BCH logarithm of the whole product) and require exact agreement.
Finally, they check that the B step is the $A \leftrightarrow B$ mirror of the A step, which is what `NScheme::step_B` relies on.
