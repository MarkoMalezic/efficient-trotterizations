"""Tests for the recursive formulae (orders 2..8).

Run from Code/Recursive_formulae:   python -m unittest discover -s pybch/tests -t . -v
"""
import random
import unittest
from fractions import Fraction

import sympy as sp

from pybch import freealg as fa
from pybch.basis import coordinates, mirror_name, mirror_sign, names_of_degree
from pybch.recursion import apply_step, recursive_formulas, to_sympy, var_name

ORDER = 8


def _mirror_var(v, inverse):
    if v in ('a', 'b'):
        return {'a': 'b', 'b': 'a'}[v]
    return var_name(mirror_name(inverse[v]))


class TestRecursion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.FA = recursive_formulas('A', ORDER)
        cls.FB = recursive_formulas('B', ORDER)

    def test_basis_mirror_symmetry(self):
        # NScheme::step_B relies on  A <-> B : K_k -> +K_{N+1-k}
        for d in (3, 5, 7, 9):
            for n in names_of_degree(d):
                self.assertEqual(mirror_sign(n), 1, n)

    def test_step_B_is_mirror_of_step_A(self):
        names = ['A', 'B'] + [n for d in range(3, ORDER + 2, 2) for n in names_of_degree(d)]
        inverse = {var_name(n): n for n in names}
        mirrored = {mirror_name(t): {tuple(sorted((_mirror_var(v, inverse), p) for v, p in m)): lin
                                     for m, lin in terms.items()}
                    for t, terms in self.FA.items()}
        self.assertEqual(mirrored, self.FB)

    def test_end_to_end_against_direct_bch(self):
        # 6-stage scheme with random rational parameters: recursion == log of the product
        rng = random.Random(7)
        steps = [('A' if i % 2 == 0 else 'B', Fraction(rng.randint(1, 9) * rng.choice([-1, 1]), rng.randint(1, 7)))
                 for i in range(6)]
        vals = {}
        for s, x in steps:
            vals = apply_step(vals, self.FA if s == 'A' else self.FB, x, s)
        factors = ([(fa.letter(s), x / 2) for s, x in reversed(steps[1:])]
                   + [(fa.letter(steps[0][0]), steps[0][1])]
                   + [(fa.letter(s), x / 2) for s, x in steps[1:]])
        direct = coordinates(fa.bch_product(factors, ORDER + 1), check=True)
        for k in set(direct) | set(vals):
            self.assertEqual(vals.get(k, 0), direct.get(k, 0), k)

    def test_known_formula_gamma1(self):
        # gammas[0] step A, as in the notebook and NCoefficients.cpp
        a, nu, sigma, alpha = sp.symbols('a nu sigma alpha')
        Q = sp.IndexedBase('gammas_Q')
        aQ, bQ = sp.symbols('alpha_Q beta_Q')
        expected = (Q[0] * a**4 * sigma + (Q[1] + Q[2]) * a**3 * nu * sigma
                    - (Q[3] + Q[4]) * a**2 * nu**2 * sigma - Q[5] * a * nu**3 * sigma
                    + aQ * a**2 * alpha - bQ * a * nu * alpha)
        self.assertEqual(sp.expand(to_sympy(self.FA)['gammas[0]'] - expected), 0)
