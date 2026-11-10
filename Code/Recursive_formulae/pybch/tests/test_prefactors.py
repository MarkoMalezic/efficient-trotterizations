"""Tests for the basis and the symmetric BCH prefactors.

Run from Code/Recursive_formulae:   python -m unittest discover -s pybch/tests -t . -v
"""
import unittest
from fractions import Fraction

import sympy as sp

from pybch import freealg as fa
from pybch.basis import coordinates, names_of_degree, poly
from pybch.bch import symmetric_bch_coefficients
from pybch.prefactors import prefactors

F = Fraction
# Values from Trotterization.nb (Methods > Symmetric compact basis > Coefficients),
# identical to src/Prefactors.cpp.
NB = {
    2: [F(-1, 24), F(1, 12)],
    4: [F(7, 5760), F(1, 480), F(1, 360), F(1, 360), F(-1, 120), F(-1, 720)],
    6: [F(-31, 967680), F(-19, 80640), F(23, 161280), F(-1, 10080), F(-1, 10080), F(-1, 20160),
        F(13, 60480), F(-1, 5040), F(-1, 3360), F(-19, 40320), F(1, 6720), F(-1, 7560), F(1, 1008),
        F(-1, 10080), F(-1, 10080), F(1, 10080), F(1, 5040), F(1, 30240)],
    8: [F(127, 154828800), F(11, 774144), F(-41, 3225600), F(17, 9676800), F(1, 302400),
        F(1, 290304), F(-1, 48384), F(1, 26880), F(-1, 13824), F(1, 15360), F(-23, 1612800),
        F(211, 4838400), F(-37, 1209600), F(7, 691200), F(47, 1209600), F(-281, 2419200),
        F(163, 1612800), F(73, 1612800), F(-43, 2419200), F(-79, 1209600), F(-1, 193536),
        F(-13, 3225600), F(-1, 25200), F(1, 11200), F(-1, 75600), F(1, 12600), F(-1, 25200),
        F(-1, 67200), F(-1, 67200), F(-31, 403200), F(11, 201600), F(-1, 1209600), F(17, 134400),
        F(-1, 19200), F(31, 604800), F(-1, 60480), F(-1, 50400), F(-1, 5400), F(1, 12600),
        F(1, 12600), F(-11, 151200), F(1, 100800), F(-29, 1209600), F(83, 1209600),
        F(-31, 1209600), F(1, 50400), F(-1, 10080), F(1, 30240), F(-1, 20160), F(1, 15120),
        F(-1, 90720), F(1, 302400), F(-1, 75600), F(1, 50400), F(-1, 60480), F(-1, 1209600)],
}
# log(exp(A) exp(B) exp(A)), order 8, F1..F5 (notebook output)
NB_FULL_ORDER8 = [F(127, 604800), F(11, 6048), F(-41, 25200), F(17, 75600), F(2, 4725)]
# Witt dimensions of the free Lie algebra on 2 letters
WITT = {3: 2, 5: 6, 7: 18, 9: 56}


class TestBasis(unittest.TestCase):
    def test_dimensions(self):
        for d, n in WITT.items():
            self.assertEqual(len(names_of_degree(d)), n)

    def test_basis_is_invertible(self):
        # coordinates() builds and inverts every multidegree block; reconstruct each element
        for d in WITT:
            for n in names_of_degree(d):
                self.assertEqual(coordinates(poly(n), check=True), {n: 1})



class TestPrefactors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.Q = symmetric_bch_coefficients(9, 'half')

    def test_orders_2_to_8_match_notebook(self):
        for order, vals in NB.items():
            names = names_of_degree(order + 1)
            self.assertEqual([self.Q[n] for n in names], vals, f'order {order}')

    def test_full_scheme_order8(self):
        Q = symmetric_bch_coefficients(9, 'full')
        self.assertEqual([Q[f'F{k}'] for k in range(1, 6)], NB_FULL_ORDER8)

    def test_sympy_output(self):
        pref = prefactors(8)
        self.assertEqual([len(pref[n]) for n in (2, 4, 6, 8)], [2, 6, 18, 56])
        self.assertIsInstance(pref[8][0], sp.Rational)
        self.assertEqual(pref[2], [sp.Rational(-1, 24), sp.Rational(1, 12)])


class TestAgainstSympy(unittest.TestCase):
    """Independent check of the word-expansion BCH with sympy's noncommutative algebra."""

    def test_bch_degree5(self):
        A, B = sp.symbols('A B', commutative=False)
        N = 5

        def _deg(t):
            return sum(p for b, p in (f.as_base_exp() for f in sp.Mul.make_args(t)) if b in (A, B))

        def trunc(e):
            return sp.Add(*[t for t in sp.Add.make_args(sp.expand(e)) if _deg(t) <= N])

        def expo(X):
            r, term = sp.S(1), sp.S(1)
            for n in range(1, N + 1):
                term = trunc(term * X / n)
                r += term
            return r

        Z = trunc(expo(A / 2) * expo(B) * expo(A / 2))
        Y = trunc(Z - 1)
        L, P = sp.S(0), sp.S(1)
        for k in range(1, N + 1):
            P = trunc(P * Y)
            L += sp.Rational((-1) ** (k + 1), k) * P
        L = sp.expand(L)
        ours = fa.bch_product([(fa.letter('A'), F(1, 2)), (fa.letter('B'), 1),
                               (fa.letter('A'), F(1, 2))], N)
        theirs = {}
        for t in sp.Add.make_args(L):
            c, rest = t.as_coeff_Mul()
            word = ''.join(str(b) * int(p) for b, p in (f.as_base_exp() for f in sp.Mul.make_args(rest)))
            theirs[word] = theirs.get(word, 0) + F(int(c.p), int(c.q))
        theirs = {w: v for w, v in theirs.items() if v}
        self.assertEqual(theirs, ours)


if __name__ == '__main__':
    unittest.main()
