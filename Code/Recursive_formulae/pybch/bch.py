"""Symmetric BCH coefficients in the compact basis.

Replaces the notebook chapter "Symmetric compact basis coefficients": instead of
expanding with the Arnal et al. algorithm (PhiCmt) and identifying each right-nested
commutator by a replacement rule, we expand the BCH series in words and solve for
the coordinates exactly.
"""
from fractions import Fraction

from . import freealg as fa
from .basis import coordinates, names_of_degree
from .basis_defs import MAX_DEGREE

# Factor sequences for the two conventions used in the notebook.
SCHEMES = {
    # log(exp(A/2) exp(B) exp(A/2)):  aQ, bQ, gQ, dQ, eQ  (notebook convention)
    'half': [('A', Fraction(1, 2)), ('B', 1), ('A', Fraction(1, 2))],
    # log(exp(A) exp(B) exp(A))
    'full': [('A', 1), ('B', 1), ('A', 1)],
}


def bch_series(scheme='half', maxdeg=MAX_DEGREE):
    factors = [(fa.letter(x), c) for x, c in SCHEMES[scheme]]
    return fa.bch_product(factors, maxdeg)


def symmetric_bch_coefficients(maxdeg=MAX_DEGREE, scheme='half'):
    """dict basis name -> Fraction for all odd degrees 1..maxdeg (zero entries included)."""
    Z = bch_series(scheme, maxdeg)
    coords = coordinates(Z)
    out = {}
    for d in range(1, maxdeg + 1, 2):
        for n in names_of_degree(d):
            out[n] = coords.get(n, Fraction(0))
    return out
