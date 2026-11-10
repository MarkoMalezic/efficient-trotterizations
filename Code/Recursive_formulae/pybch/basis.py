"""The symmetric compact basis (A, B, C_k, D_k, E_k, F_k, G_k) and coordinates in it.

`coordinates(p)` expresses any Lie polynomial p of odd degree <= 9 in the basis.
This is what the notebook does by hand with replacement rules (chapter "Symmetric compact
basis coefficients") and with TransformWord (chapter "Hall word transformations").
"""
from fractions import Fraction
from functools import lru_cache

from . import freealg as fa
from .basis_defs import LEVELS, MAX_DEGREE, RIGHT_NESTED

LEVEL_OF_DEGREE = {1: None, **{d: L for L, d in LEVELS.items()}}


ALL_NAMES = ['A', 'B'] + sorted(RIGHT_NESTED, key=lambda n: (LEVELS[n[0]], int(n[1:])))


def degree(name):
    return 1 if name in ('A', 'B') else LEVELS[name[0]]


def names_of_degree(d):
    if d == 1:
        return ['A', 'B']
    return [n for n in ALL_NAMES if n not in ('A', 'B') and LEVELS[n[0]] == d]


def word_of(name):
    """(sign, right-nested word) of a basis element."""
    if name in ('A', 'B'):
        return 1, name
    return RIGHT_NESTED[name]


@lru_cache(maxsize=None)
def _poly_cached(name):
    s, w = word_of(name)
    return tuple(fa.scale(fa.right_nested_word(w), s).items())


def poly(name):
    """Basis element expanded in the free associative algebra (words in A, B)."""
    return dict(_poly_cached(name))


def multidegree(name):
    _, w = word_of(name)
    return w.count('A'), w.count('B')


class _Block:
    """Basis elements of one multidegree (nA, nB) with a precomputed left inverse.

    We pick a set of 'pivot' words such that the square matrix M[pivot, basis] is
    invertible; then coords(p) = M^-1 p[pivot] for any p in the span.
    `extra` are additional spanning polynomials (e.g. of an ideal) whose coordinates
    are computed but discarded.
    """

    def __init__(self, names, extra=()):
        self.names = list(names) + [None] * len(extra)
        polys = [poly(n) for n in names] + list(extra)
        names = self.names
        words = sorted(set(w for p in polys for w in p))
        n = len(names)
        # Row-reduce the transposed matrix (basis x words) to find independent words.
        rows = [[p.get(w, Fraction(0)) for w in words] for p in polys]
        piv_cols = []
        r = 0
        work = [row[:] for row in rows]
        for c in range(len(words)):
            pr = next((i for i in range(r, n) if work[i][c] != 0), None)
            if pr is None:
                continue
            work[r], work[pr] = work[pr], work[r]
            inv = 1 / work[r][c]
            work[r] = [x * inv for x in work[r]]
            for i in range(n):
                if i != r and work[i][c] != 0:
                    f = work[i][c]
                    work[i] = [a - f * b for a, b in zip(work[i], work[r])]
            piv_cols.append(c)
            r += 1
            if r == n:
                break
        if r < n:
            raise ValueError(f'basis elements {names} are linearly dependent (rank {r} < {n})')
        self.pivots = [words[c] for c in piv_cols]
        # M[i][j] = coefficient of pivot word i in basis element j ; invert it.
        M = [[polys[j].get(w, Fraction(0)) for j in range(n)] for w in self.pivots]
        self.Minv = _invert(M)

    def coordinates(self, p):
        v = [p.get(w, Fraction(0)) for w in self.pivots]
        return {self.names[i]: c for i, row in enumerate(self.Minv)
                if self.names[i] is not None
                and (c := sum((a * b for a, b in zip(row, v) if a and b), Fraction(0)))}


def _invert(M):
    n = len(M)
    aug = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        pr = next(i for i in range(c, n) if aug[i][c] != 0)
        aug[c], aug[pr] = aug[pr], aug[c]
        inv = 1 / aug[c][c]
        aug[c] = [x * inv for x in aug[c]]
        for i in range(n):
            if i != c and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [a - f * b for a, b in zip(aug[i], aug[c])]
    return [row[n:] for row in aug]


@lru_cache(maxsize=None)
def _block(nA, nB):
    names = [n for n in names_of_degree(nA + nB) if multidegree(n) == (nA, nB)]
    return _Block(names)


def coordinates(p, check=False):
    """Coordinates of a homogeneous-or-not Lie polynomial p in the basis.

    p may contain several degrees/multidegrees; it must lie in the span of the
    basis (odd degrees 1..MAX_DEGREE).  With check=True the reconstruction is verified.
    """
    groups = {}
    for w, v in p.items():
        groups.setdefault((w.count('A'), w.count('B')), {})[w] = v
    out = {}
    for (nA, nB), sub in groups.items():
        if (nA + nB) % 2 == 0 or nA + nB > MAX_DEGREE:
            raise ValueError(f'multidegree {(nA, nB)} is not covered by the basis')
        if nA + nB == 1:
            out.update(sub)
            continue
        out.update(_block(nA, nB).coordinates(sub))
    if check:
        rec = {}
        for n, c in out.items():
            rec = fa.add(rec, poly(n), c)
        if fa.add(rec, p, -1):
            raise ValueError('polynomial is not in the span of the basis')
    return out


def _independent(polys, start=()):
    """Greedy selection: indices of polys independent of `start` and of each other."""
    echelon = []  # list of (pivot word, poly) in reduced form

    def reduce(p):
        p = dict(p)
        for w, b in echelon:
            if w in p:
                p = fa.add(p, b, -p[w] / b[w])
        return p

    for p in start:
        r = reduce(p)
        if r:
            echelon.append((max(r), r))
    chosen = []
    for i, p in enumerate(polys):
        r = reduce(p)
        if r:
            echelon.append((max(r), r))
            chosen.append(i)
    return chosen


def from_coordinates(coords):
    rec = {}
    for n, c in coords.items():
        rec = fa.add(rec, poly(n), c)
    return rec


# ---------------------------------------------------------------- A <-> B mirror symmetry
def mirror_name(name):
    """Partner of a basis element under A <-> B as used by NScheme::step_B:
    the coefficient arrays are reversed, i.e. K_k <-> K_{N+1-k}; A <-> B; C1 <-> C2."""
    if name in ('A', 'B'):
        return {'A': 'B', 'B': 'A'}[name]
    N = len(names_of_degree(degree(name)))
    return f'{name[0]}{N + 1 - int(name[1:])}'


def swap_letters(p):
    t = str.maketrans('AB', 'BA')
    return {w.translate(t): v for w, v in p.items()}


def mirror_sign(name):
    """+1 or -1 if swap_letters(K) = +-K_mirror, otherwise None (symmetry broken)."""
    p, q = swap_letters(poly(name)), poly(mirror_name(name))
    if p == q:
        return 1
    if p == fa.scale(q, -1):
        return -1
    return None
