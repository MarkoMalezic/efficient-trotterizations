"""Recursive formulae for the basis coefficients of a symmetric Trotter-Suzuki scheme.

Replaces the notebook chapter "Recursive formula derivation".

The scheme exponent after i stages is

    Phi = nu A + sigma B + alpha C1 + beta C2 + sum_k gamma_k D_k + sum_k delta_k E_k
          + sum_k epsilon_k F_k                                          (h-powers implicit)

and a stage with the A (or B) operator updates it as

    step A:  Phi' = log( exp(a A/2) exp(Phi) exp(a A/2) )
    step B:  Phi' = log( exp(b B/2) exp(Phi) exp(b B/2) )

Writing the symmetric BCH formula as

    log(exp(X/2) exp(Y) exp(X/2)) = X + Y + sum_K Q_K K(X, Y)

(Q_K = alpha_Q, beta_Q, gammas_Q[k], ... from prefactors.py; K(X, Y) = basis element with
A -> X, B -> Y) we substitute X = a A, Y = Phi, expand every K multilinearly over the
components of Phi (the notebook's term-by-term expansion, e.g. C1[a A, Phi1 + Phi3 + ...]),
and express each resulting nested commutator in the basis (the notebook's Hall word
transformations).

Each term of the result comes from exactly one K, so the coefficient of every
(target, monomial) is a linear combination of the prefactors Q_K with rational
coefficients.  We keep it in that form, exactly like the notebook and NCoefficients.cpp,
e.g.  gamma1 += (gammas_Q[1] + gammas_Q[2]) a^3 nu sigma.

Usage:
    python -m pybch.recursion --order 4             # increments of all coefficients up to order 4, step A
    python -m pybch.recursion --order 8 --step B    # step B, up to order 8
    python -m pybch.recursion --order 6 --numeric   # prefactors Q_K replaced by their values
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction

import sympy as sp

from . import freealg as fa
from .basis import ALL_NAMES, coordinates, degree, names_of_degree, poly, word_of
from .bch import symmetric_bch_coefficients

# ---------------------------------------------------------------- naming
# Scheme coefficients and prefactors follow NCoefficients.cpp (0-based arrays).
LEVEL_ARRAY = {'D': 'gammas', 'E': 'deltas', 'F': 'epsilons'}
STEP_VAR = {'A': 'a', 'B': 'b'}


def var_name(basis_name):
    """Name of the scheme coefficient multiplying a basis element in Phi, e.g. 'gammas[0]'."""
    special = {'A': 'nu', 'B': 'sigma', 'C1': 'alpha', 'C2': 'beta'}
    if basis_name in special:
        return special[basis_name]
    return f'{LEVEL_ARRAY[basis_name[0]]}[{int(basis_name[1:]) - 1}]'


def prefactor_name(basis_name):
    """Name of the BCH prefactor Q_K, e.g. 'gammas_Q[0]'."""
    if basis_name in ('C1', 'C2'):
        return {'C1': 'alpha_Q', 'C2': 'beta_Q'}[basis_name]
    return f'{LEVEL_ARRAY[basis_name[0]]}_Q[{int(basis_name[1:]) - 1}]'


def max_degree(order):
    return order + 1


# ---------------------------------------------------------------- derivation
class _Expander:
    def __init__(self):
        self._nest = {}
        self._coords = {}

    def nest(self, names):
        """Right-nested commutator of basis elements, memoised on suffixes."""
        if names in self._nest:
            return self._nest[names]
        if len(names) == 1:
            r = poly(names[0])
        else:
            inner = self.nest(names[1:])
            r = fa.comm(poly(names[0]), inner) if inner else {}
        self._nest[names] = r
        return r

    def coords(self, names):
        if names not in self._coords:
            p = self.nest(names)
            self._coords[names] = coordinates(p) if p else {}
        return self._coords[names]


def recursive_formulas(step='A', order=8):
    """Increments of all coefficients up to `order` for one A (or B) stage.

    Returns {target basis name: {monomial: {K: Fraction}}} where
      * monomial is a sorted tuple of (variable, power) in the step parameter ('a' or 'b')
        and the old coefficients ('nu', 'sigma', 'alpha', 'gammas[0]', ...);
      * {K: c} is the linear form  sum_K c * Q_K  in the prefactors.
    The trivial parts (old value of the coefficient itself, and nu += a / sigma += b)
    are not included: these are the '+=' terms of NCoefficients.
    """
    maxdeg = max_degree(order)
    comps = [n for n in ALL_NAMES if degree(n) <= maxdeg]
    by_deg = defaultdict(list)
    for n in comps:
        by_deg[degree(n)].append(n)
    x = STEP_VAR[step]
    ex = _Expander()
    res = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))

    for d in range(3, maxdeg + 1, 2):
        for K in names_of_degree(d):
            sign, word = word_of(K)
            slots = ['X' if ch == 'A' else 'Y' for ch in word]
            nX = slots.count('X')
            budget = maxdeg - nX

            def assign(i, chosen, used):
                if i == len(slots):
                    yield tuple(chosen)
                    return
                if slots[i] == 'X':
                    chosen.append(None)
                    yield from assign(i + 1, chosen, used)
                    chosen.pop()
                    return
                rest = sum(1 for s in slots[i + 1:] if s == 'Y')
                for dd, names in by_deg.items():
                    if used + dd + rest > budget:
                        continue
                    for n in names:
                        chosen.append(n)
                        yield from assign(i + 1, chosen, used + dd)
                        chosen.pop()

            for choice in assign(0, [], 0):
                names = tuple(step if c is None else c for c in choice)
                cs = ex.coords(names)
                if not cs:
                    continue
                mono = Counter({x: nX})
                for c in choice:
                    if c is not None:
                        mono[var_name(c)] += 1
                mono = tuple(sorted(mono.items()))
                for target, val in cs.items():
                    res[target][mono][K] += sign * val

    out = {}
    for t, terms in res.items():
        clean = {}
        for m, lin in terms.items():
            lin = {K: c for K, c in lin.items() if c}
            if lin:
                clean[m] = lin
        if clean:
            out[t] = clean
    return out


# ---------------------------------------------------------------- evaluation / sympy
def numeric(formulas, Q=None):
    """Replace the prefactor linear forms by numbers: {target: {monomial: Fraction}}."""
    if Q is None:
        Q = symmetric_bch_coefficients()
    out = {}
    for t, terms in formulas.items():
        num = {}
        for m, lin in terms.items():
            v = sum((c * Q[K] for K, c in lin.items()), Fraction(0))
            if v:
                num[m] = v
        out[t] = num
    return out


def apply_step(values, formulas, step_value, step='A', Q=None):
    """One stage evaluated exactly.  values: {basis name: Fraction}; returns new values."""
    num = numeric(formulas, Q)
    env = {var_name(n): v for n, v in values.items()}
    env[STEP_VAR[step]] = Fraction(step_value)
    out = dict(values)
    for target, terms in num.items():
        s = Fraction(0)
        for mono, c in terms.items():
            t = c
            for v, p in mono:
                t *= env.get(v, 0) ** p
                if not t:
                    break
            s += t
        out[target] = out.get(target, 0) + s
    out[step] = out.get(step, 0) + Fraction(step_value)
    return {k: v for k, v in out.items() if v}


_SYMS = {}


def symbol(name):
    """sympy object for 'nu', 'a', 'alpha_Q', 'gammas[0]', 'gammas_Q[3]', ..."""
    if name not in _SYMS:
        if '[' in name:
            base, idx = name[:-1].split('[')
            _SYMS[name] = sp.IndexedBase(base)[int(idx)]
        else:
            _SYMS[name] = sp.Symbol(name)
    return _SYMS[name]


def to_sympy(formulas, Q=None):
    """{target var name (e.g. 'gammas[0]'): sympy expression of the increment}.

    The prefactors stay symbolic (alpha_Q, gammas_Q[k], ...) unless their values
    Q = {basis name: Fraction} are given, e.g. Q = symmetric_bch_coefficients().
    """
    def prefactor(K):
        if Q is None:
            return symbol(prefactor_name(K))
        return sp.Rational(Q[K].numerator, Q[K].denominator)

    out = {}
    for t, terms in formulas.items():
        expr = sp.S(0)
        for m, lin in terms.items():
            q = sum((sp.Rational(c.numerator, c.denominator) * prefactor(K)
                     for K, c in lin.items()), sp.S(0))
            mon = sp.Mul(*[symbol(v) ** p for v, p in m])
            expr += q * mon
        out[var_name(t)] = expr
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--order', type=int, default=4, choices=[2, 4, 6, 8])
    ap.add_argument('--step', default='A', choices=['A', 'B'])
    ap.add_argument('--numeric', action='store_true', help='insert the values of the prefactors')
    args = ap.parse_args()
    F = recursive_formulas(args.step, args.order)
    exprs = to_sympy(F, symmetric_bch_coefficients() if args.numeric else None)
    x, v = (('a', 'nu'), ('b', 'sigma'))[args.step == 'B']
    print(f'# step {args.step} up to order {args.order}: every coefficient is incremented by the expression'
          f' below (old values on the right-hand side); in addition {v} += {x}')
    for n in ALL_NAMES:
        t = var_name(n)
        if t in exprs:
            print(f'{t} += {exprs[t]}')


if __name__ == '__main__':
    main()
