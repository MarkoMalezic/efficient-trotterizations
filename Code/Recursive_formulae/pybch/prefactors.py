"""Basis prefactors of the symmetric BCH formula (the `Prefactors` of the C++ code).

    log(exp(A/2) exp(B) exp(A/2)) = A + B + alpha_Q C1 + beta_Q C2 + sum_k gammas_Q[k] D_k
                                    + sum_k deltas_Q[k] E_k + sum_k epsilons_Q[k] F_k
                                    + O(h^11)

Notebook names: aQ, bQ, gQ, dQ, eQ.  Order n of the scheme <-> degree n+1 of the
basis elements (order 2: C, 4: D, 6: E, 8: F).

Usage:
    python -m pybch.prefactors                 # print all orders
    python -m pybch.prefactors --order 8       # only order 8
    python -m pybch.prefactors --scheme full   # log(exp(A) exp(B) exp(A)) instead
    python -m pybch.prefactors --write         # also write pybch/data/prefactors_<scheme>.{json,m}
"""
import argparse
import json
import os
import time

import sympy as sp

from .basis import names_of_degree
from .bch import symmetric_bch_coefficients

ORDER_TO_LEVEL = {2: 'C', 4: 'D', 6: 'E', 8: 'F'}
CPP_NAMES = {2: ('alpha', 'beta'), 4: 'gammas', 6: 'deltas', 8: 'epsilons'}
NB_NAMES = {4: 'gQ', 6: 'dQ', 8: 'eQ'}


def prefactors(max_order=8, scheme='half'):
    """{order: [sympy.Rational, ...]} in basis order (C1, C2 | D1..D6 | E1..E18 | F1..F56)."""
    Q = symmetric_bch_coefficients(max_order + 1, scheme)
    return {n: [sp.Rational(Q[k].numerator, Q[k].denominator) for k in names_of_degree(n + 1)]
            for n in range(2, max_order + 1, 2)}


def to_mathematica(pref):
    lines = []
    for n, vals in pref.items():
        if n == 2:
            lines.append(f'aQ={sp.mathematica_code(vals[0])};')
            lines.append(f'bQ={sp.mathematica_code(vals[1])};')
        else:
            lines.append(f'{NB_NAMES[n]}={{' + ','.join(sp.mathematica_code(v) for v in vals) + '};')
    return '\n'.join(lines) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--order', type=int, choices=[2, 4, 6, 8])
    ap.add_argument('--scheme', default='half', choices=['half', 'full'])
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    t = time.time()
    pref = prefactors(args.order or 8, args.scheme)
    if args.order:
        pref = {args.order: pref[args.order]}
    for n, vals in pref.items():
        level = ORDER_TO_LEVEL[n]
        print(f'# order {n}: {len(vals)} prefactors ({level}1..{level}{len(vals)})')
        for k, v in enumerate(vals, 1):
            print(f'  {level}{k:<4} {v}')
    print(f'# computed in {time.time() - t:.1f}s')
    if args.write:
        d = os.path.join(os.path.dirname(__file__), 'data')
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f'prefactors_{args.scheme}.json'), 'w') as f:
            json.dump({str(n): [str(v) for v in vals] for n, vals in pref.items()}, f, indent=1)
        with open(os.path.join(d, f'prefactors_{args.scheme}.m'), 'w') as f:
            f.write(to_mathematica(pref))
        print(f'# written to {d}')


if __name__ == '__main__':
    main()
