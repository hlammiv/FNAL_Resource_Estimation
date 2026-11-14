"""Explicit JW-algebra check of hamiltonian_strings_count: a generic rank-0 one- plus two-body Hamiltonian with
random coefficients (one-body hops between equal l, j, m, the oscillator case), strings collected exactly.
Usage (repository root): python3 applications/ch03_mu2e_0nubb/pauli/explicit_check.py sd emax1 emax2 mo_ru
Output 2026-11-13: sd 12900, emax1 1800, emax2 57500, mo_ru 114909, each equal to the structured count."""
import itertools, pathlib, random, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from estimates import ch03_mu2e_0nubb as c
def explicit(orb_p, orb_n, seed=1):
    rnd = random.Random(seed)
    orbits = tuple(orb_p) + tuple(orb_n)
    S = [(o, tj, m2, "p") for o, (_, _, tj) in enumerate(orb_p) for m2 in range(-tj, tj+1, 2)] + \
        [(o+len(orb_p), tj, m2, "n") for o, (_, _, tj) in enumerate(orb_n) for m2 in range(-tj, tj+1, 2)]
    st = [s[:3] for s in S]
    pairs = list(itertools.combinations(range(len(S)), 2))
    info = {pr: c.pair_quantum_numbers(orbits, st[pr[0]], st[pr[1]]) + ("".join(sorted(S[pr[0]][3]+S[pr[1]][3])),) for pr in pairs}
    acc = {}
    def add(terms, cf):
        for cc, x, z in terms:
            acc[(x, z)] = acc.get((x, z), 0.0) + cf*cc
    for ia in range(len(pairs)):
        for ib in range(ia, len(pairs)):
            A, B = info[pairs[ia]], info[pairs[ib]]
            if A[3] != B[3] or not (A[0]==B[0] and A[1]==B[1] and (A[2] & B[2])): continue
            (i,j),(k,l) = pairs[ia], pairs[ib]
            cf = rnd.uniform(0.5,1.5)
            add(c.jw_quartic(i,j,l,k), cf)
            if ia != ib: add(c.jw_quartic(k,l,j,i), cf)
    for p in range(len(S)):
        for q in range(p, len(S)):
            op, oq = orbits[S[p][0]], orbits[S[q][0]]
            if S[p][3]==S[q][3] and op[1]==oq[1] and op[2]==oq[2] and S[p][2]==S[q][2]:
                cf = rnd.uniform(0.5,1.5)
                add(c._pmul(c._ladder(p,True), c._ladder(q,False)), cf)
                if p != q: add(c._pmul(c._ladder(q,True), c._ladder(p,False)), cf)
    nz = {k for k, v in acc.items() if abs(v) > 1e-12}
    return len(nz - {(0,0)})
for k in sys.argv[1:]:
    sp = c.step_space(k)
    print(k, explicit(*sp), c.hamiltonian_strings_count(*sp)["gamma"], flush=True)
