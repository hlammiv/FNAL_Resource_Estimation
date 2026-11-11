"""As run_gaps.py, with a J = J0 projected start vector (Lowdin projector over J = J0+1 .. Jmax), so the Krylov
space starts in the J0 sector; the J^2 shift keeps round-off leakage from J != J0 states out of the low end."""
from sm import *
import json, sys
orbs = {'f7p3': {(0, 3, 7), (1, 1, 3)}, 'f7p3p1': {(0, 3, 7), (1, 1, 3), (1, 1, 1)},
        'pf': {(0, 3, 7), (1, 1, 3), (1, 1, 1), (0, 3, 5)}, 'sd': {(0, 2, 5), (1, 0, 1), (0, 2, 3)}}
nucs = {'48Ca': (0, 8, 48), '48Ti': (2, 6, 48), '27Al': (5, 6, 27)}
snt, space, nuc, M2, J0, k = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
Z, N, A = nucs[nuc]
t0 = time.time()
m = Model(snt + '.snt', orbs[space], Z, N, A)
m.penalty = (1.0, J0)
m.setup_sector(M2)
rng = np.random.default_rng(1)
Psi = rng.standard_normal(m.mask.shape) * m.mask
jmax2 = sum(sorted(m.m2[-1], reverse=True)[:Z]) + sum(sorted(m.m2[1], reverse=True)[:N])
J0x2 = int(round(2 * J0))
for j2 in range(J0x2 + 2, int(jmax2) + 1, 2):
    j = j2 / 2
    Psi = m.J2_full(Psi) - j * (j + 1) * Psi
    Psi /= np.linalg.norm(Psi)
v0 = Psi[m.mask]
print('jmax', jmax2 / 2, 'start <J2>', m.J2(v0 / np.linalg.norm(v0)), 't', round(time.time() - t0, 1), flush=True)
dim = int(m.mask.sum())
count = [0]


def mv(x):
    count[0] += 1
    P = np.zeros(m.mask.shape)
    P[m.mask] = x
    return m.matvec_full(P)[m.mask]


op = LinearOperator((dim, dim), matvec=mv, dtype=float)
w, v = eigsh(op, k=k, which='SA', tol=1e-10, ncv=20, v0=v0)
m.penalty = None
E, Js = [], []
for i in range(v.shape[1]):
    P = np.zeros(m.mask.shape)
    P[m.mask] = v[:, i]
    E.append(float(np.sum(P * m.matvec_full(P))))
    j2 = m.J2(v[:, i])
    Js.append((-1 + math.sqrt(1 + 4 * j2)) / 2)
o = np.argsort(E)
E = [E[i] for i in o]
Js = [Js[i] for i in o]
print(snt, space, nuc, 'M2', M2, 'dim', dim, 'J0', J0, 'E', [round(e, 5) for e in E], 'Ex', [round(e - E[0], 4) for e in E],
      'J', [round(j, 3) for j in Js], 'matvecs', count[0], 't', round(time.time() - t0, 1), flush=True)
json.dump(dict(snt=snt, space=space, nuc=nuc, M2=M2, dim=dim, J0=J0, E=E, J=Js, method='projected start',
               matvecs=count[0]), open(f'gapP_{snt}_{space}_{nuc}_J{J0}.json', 'w'), indent=1)
