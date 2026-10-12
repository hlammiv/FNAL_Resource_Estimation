"""Exploratory two-channel contact benchmark; does not edit the manuscript.

Exact relative-coordinate spectrum, checked against full two-position dynamics.
No physical cross section or quantum advantage is asserted.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.linalg import eigh

OUT = Path(__file__).resolve().parent
t, a, binding, scattering = 10.5794, 1.4, 2.224575, -23.76

def contacts(n):
    x, w = roots_legendre(n)
    k = (x + 1) * np.pi / 2
    w = w / 2
    ek = 4*t*(3-np.cos(k[:,None,None])-np.cos(k[None,:,None])-np.cos(k[None,None,:]))
    weights = w[:,None,None]*w[None,:,None]*w[None,None,:]
    gt = np.sum(weights/(-binding-ek))
    g0 = np.sum(weights/(-ek))
    return float(1/gt), float(1/(g0+a/(8*np.pi*t*scattering)))

ct, cs = contacts(160)
output = dict(parameters=dict(t=t, a_fm=a, triplet_C_MeV=ct, singlet_C_MeV=cs,
                              infinite_volume_binding_MeV=binding, singlet_a_fm=scattering),
              quadrature_checks={str(n): contacts(n) for n in (60,100,160)}, volumes=[])

for L in (4,6):
    m = L**3
    xyz = np.array(np.unravel_index(np.arange(m),(L,L,L))).T
    K = np.eye(m)*(12*t)
    for axis in range(3):
        for sign in (-1,1):
            other = xyz.copy(); other[:,axis] = (other[:,axis]+sign)%L
            ix = np.ravel_multi_index(other.T,(L,L,L))
            K[np.arange(m),ix] -= 2*t
    ht = K.copy(); ht[0,0] += ct
    hs = K.copy(); hs[0,0] += cs
    et, vt = eigh(ht); es, vs = eigh(hs)
    psi = vt[:,0]; psi *= np.sign(psi.sum())
    weights = (vs.T@psi)**2
    entry = dict(L=L, initial_energy=float(et[0]), final_lowest_energy=float(es[0]),
                 transitions=[dict(energy=float(e-et[0]), strength=float(w))
                              for e,w in zip(es,weights) if w>1e-6])
    if L == 6:
        # Translation invariant wavefunction psi(r1-r2)/sqrt(number of sites).
        diff=(xyz[:,None,:]-xyz[None,:,:])%L
        rel=np.ravel_multi_index(diff.reshape(-1,3).T,(L,L,L)).reshape(m,m)
        source=psi[rel]/np.sqrt(m)
        permutations=[]
        for axis in range(3):
            for parity in (0,1):
                dest=xyz.copy()
                dest[:,axis]=(dest[:,axis]+np.where(xyz[:,axis]%2==parity,1,-1))%L
                permutations.append(np.ravel_multi_index(dest.T,(L,L,L)))
        def match(v, p, dt):
            c,s=np.cos(t*dt),1j*np.sin(t*dt)
            v=c*v+s*v[p,:]
            return c*v+s*v[:,p]
        def step(v, dt):
            for p in permutations: v=match(v,p,dt/2)
            v[np.diag_indices(m)]*=np.exp(-1j*cs*dt)
            for p in permutations[::-1]: v=match(v,p,dt/2)
            return v*np.exp(-1j*(12*t-et[0])*dt)
        checks=[]
        for dt in (.005,.0025,.00125):
            v=source.astype(complex)
            for j in range(1,round(.5/dt)+1):
                v=step(v,dt)
                tau=j*dt
                if any(abs(tau-x)<1e-9 for x in (.05,.1,.25,.5)):
                    exact=np.sum(weights*np.exp(-1j*(es-et[0])*tau))
                    approx=np.vdot(source,v)
                    checks.append(dict(tau=tau, dt=dt, steps=j, exact=[exact.real,exact.imag],
                                       trotter=[approx.real,approx.imag], bias=float(abs(exact-approx)),
                                       norm=float(np.linalg.norm(v))))
        entry['correlator_checks']=checks
        np.savez(OUT/'two_body_states.npz', initial_relative=psi, final_energies=es,
                 strengths=weights, initial_full=source)
    output['volumes'].append(entry)

(OUT/'two_body_results.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
