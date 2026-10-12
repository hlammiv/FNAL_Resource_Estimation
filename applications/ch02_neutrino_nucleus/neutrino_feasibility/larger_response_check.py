"""Observable check of deterministic preparation; also verifies momentum reduction.
Reuses the explicitly assembled small sparse Hamiltonians from larger_box.py.
"""
from pathlib import Path
import json,sys
import numpy as np
from scipy.sparse.linalg import expm_multiply
import larger_box as b
P=Path(__file__).resolve().parent
trial=np.load(P/'unitary_state_8.npy'); exact=b.psi; E=b.e[0]
rows=[]
for tau,steps in [(0,2),(.005,2),(.01,2),(.02,4),(.03,6)]:
 target=b.current*exact; prep=b.current*trial
 truth=np.vdot(target,expm_multiply(-1j*tau*b.Hq,target))*np.exp(1j*E*tau)
 prepared=np.vdot(prep,expm_multiply(-1j*tau*b.Hq,prep))*np.exp(1j*E*tau)
 approx=np.vdot(prep,b.split(prep,tau,steps))*np.exp(1j*E*tau)
 # More Trotter steps at the final point validate convergence separately from prep bias.
 twice=np.vdot(prep,b.split(prep,tau,2*steps))*np.exp(1j*E*tau)
 rows.append(dict(tau_MeV_inv=tau,time_fm_c=tau*197.3269804,steps=steps,exact=[truth.real,truth.imag],prepared_exact=[prepared.real,prepared.imag],circuit=[approx.real,approx.imag],state_bias=abs(prepared-truth),trotter_error=abs(approx-prepared),trotter_error_double_steps=abs(twice-prepared),total_bias=abs(approx-truth),free=[.5*np.cos(3*b.t*tau),-.5*np.sin(3*b.t*tau)]))
# Momentum reduction verification by applying full 27^4 Hamiltonian to a lifted arbitrary vector.
# This does not simulate all 36 padded/workspace qubits. It checks the physical legal sector.
rng=np.random.default_rng(42);v=rng.normal(size=b.N)+1j*rng.normal(size=b.N); v/=np.linalg.norm(v)
# Check 100 full-coordinate rows independently against the reduced Hq action.
red=b.Hq@v; maxerr=0.
def value(pos):
 rel=(pos[:3]-pos[3])%3;sites=rel[:,0]+3*rel[:,1]+9*rel[:,2]
 return np.exp(1j*b.q*pos[3,0])*v[sites[0]+27*sites[1]+729*sites[2]]/np.sqrt(27)
for _ in range(100):
 pos=rng.integers(0,3,size=(4,3)); val=0j
 for f in range(4):
  for axis in range(3):
   for direction in [-1,1]:
    new=pos.copy();new[f,axis]=(new[f,axis]+direction)%3;val-=b.t*value(new)
 pairs=sum(np.array_equal(pos[a],pos[c]) for a,c in __import__('itertools').combinations(range(4),2))
 triples=sum(np.array_equal(pos[a],pos[c]) and np.array_equal(pos[c],pos[d]) for a,c,d in __import__('itertools').combinations(range(4),3))
 val+=(b.C0*pairs+b.D0*triples)*value(pos)
 rel=(pos[:3]-pos[3])%3;sites=rel[:,0]+3*rel[:,1]+9*rel[:,2];j=sites[0]+27*sites[1]+729*sites[2]
 check=np.exp(1j*b.q*pos[3,0])*red[j]/np.sqrt(27);maxerr=max(maxerr,abs(val-check))
assert maxerr<1e-10
out=dict(rows=rows,momentum_sector_max_error=maxerr)
(P/'larger_unitary_response.json').write_text(json.dumps(out,indent=2)+'\n');print('UNITARY RESPONSE',json.dumps(out,indent=2))
