"""Exploratory SU(4) four-nucleon 2^3 periodic-box correlator; not a physical He prediction.
Run from any directory. Writes results alongside this script. Published source:
arxiv:1911.06368 Eq.(lattice Hamiltonian), t/C0/D0 table at a=1.4 fm.
Synthesis counts use the existing report conventions, NOT compiled Clifford+T circuits.
"""
from pathlib import Path
import sys, json, math, itertools
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh, expm_multiply
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from estimates.common import t_per_rotation, eps_per_rotation
OUT=Path(__file__).resolve().parent
N=4096; ix=np.arange(N); pos=np.array([(ix >> (3*f)) & 7 for f in range(4)])
t=10.5794; c0=-98.2265511; d0=127.839693
V=np.zeros(N)
for a,b in itertools.combinations(range(4),2): V+=c0*(pos[a]==pos[b])
for a,b,c in itertools.combinations(range(4),3): V+=d0*((pos[a]==pos[b])&(pos[b]==pos[c]))
# Keep both directed periodic bonds: L=2 nearest neighbours coincide.
K=coo_matrix((np.full(12*N,-2*t),(np.tile(ix,12),np.concatenate([ix^(1<<b) for b in range(12)]))),shape=(N,N)).tocsr()
H=K+diags(V)
e,vec=eigsh(H,k=2,which='SA',tol=1e-12,v0=np.ones(N)); order=np.argsort(e); e=e[order]; psi=vec[:,order[0]]
if psi.sum()<0: psi=-psi
assert psi.min()>0
# Real positive state: 9 relative-coordinate bits, plus a uniform 3-bit origin.
r=pos[:3]^pos[3]; relative=r[0]+8*r[1]+64*r[2]
phi=psi[:512]*np.sqrt(8)  # fourth particle at origin
reconstruction=phi[relative]/np.sqrt(8)
load_error=float(np.linalg.norm(psi-reconstruction)); assert load_error<1e-10
# Explicit uniformly-controlled real-state loader angles, MSB first.
angles=[]; loaded=np.ones(1)
for level in range(9):
 weights=(phi**2).reshape(2**level,2,-1).sum(axis=2)
 theta=2*np.arctan2(np.sqrt(weights[:,1]),np.sqrt(weights[:,0])); angles.extend(theta)
 loaded=np.column_stack((loaded*np.cos(theta/2),loaded*np.sin(theta/2))).ravel()
assert np.linalg.norm(loaded-phi)<1e-10
np.save(OUT/'relative_state.npy',phi); np.save(OUT/'loader_multiplexor_angles.npy',np.array(angles))
# Walsh coefficients enumerate all diagonal Pauli strings in V, independently of pair counting.
def fwht(x):
 a=np.array(x,copy=True); h=1
 while h<len(a):
  z=a.reshape(-1,2*h); u=z[:,:h].copy(); v=z[:,h:].copy(); z[:,:h]=u+v; z[:,h:]=u-v; h*=2
 return a
coeff=fwht(V)/N; masks=np.flatnonzero(np.abs(coeff)>1e-10); nV=int(np.count_nonzero(masks)); identity=float(coeff[0])
assert np.max(np.abs(fwht(coeff)-V))<1e-10
# Proton species f=0,1. q=(pi/a,0,0): rho(q)=Z_0+Z_3; normalized J=rho/2.
z0=1-2*((ix>>0)&1); z1=1-2*((ix>>3)&1); J=(z0+z1)/2
ref=np.ones(N)/np.sqrt(N)
# Symmetric split e^-iKdt/2 e^-iVdt e^-iKdt/2, exact within each commuting layer.
def kinetic(v,dt):
 for b in range(12): v=np.cos(2*t*dt)*v+1j*np.sin(2*t*dt)*v[ix^(1<<b)]
 return v
def split(v,tau,steps):
 dt=tau/steps; v=kinetic(v.copy().astype(complex),dt/2)
 for step in range(steps):
  v*=np.exp(-1j*V*dt)
  v=kinetic(v,dt if step<steps-1 else dt/2)
 return v
rows=[]
for tau in [0.,.0025,.005,.01,.02]:
 target=J*psi; exact_state=expm_multiply((-1j*tau)*H,target)
 exact=np.vdot(target,exact_state)*np.exp(1j*e[0]*tau)
 # Same-H uncorrelated-state comparison is a transition amplitude, not a ground-state response.
 uncor=np.vdot(J*ref,expm_multiply((-1j*tau)*H,J*ref))*np.exp(1j*e[0]*tau)
 free=.5*np.exp(-1j*4*t*tau) # free ground-state response, own E0 removed
 for steps in [1,2,3,4,6,8]:
  approx=np.vdot(target,split(target,tau,steps))*np.exp(1j*e[0]*tau)
  # merged K boundaries: 12*(steps+1) rotations; V: steps*nV.
  # Every controlled Pauli rotation -> TWO elementary rotations and Clifford gates.
  # Identity of V -> one ancilla phase rotation per V layer; +1 for E0 phase.
  nr=511+2*(12*(steps+1)+nV*steps)+steps+1
  eps=eps_per_rotation(.005,nr,'coherent')
  tr=nr*t_per_rotation(eps,'rs')
  rus=nr*t_per_rotation(eps_per_rotation(.01,nr),'rus')
  rows.append(dict(tau_MeV_inv=tau,time_fm_c=tau*197.3269804,steps=steps,exact=[exact.real,exact.imag],approx=[approx.real,approx.imag],abs_trotter_error=abs(approx-exact),free=[free.real,free.imag],uncorrelated_same_H=[uncor.real,uncor.imag],rotations=nr,T_RS_leading_estimate=tr,T_report_RUS_expected=rus))
# Dense/sparse physics and estimator checks.
residual=float(np.linalg.norm(H@psi-e[0]*psi)); assert residual<1e-8
# Four Hadamard circuits directly measure Z_a U Z_b; uniform random pair yields C/4.
tau=.01; terms=[]
for za in [z0,z1]:
 for zb in [z0,z1]: terms.append(np.vdot(za*psi,expm_multiply(-1j*tau*H,zb*psi)))
assert abs(sum(terms)/4-np.vdot(J*psi,expm_multiply(-1j*tau*H,J*psi)))<1e-10
result=dict(model='SU(4), one particle per spin-isospin species, 2x2x2 periodic box',a_fm=1.4,q_MeV=math.pi*197.3269804/1.4,t_MeV=t,C0_MeV=c0,D0_MeV=d0,energy_without_kinetic_constant_MeV=float(e[0]),energy_with_kinetic_constant_MeV=float(e[0]+24*t),gap_MeV=float(e[1]-e[0]),residual=residual,loader_error=load_error,loader_rotations=511,logical_qubits_deterministic=13,V_nonidentity_Z_strings=nV,V_identity_MeV=identity,static_normalized_strength=float(np.vdot(J*psi,J*psi)),free_static_strength=.5,rows=rows)
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
for row in rows:
 if row['steps'] in [2,4,8] and row['tau_MeV_inv'] in [.005,.01,.02]: print(row)
