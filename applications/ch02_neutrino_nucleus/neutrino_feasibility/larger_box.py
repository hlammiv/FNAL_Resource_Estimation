"""3^3 four-species SU(4) study in the zero-momentum relative-coordinate sector.
Exploratory analysis only; manuscript inputs are untouched.
"""
from pathlib import Path
import json, itertools, time
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import eigsh,expm_multiply
from scipy.optimize import minimize
P=Path(__file__).resolve().parent
L=3; M=27; N=M**3; t=10.5794; C0=-98.2265511; D0=127.839693
ix=np.arange(N); sites=np.array([ix%M,(ix//M)%M,ix//M**2,np.zeros(N,dtype=int)])
xyz=np.array([sites%L,(sites//L)%L,sites//L**2]).transpose(1,0,2)
features=[]; labels=[]
for d in range(3):
 features.append(sum((np.count_nonzero(xyz[a]!=xyz[b],axis=0)==d).astype(float) for a,b in itertools.combinations(range(4),2))); labels.append('pairs_distance'+str(d))
features.append(sum(((sites[a]==sites[b])&(sites[b]==sites[c])).astype(float) for a,b,c in itertools.combinations(range(4),3)));labels.append('triples_onsite')
features.append(np.all(sites==sites[0],axis=0).astype(float));labels.append('quadruples_onsite')
F=np.array(features).T
V=C0*F[:,0]+D0*F[:,3]
# Map absolute single-particle shifts into r_i=x_i-x_3. Each outgoing bond has amplitude -t.
perms=[]; particle_shifts=[]
for f in range(4):
 for axis in range(3):
  for direction in [-1,1]:
   new=xyz[:3].copy()
   if f<3: new[f,axis]=(new[f,axis]+direction)%3
   else: new[:,axis]=(new[:,axis]-direction)%3
   ns=new[:,0]+3*new[:,1]+9*new[:,2]
   perms.append(ns[0]+27*ns[1]+729*ns[2]); particle_shifts.append((f,axis,direction))
perms=np.array(perms)
K=coo_matrix((np.full(24*N,-t),(np.tile(ix,24),perms.ravel())),shape=(N,N)).tocsr(); H=K+diags(V)
print('diagonalizing relative sector',N,flush=True)
e,vs=eigsh(H,k=2,which='SA',tol=1e-11,v0=np.ones(N)); ii=np.argsort(e);e=e[ii];psi=vs[:,ii[0]];psi*=np.sign(psi.sum())
assert psi.min()>0
print('E',e[0]+24*t,flush=True)
def objective(theta):
 log=F@theta; state=np.exp(log-log.max());state/=np.linalg.norm(state)
 hs=H@state; energy=state@hs
 grad=2*((F*state[:,None]).T@(hs-energy*state))
 return energy,grad
opt=minimize(objective,[.3,.1,.03,-.1,0],jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':300})
log=F@opt.x; trial=np.exp(log-log.max());trial/=np.linalg.norm(trial)
# Filter each pair/triple/quad independently, with its own maximal eigenvalue shifted to zero.
# pair distance=3 coefficient is zero. triples and quad filters have two values {0, coefficient}.
shift=6*max(0.,*opt.x[:3])+4*max(0.,opt.x[3])+max(0.,opt.x[4])
success=float(np.mean(np.exp(2*(log-shift))))
print('opt',opt.fun+24*t,'fidelity',float((psi@trial)**2),'success',success,'theta',opt.x,flush=True)
# J(q)=(e^iq x_p0+e^iq x_p1)/2. Relative state carries total momentum q after insertion.
q=2*np.pi/3; current=(np.exp(1j*q*xyz[0,0])+np.exp(1j*q*xyz[1,0]))/2
# Wavefunction after current is exp(iq x_3)*g(r). Fourth-particle shifts gain exp(iq direction).
phases=np.array([np.exp(1j*q*d) if f==3 and a==0 else 1 for f,a,d in particle_shifts])
Kq=coo_matrix((np.repeat(-t*phases,N),(np.tile(ix,24),perms.ravel())),shape=(N,N)).tocsr(); Hq=Kq+diags(V)
assert np.max(np.abs((Hq-Hq.conj().T).data),initial=0)<1e-12
rows=[]
# K and V split with K applied exactly by sparse expm; circuit counts handled separately.
def split(v,tau,n):
 dt=tau/n; v=expm_multiply(-1j*dt/2*Kq,v)
 for j in range(n):
  v*=np.exp(-1j*dt*V)
  v=expm_multiply(-1j*dt*(1 if j<n-1 else .5)*Kq,v)
 return v
for tau in [0,.0025,.005,.01,.02]:
 target=current*psi; approximate=current*trial
 exact=np.vdot(target,expm_multiply(-1j*tau*Hq,target))*np.exp(1j*e[0]*tau)
 prepared=np.vdot(approximate,expm_multiply(-1j*tau*Hq,approximate))*np.exp(1j*e[0]*tau)
 for steps in [1,2,4]:
  sim=np.vdot(approximate,split(approximate,tau,steps))*np.exp(1j*e[0]*tau)
  rows.append(dict(tau=tau,steps=steps,exact=[exact.real,exact.imag],trial_exact=[prepared.real,prepared.imag],circuit=[sim.real,sim.imag],state_bias=abs(prepared-exact),trotter_error=abs(sim-prepared),total_bias=abs(sim-exact)))
print('finished responses',flush=True)
out=dict(lattice='3x3x3',dimension=N,parameters=dict(t=t,C0=C0,D0=D0,a_fm=1.4),ground_energy_MeV=float(e[0]+24*t),ground_residual=float(np.linalg.norm(H@psi-e[0]*psi)),variational_energy_MeV=float(opt.fun+24*t),optimizer_success=bool(opt.success),optimizer_message=str(opt.message),gradient_norm=float(np.linalg.norm(opt.jac)),features=labels,theta=opt.x.tolist(),state_fidelity=float((psi@trial)**2),filter_success_probability=success,exact_static_strength=float(np.vdot(current*psi,current*psi).real),trial_static_strength=float(np.vdot(current*trial,current*trial).real),rows=rows)
(P/'larger_box_results.json').write_text(json.dumps(out,indent=2)+'\n');np.savez(P/'larger_box_states.npz',exact=psi,trial=trial,theta=opt.x)
print(json.dumps(out,indent=2),flush=True)
