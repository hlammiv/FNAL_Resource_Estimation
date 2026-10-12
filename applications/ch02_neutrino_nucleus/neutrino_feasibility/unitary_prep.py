"""Optimize deterministic kinetic/contact pulse preparation on the 3^3 model.
Explicit classical variational experiment, no manuscript mutation.
"""
from pathlib import Path
import json,itertools,time
import numpy as np
from scipy.optimize import minimize
P=Path(__file__).resolve().parent; data=np.load(P/'larger_box_states.npz'); exact=data['exact'];N=27**3
idx=np.indices((3,)*9); # axis order grouped as particle2 (z,y,x), particle1, particle0
p=[idx[6:9],idx[3:6],idx[0:3],np.zeros_like(idx[0:3])]
f2=sum(np.all(p[a]==p[b],axis=0).astype(float) for a,b in itertools.combinations(range(4),2))
f3=sum((np.all(p[a]==p[b],axis=0)&np.all(p[b]==p[c],axis=0)).astype(float) for a,b,c in itertools.combinations(range(4),3))
V=(-98.2265511*f2+127.839693*f3)/100
K=-2*10.5794*(sum(np.cos(2*np.pi*x/3).sum(axis=0) for x in p[:3])+np.cos(2*np.pi*(p[0]+p[1]+p[2])/3).sum(axis=0))/100
psi0=np.ones((3,)*9,dtype=complex)/np.sqrt(N)
def fft(v):return np.fft.fftn(v,norm='ortho')
def ifft(v):return np.fft.ifftn(v,norm='ortho')
def apply(v,angle,kind):return ifft(np.exp(-1j*angle*K)*fft(v)) if kind else np.exp(-1j*angle*V)*v
def generator(v,kind):return ifft(K*fft(v)) if kind else V*v
def state(theta):
 s=psi0
 for i,a in enumerate(theta):s=apply(s,a,i%2)
 return s
def objective(theta):
 states=[psi0]
 for i,a in enumerate(theta):states.append(apply(states[-1],a,i%2))
 s=states[-1]; back=generator(s,1)+generator(s,0);E=np.vdot(s,back).real
 grad=np.empty(len(theta))
 for i in reversed(range(len(theta))):
  grad[i]=2*np.vdot(back,-1j*generator(states[i+1],i%2)).real
  back=apply(back,-theta[i],i%2)
 return E,grad
# Explicitly check FFT kinetic against the known ground energy (ordering/mapping check).
h=generator(exact.reshape((3,)*9),0)+generator(exact.reshape((3,)*9),1)
E=float(np.vdot(exact,h.ravel()).real*100)
assert np.linalg.norm(h.ravel()-E/100*exact)<1e-9
# Finite-difference independent gradient check.
th=np.array([.12,.2,.21,.31]);val,g=objective(th)
fd=np.array([(objective(th+np.eye(4)[k]*1e-6)[0]-objective(th-np.eye(4)[k]*1e-6)[0])/2e-6 for k in range(4)])
assert np.max(abs(fd-g))<1e-7
rng=np.random.default_rng(20261003); results=[];previous=None
for layers in [2,4,6,8]:
 seeds=[]
 if previous is not None:
  seeds.append(np.r_[previous,[.02,.02]*(layers-len(previous)//2)])
 for total in [2.,5.]:
  seed=np.empty(2*layers);seed[0::2]=total/layers*(np.arange(layers)+.5)/layers;seed[1::2]=total/layers;seeds.append(seed)
 best=None
 for j,seed in enumerate(seeds):
  start=time.monotonic();opt=minimize(objective,seed,jac=True,method='L-BFGS-B',bounds=[(-8,8)]*len(seed),options={'maxiter':220,'ftol':1e-11,'gtol':1e-6,'maxls':30})
  s=state(opt.x).ravel();fid=float(abs(np.vdot(exact,s))**2)
  print(layers,j,'E',opt.fun*100+24*10.5794,'fid',fid,'iters',opt.nit,'sec',time.monotonic()-start,flush=True)
  if best is None or opt.fun<best.fun:best=opt
 previous=best.x;s=state(best.x).ravel()
 entry=dict(layers=layers,theta=best.x.tolist(),energy_MeV=best.fun*100+24*10.5794,fidelity=float(abs(np.vdot(exact,s))**2),success=bool(best.success),message=str(best.message),gradient_norm=float(np.linalg.norm(best.jac)))
 results.append(entry); np.save(P/f'unitary_state_{layers}.npy',s)
 (P/'unitary_prep_results.json').write_text(json.dumps(results,indent=2)+'\n')
