"""Independent small-box occupation-space check of the relative-spin calculation."""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
from scipy.linalg import eigh
from triton import Triton

P=Path(__file__).resolve().parent
p=json.loads((P/'triton_results.json').read_text())['parameters']
ct,cs,D,t=p['Ct'],p['Cs'],p['D'],p['t'];L=2;M=8
ann=[]
for j in range(4):
    a=np.zeros((16,16))
    for b in range(16):
        if b>>j&1:a[b^(1<<j),b]=(-1)**((b&((1<<j)-1)).bit_count())
    ann.append(a)
pair=lambda i,j:ann[j]@ann[i]
singlets=[pair(0,1),pair(2,3),(pair(0,3)-pair(1,2))/np.sqrt(2)]
triplets=[pair(0,2),pair(1,3),(pair(0,3)+pair(1,2))/np.sqrt(2)]
V=sum(cs*a.T@a for a in singlets)+sum(ct*a.T@a for a in triplets)
for b in range(16):
    n=b.bit_count();V[b,b]+=D*n*(n-1)*(n-2)/6
basis=[]
for proton in range(M):
    for nu in range(M):
        for nd in range(M):basis.append((1<<(4*proton))|(1<<(4*nu+2))|(1<<(4*nd+3)))
    for nu1,nu2 in combinations(range(M),2):
        basis.append((1<<(4*proton+1))|(1<<(4*nu1+2))|(1<<(4*nu2+2)))
lookup={b:i for i,b in enumerate(basis)};assert len(lookup)==736
H=np.eye(len(basis))*18*t
xyz=np.array(np.unravel_index(np.arange(M),(L,L,L))).T
for col,b in enumerate(basis):
    for site in range(M):
        local=(b>>(4*site))&15
        for out in np.flatnonzero(abs(V[:,local])>1e-12):
            target=(b&~(15<<(4*site)))|(int(out)<<(4*site))
            H[lookup[target],col]+=V[out,local]
        for flavor in range(4):
            j=4*site+flavor
            if not (b>>j&1):continue
            removed=b^(1<<j);sign0=(-1)**((b&((1<<j)-1)).bit_count())
            for axis in range(3):
                for direction in (-1,1):
                    dest=xyz[site].copy();dest[axis]=(dest[axis]+direction)%L
                    other=np.ravel_multi_index(dest,(L,L,L));i=int(4*other+flavor)
                    if removed>>i&1:continue
                    sign1=(-1)**((removed&((1<<i)-1)).bit_count())
                    H[lookup[removed|(1<<i)],col]-=t*sign0*sign1
assert np.max(abs(H-H.T))<1e-10
E,U=eigh(H);ground=U[:,0]
model=Triton(L,ct,cs,D);e,psi,res=model.ground()
energy_error=abs(e-E[0]);assert energy_error<1e-8
phases=[]
for b in basis:
    proton=next(site for site in range(M) if (b>>(4*site))&3)
    phases.append(np.exp(1j*np.pi*xyz[proton,0]))
source=np.array(phases)*ground
w=abs(U.T@source)**2
qmodel=Triton(L,ct,cs,D,q=np.pi)
eye=np.eye(qmodel.n)
Hq=np.column_stack([qmodel.mv(v.astype(complex)) for v in eye])
eq,Uq=eigh(Hq);wq=abs(Uq.conj().T@psi)**2
errors=[]
for tau in (.01,.05,.1,.2):
    full=np.sum(w*np.exp(-1j*(E-E[0])*tau))
    relative=np.sum(wq*np.exp(-1j*(eq-e)*tau))
    errors.append(float(abs(full-relative)))
assert max(errors)<1e-9,errors
result=dict(full_occupation_dimension=len(basis),relative_working_dimension=model.n,
            ground_energy_full=float(E[0]),ground_energy_relative=e,energy_error=float(energy_error),
            correlator_errors=errors)
(P/'triton_crosscheck.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
