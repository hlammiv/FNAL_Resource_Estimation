"""0+ (or J0) spectrum with a (J^2 - J0(J0+1))^2 penalty: the lowest states are all J = J0; energies = H (penalty 0)."""
from sm import *
import json, sys
orbs={'f7p3':{(0,3,7),(1,1,3)},'f7p3p1':{(0,3,7),(1,1,3),(1,1,1)},'pf':{(0,3,7),(1,1,3),(1,1,1),(0,3,5)},
      'sd':{(0,2,5),(1,0,1),(0,2,3)}}
nucs={'48Ca':(0,8,48),'48Ti':(2,6,48),'27Al':(5,6,27)}
snt,space,nuc,M2,J0,k=sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6])
Z,N,A=nucs[nuc]
t0=time.time()
m=Model(snt+'.snt',orbs[space],Z,N,A)
dim,w,v=m.lowest(M2,k=k,tol=1e-10,penalty=(1.0,J0) if J0>=0 else None)
m.penalty=None
E=[];Js=[]
for i in range(v.shape[1]):
    Psi=np.zeros(m.mask.shape); Psi[m.mask]=v[:,i]
    E.append(float(np.sum(Psi*m.matvec_full(Psi))))
    j2=m.J2(v[:,i]); Js.append((-1+math.sqrt(1+4*j2))/2)
o=np.argsort(E); E=[E[i] for i in o]; Js=[Js[i] for i in o]
print(snt,space,nuc,'M2',M2,'dim',dim,'J0',J0,'E',[round(e,5) for e in E],'Ex',[round(e-E[0],4) for e in E],'J',[round(j,3) for j in Js],'t',round(time.time()-t0,1),flush=True)
json.dump(dict(snt=snt,space=space,nuc=nuc,M2=M2,dim=dim,J0=J0,E=E,J=Js),open(f'gap_{snt}_{space}_{nuc}_J{J0}.json','w'),indent=1)
