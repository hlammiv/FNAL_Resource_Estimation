"""Three mobile nucleons with spin-dependent contacts; finite-volume benchmark.

Particles: proton 0, neutrons 1 and 2; spin basis labels the down-spin particle.
Two relative positions, r_n1-r_p and r_n2-r_p. Neutron exchange is antisymmetric.
All reported energies include the positive lattice kinetic diagonal.
"""
from pathlib import Path
import json, argparse
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parent

class Triton:
    def __init__(self,L,ct,cs,D=0,q=0):
        self.L=L; self.m=L**3; self.n=3*self.m**2; self.D=D; self.q=q
        self.t=10.5794; self.ct=ct; self.cs=cs
        self.xyz=np.array(np.unravel_index(np.arange(self.m),(L,L,L))).T
        self.shifts=[]
        for axis in range(3):
            for sign in (-1,1):
                dest=self.xyz.copy();dest[:,axis]=(dest[:,axis]+sign)%L
                self.shifts.append((axis,sign,np.ravel_multi_index(dest.T,(L,L,L))))
        self.pair=[]
        for i,j in ((0,1),(0,2),(1,2)):
            swap=np.eye(3);swap[[i,j]]=swap[[j,i]]
            ps=(np.eye(3)-swap)/2;pt=(np.eye(3)+swap)/2
            self.pair.append(cs*ps+(ct*pt if i==0 else 0))
        self.w=np.array([0.,1.,-1.])/np.sqrt(2)
        self.diag=np.arange(self.m)
        kk=2*np.pi*self.xyz/L
        e1=2*self.t*np.sum(1-np.cos(kk),axis=1)
        total=np.array([q,0.,0.])-kk[:,None,:]-kk[None,:,:]
        self.ek=e1[:,None]+e1[None,:]+2*self.t*np.sum(1-np.cos(total),axis=2)
        self.op=LinearOperator((self.n,self.n),matvec=self.mv,dtype=np.complex128 if q else np.float64)

    def project(self,v):
        v=v.reshape(3,self.m,self.m)
        return ((v-v[[0,2,1]].transpose(0,2,1))/2).ravel()

    def mv(self,v):
        v=v.reshape(3,self.m,self.m)
        y=18*self.t*v.copy()
        for axis,sign,p in self.shifts:
            y-=self.t*v[:,p,:]
            y-=self.t*v[:,:,p]
            phase=np.exp(-1j*self.q*sign) if self.q and axis==0 else 1.
            y-=self.t*phase*v[:,p,:][:,:,p]
        y[:,0,:]+=self.pair[0]@v[:,0,:]
        y[:,:,0]+=self.pair[1]@v[:,:,0]
        y[:,self.diag,self.diag]+=self.pair[2]@v[:,self.diag,self.diag]
        y[:,0,0]+=self.D*v[:,0,0]
        return y.ravel()

    def ground(self,v0=None):
        if v0 is None:
            v0=np.broadcast_to(self.w[:,None,None],(3,self.m,self.m)).copy().ravel()
        e,v=eigsh(self.op,k=1,which='SA',v0=self.project(v0),tol=2e-10)
        psi=self.project(v[:,0]);psi/=np.linalg.norm(psi)
        residual=np.linalg.norm(self.mv(psi)-e[0]*psi)
        assert residual<1e-6, residual
        return float(e[0]),psi,float(residual)

    def kinetic(self,v,dt):
        shape=(3,)+(self.L,)*6
        f=np.fft.fftn(v.reshape(shape),axes=tuple(range(1,7)))
        f*=np.exp(-1j*dt*self.ek).reshape((1,)+(self.L,)*6)
        return np.fft.ifftn(f,axes=tuple(range(1,7))).reshape(3,self.m,self.m)

    def potential(self,v,dt):
        v=v.copy(); origin=v[:,0,0].copy()
        for pair,where in zip(self.pair,(0,1,2)):
            e,u=np.linalg.eigh(pair); gate=(u*np.exp(-1j*dt*e))@u.T
            if where==0: v[:,0,1:]=gate@v[:,0,1:]
            elif where==1: v[:,1:,0]=gate@v[:,1:,0]
            else: v[:,self.diag[1:],self.diag[1:]]=gate@v[:,self.diag[1:],self.diag[1:]]
        e,u=np.linalg.eigh(sum(self.pair)+self.D*np.eye(3))
        v[:,0,0]=(u*np.exp(-1j*dt*e))@u.T@origin
        return v

    def step(self,v,dt):
        return self.potential(self.kinetic(self.potential(v,dt/2),dt),dt/2)

def main():
    params=json.loads((OUT/'two_body_results.json').read_text())['parameters']
    ct,cs=params['triplet_C_MeV'],params['singlet_C_MeV']
    model=Triton(6,ct,cs)
    last=None; fits=[]
    def objective(D):
        nonlocal last
        model.D=D;e,last,res=model.ground(last)
        fits.append(dict(D=D,E=e,residual=res))
        print('fit',D,e,flush=True)
        return e+8.4818
    D=brentq(objective,-100,200,xtol=1e-7)
    model.D=D; e,psi,res=model.ground(last)
    qmodel=Triton(6,ct,cs,D,q=2*np.pi/6)
    eq,psiq,resq=qmodel.ground(psi)
    elastic=float(abs(np.vdot(psiq,psi))**2)
    # Independent real-space versus Fourier kinetic action check.
    rng=np.random.default_rng(810)
    test=qmodel.project(rng.normal(size=model.n)+1j*rng.normal(size=model.n))
    test/=np.linalg.norm(test)
    shape=(3,)+(6,)*6
    fftk=np.fft.ifftn(np.fft.fftn(test.reshape(shape),axes=tuple(range(1,7)))*qmodel.ek.reshape((1,)+(6,)*6),axes=tuple(range(1,7))).ravel()
    saved=qmodel.pair;qmodel.pair=[np.zeros((3,3))]*3;savedD=qmodel.D;qmodel.D=0
    kcheck=float(np.linalg.norm(fftk-qmodel.mv(test)))
    qmodel.pair=saved;qmodel.D=savedD
    assert kcheck<1e-10,kcheck
    # Krylov projection supplies independently converged exact-reference correlators.
    def lanczos(n):
        start=psi.astype(complex); Q=[];alphas=[];betas=[];v=start.copy();prev=np.zeros_like(v);beta=0
        for j in range(n):
            Q.append(v.copy());r=qmodel.mv(v)-beta*prev;alpha=np.vdot(v,r).real;r-=alpha*v
            # Full reorthogonalization; only a small reference Krylov basis is needed.
            for old in Q: r-=np.vdot(old,r)*old
            b=np.linalg.norm(r);alphas.append(alpha)
            if b<1e-11 or j==n-1: break
            betas.append(b);prev,v=v,r/b;beta=b
        tri=np.diag(alphas)+np.diag(betas,1)+np.diag(betas,-1)
        energies,U=np.linalg.eigh(tri)
        return energies,U[0]**2
    en,w=lanczos(100);en2,w2=lanczos(140)
    times=[.025,.05,.1,.2]
    reference=lambda tau: np.sum(w2*np.exp(-1j*(en2-e)*tau))
    krylov_error=max(abs(reference(tau)-np.sum(w*np.exp(-1j*(en-e)*tau))) for tau in times)
    assert krylov_error<1e-8,krylov_error
    rows=[]
    for dt in (.004,.002,.001):
        state=psi.astype(complex).reshape(3,216,216)
        for step in range(1,round(.2/dt)+1):
            state=qmodel.step(state,dt)
            tau=step*dt
            if any(abs(tau-target)<1e-9 for target in times):
                exact=reference(tau);approx=np.vdot(psi,state.ravel())*np.exp(1j*e*tau)
                pole=elastic*np.exp(-1j*(eq-e)*tau)
                rows.append(dict(tau=tau,dt=dt,steps=step,exact=[exact.real,exact.imag],
                                 approximate=[approx.real,approx.imag],bias=float(abs(approx-exact)),
                                 inelastic_exact=[(exact-pole).real,(exact-pole).imag],
                                 norm=float(np.linalg.norm(state))))
        print('evolution finished dt',dt,flush=True)
    small=Triton(4,ct,cs,D);e4,_,r4=small.ground()
    larger=Triton(8,ct,cs,D);e8,_,r8=larger.ground()
    result=dict(model='mobile pnn, two spin channels plus three-body contact',L=6,system_qubits=864,
                parameters=dict(t=model.t,Ct=ct,Cs=cs,D=D,a_fm=1.4),
                calibration='D fitted to E=-8.4818 MeV IN THE L=6 BOX, not an infinite-volume matching',
                fit_history=fits,ground_energy=e,ground_residual=res,moving_ground_energy=eq,
                moving_ground_residual=resq,elastic_weight=elastic,inelastic_strength=1-elastic,
                L4_energy_same_couplings=e4,L8_energy_same_couplings=e8,L8_ground_residual=r8,
                kinetic_check=kcheck,krylov_convergence=krylov_error,rows=rows)
    (OUT/'triton_results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(OUT/'triton_states.npz',ground=psi,moving_ground=psiq,krylov_energies=en2,krylov_weights=w2)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__': main()
