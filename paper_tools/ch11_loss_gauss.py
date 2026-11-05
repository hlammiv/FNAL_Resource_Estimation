"""Ch. 11 photon loss vs Gauss's law, one plaquette (reproduces app02_qudit_platform.tex, first section).

Reduces the plaquette to the pair (g1, u); physical states are class functions of the holonomy g1*u.
Applies single-photon loss a|n> = sqrt(n)|n-1> on link g1 with a given level->element assignment and
returns the exact min and max over physical states of ||P a psi||^2 / ||a psi||^2, i.e. the probability
that a Gauss-law check misses the loss.  Expected output: Z_8 in group order 0.811 for all states;
2T max 0.866 for every assignment, min from ~1e-9 to ~0.023 depending on assignment.
"""
import numpy as np, itertools, scipy.linalg as sl
rng=np.random.default_rng(7)
def bt():
    q=[]
    for i in range(4):
        for s in (1,-1):
            v=[0,0,0,0]; v[i]=s; q.append(tuple(v))
    for s in itertools.product((0.5,-0.5),repeat=4): q.append(s)
    def qm(a,b):
        a0,a1,a2,a3=a; b0,b1,b2,b3=b
        return (a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0)
    key=lambda x:tuple(round(c*2) for c in x); idx={key(x):i for i,x in enumerate(q)}
    M=np.array([[idx[key(qm(a,b))] for b in q] for a in q]); I=np.array([idx[key((a[0],-a[1],-a[2],-a[3]))] for a in q])
    return M,I
def zn(N):
    M=np.array([[(a+b)%N for b in range(N)] for a in range(N)]); I=np.array([(-a)%N for a in range(N)]); return M,I
def analyze(M,I,el_at):
    n=len(M)
    cls=-np.ones(n,int); c=0
    for g in range(n):
        if cls[g]<0:
            for h in range(n): cls[M[M[h,g],I[h]]]=c
            c+=1
    # reduced space (g1,u), hol=g1 u
    hol=M[np.arange(n)[:,None],np.arange(n)[None,:]]
    B=np.zeros((n*n,c))
    for k in range(c): B[:,k]=(cls[hol]==k).ravel()
    B/=np.linalg.norm(B,axis=0)          # orthonormal physical basis
    A=np.zeros((n*n,n*n))                 # loss on g1
    lev=np.empty(n,int); lev[el_at]=np.arange(n)
    for L in range(1,n):
        for u in range(n):
            A[el_at[L-1]*n+u, el_at[L]*n+u]=np.sqrt(L)
    AB=A@B; PA=B.T@AB
    num=PA.T@PA; den=AB.T@AB
    w=sl.eigh(num,den,eigvals_only=True)
    return w.min(),w.max(),c
M,I=zn(8); print("Z8 order",analyze(M,I,np.arange(8)))
for t in range(3): print("Z8 rand",analyze(M,I,rng.permutation(8)))
M,I=bt()
mins=[];maxs=[]
r=analyze(M,I,np.arange(24)); print("2T list order",r)
for t in range(300):
    a,b,_=analyze(M,I,rng.permutation(24)); mins.append(a); maxs.append(b)
print("2T random assignments: min over states range",min(mins),max(mins)," max over states range",min(maxs),max(maxs))
