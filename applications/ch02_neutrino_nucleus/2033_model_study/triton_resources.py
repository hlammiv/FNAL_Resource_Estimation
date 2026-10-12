"""Conservative constructive counts and local/preparation checks for triton."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from qiskit.quantum_info import SparsePauliOp

OUT=Path(__file__).resolve().parent
data=json.loads((OUT/'triton_results.json').read_text())
state=np.load(OUT/'triton_states.npz')['ground'].reshape(3,216,216)
p=data['parameters'];ct,cs,D,t=p['Ct'],p['Cs'],p['D'],p['t']

# Full four-mode on-site Hamiltonian, with fermionic signs.
ann=[]
for j in range(4):
    op=np.zeros((16,16))
    for b in range(16):
        if b & (1<<j):op[b^(1<<j),b]=(-1)**((b & ((1<<j)-1)).bit_count())
    ann.append(op)
pair=lambda i,j:ann[j]@ann[i]
singlets=[pair(0,1),pair(2,3),(pair(0,3)-pair(1,2))/np.sqrt(2)]
triplets=[pair(0,2),pair(1,3),(pair(0,3)+pair(1,2))/np.sqrt(2)]
V=sum(cs*x.T@x for x in singlets)+sum(ct*x.T@x for x in triplets)
for b in range(16):
    n=b.bit_count();V[b,b]+=D*(n*(n-1)*(n-2)//6)
terms=[(s,float(c.real)) for s,c in SparsePauliOp.from_operator(V).simplify().to_list() if s!='IIII']
diag=np.diag(np.diag(V));off=V-diag
assert np.max(abs(diag@off-off@diag))<1e-9
offterms=[SparsePauliOp(s).to_matrix() for s,c in terms if 'X' in s or 'Y' in s]
for x in offterms:
    for y in offterms:assert np.max(abs(x@y-y@x))<1e-12
assert abs(V[13,13]-(1.5*(ct+cs)+D))<1e-10

# Exact real-state loader on 20 bits: two 9-bit relative positions, two spin bits.
xyz=np.array(np.unravel_index(np.arange(216),(6,6,6))).T
code=xyz[:,0]*64+xyz[:,1]*8+xyz[:,2]
indices=np.arange(3)[:,None,None]*512**2+code[None,:,None]*512+code[None,None,:]
amplitudes=np.zeros(2**20);amplitudes[indices.ravel()]=state.real.ravel()
assert np.linalg.norm(state.imag if np.iscomplexobj(state) else 0)<1e-10
angle_levels=[];layer=amplitudes.copy()
for level in range(20):
    pairs=layer.reshape(-1,2)
    angle_levels.append(2*np.arctan2(pairs[:,1],pairs[:,0]))
    layer=np.linalg.norm(pairs,axis=1)
reconstructed=layer.copy()
for theta in angle_levels[::-1]:
    reconstructed=np.stack((reconstructed*np.cos(theta/2),reconstructed*np.sin(theta/2)),axis=1).ravel()
loader_error=float(np.linalg.norm(reconstructed-amplitudes))
assert loader_error<1e-12
exchange_error=float(np.linalg.norm(state+state[[0,2,1]].transpose(0,2,1)))
spin32_component=float(np.linalg.norm(state.sum(axis=0)))
assert exchange_error<1e-10
assert spin32_component<1e-8
assert np.max(abs(state[0,np.arange(216),np.arange(216)]))<1e-10
np.savez_compressed(OUT/'triton_loader_angles.npz',**{f'level_{i}':x for i,x in enumerate(angle_levels[::-1])})

# Fourier decomposition: 15 complex Givens and 6 diagonal phases per 6-mode block.
F=np.exp(2j*np.pi*np.arange(6)[:,None]*np.arange(6)[None,:]/6)/np.sqrt(6)
A=F.copy();rotations=0
for col in range(5):
    for row in range(5,col,-1):
        a,b=A[row-1,col],A[row,col];norm=np.hypot(abs(a),abs(b))
        if norm<1e-14:continue
        phase=a/abs(a) if abs(a)>1e-14 else 1.
        c=abs(a)/norm;s=phase*b.conjugate()/norm
        G=np.array([[c,s],[-s.conjugate(),c]])
        A[[row-1,row],:]=G@A[[row-1,row],:];rotations+=1
fourier_error=float(np.linalg.norm(A-np.diag(np.diag(A))))
assert fourier_error<1e-12
assert rotations<=15

# Finite-volume d+n cluster threshold, only a diagnostic, not an exclusive yield.
k=2*np.pi*xyz/6; eps=2*t*np.sum(1-np.cos(k),axis=1)
q=np.array([2*np.pi/6,0,0]);thresholds=[]
for P in k:
    free=eps+2*t*np.sum(1-np.cos(P-k),axis=1)
    f=lambda E:np.mean(1/(E-free))-1/ct
    E=brentq(f,float(free.min()-1000),float(free.min()-1e-10))
    thresholds.append(E+2*t*np.sum(1-np.cos(q-P)))

M=216
allocations=dict(occupation=864,coordinates=27,spin_labels_and_decode=5,
                 clean_arithmetic_scratch=10,neutron_sort_flag=1,hadamard=1)
# At most 63 transpositions for each permutation of a 6-bit word; each at most
# 11 five-controlled NOTs, at 7 Toffolis each. Six relative-to-absolute adds.
toff=dict(translations=6*63*11*7,occupation_write=3*216*2*17,
          neutron_label_erase=432*20,neutron_sort=256,spin_decode=32)
loader_rot=2**20-1+21 # plus three uniform six-valued COM coordinates
current_rot=432 # conservative: phase on each proton mode, including zero/Clifford angles
kinetic_rot=2*4*(3*M//6)*(15*4+6)+2*4*M
contact_rot=2*M*len(terms)
rows=[]
for steps in (25,50,100,200):
    # V-half K V ... K V-half has steps+1 contact layers.
    nrot=loader_rot+current_rot+steps*kinetic_rot+(steps+1)*contact_rot+1
    leading=float(3*nrot*np.log2(nrot/.005)+7*sum(toff.values()))
    rows.append(dict(steps=steps,tau_at_dt002=steps*.002,rotations=nrot,
                     leading_T=leading,plus_10_per_rotation_T=leading+10*nrot))
out=dict(allocated_qubits=allocations,total_qubits=sum(allocations.values()),
         loader_rotations=loader_rot,current_rotations=current_rot,
         preparation_toffoli_allowances=toff,total_toffoli_allowance=sum(toff.values()),
         kinetic_rotations_per_step=kinetic_rot,contact_rotations_per_layer=contact_rot,
         onsite_nonidentity_paulis=terms,rows=rows,
         validation=dict(loader_amplitude_error=loader_error,neutron_exchange_error=exchange_error,
                         spin_three_halves_component=spin32_component,fourier_QR_error=fourier_error,
                         finite_volume_dn_threshold_at_q=float(min(thresholds))))
(OUT/'triton_resources.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
