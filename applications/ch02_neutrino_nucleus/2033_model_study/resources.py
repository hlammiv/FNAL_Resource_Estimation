"""Construction-level counts, not compiled Clifford+T guarantees."""
import json
from pathlib import Path
import numpy as np
from qiskit.quantum_info import SparsePauliOp

out=Path(__file__).resolve().parent
data=json.loads((out/'two_body_results.json').read_text())
ct=data['parameters']['triplet_C_MeV']; cs=data['parameters']['singlet_C_MeV']
ann=[]
for j in range(4):
    op=np.zeros((16,16))
    for b in range(16):
        if b & (1<<j): op[b^(1<<j),b]=(-1)**((b & ((1<<j)-1)).bit_count())
    ann.append(op)
pair=lambda i,j: ann[j]@ann[i]
singlets=[pair(0,1),pair(2,3),(pair(0,3)-pair(1,2))/np.sqrt(2)]
triplets=[pair(0,2),pair(1,3),(pair(0,3)+pair(1,2))/np.sqrt(2)]
V=sum(cs*p.T@p for p in singlets)+sum(ct*p.T@p for p in triplets)
paulis=SparsePauliOp.from_operator(V).simplify()
terms=[(p,float(c.real)) for p,c in paulis.to_list() if p!='IIII']
diag=np.diag(np.diag(V)); off=V-diag
assert np.max(abs(diag@off-off@diag))<1e-9
offterms=[SparsePauliOp(p).to_matrix() for p,c in terms if 'X' in p or 'Y' in p]
for x in offterms:
    for y in offterms: assert np.max(abs(x@y-y@x))<1e-12
# Two particles in pp have exactly Cs on contact.
assert abs(V[3,3]-cs)<1e-10
assert abs(V[5,5]-ct)<1e-10

m=216; steps=200
# 3M bonds/flavor; 4 flavors; 2 Pauli rotations/hop; forward/backward;
# controlled Pauli rotation decomposes into 2 elementary rotations.
kinetic_rotations_per_step=3*m*4*2*2*2
contact_rotations_per_step=m*len(terms)*2
# Constants and the reference energy combine into one controlled phase per shot.
prep_rotations=2**16-1
rotations=steps*(kinetic_rotations_per_step+contact_rotations_per_step)+prep_rotations+1
# Eight controls: 6 clean scratch and 13 Toffolis. Two one-hot outputs.
toffolis=2*m*13
leading=rotations*3*np.log2(rotations/.005)+7*toffolis
results=dict(system_qubits=4*m,position_loader_qubits=16,loader_scratch=6,hadamard_ancilla=1,
             conservative_peak_qubits=4*m+16+6+1,steps=steps,tau=.5,
             onsite_pauli_terms=terms,kinetic_rotations_per_step=kinetic_rotations_per_step,
             contact_rotations_per_step=contact_rotations_per_step,
             preparation_rotations=prep_rotations,total_rotations=rotations,
             toffolis=toffolis,leading_T=float(leading),illustrative_plus_10_T_per_rotation=float(leading+10*rotations))
(out/'resources.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
