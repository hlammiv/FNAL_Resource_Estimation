"""Build and validate one full 13-qubit Hadamard-test circuit; export QASM.
Requires study.py outputs. Arbitrary rotations remain unsynthesized.
"""
from pathlib import Path
import json,itertools
import numpy as np
from qiskit import QuantumCircuit,transpile,qasm2
from qiskit.circuit.library import UCRYGate
from qiskit.quantum_info import Statevector
P=Path(__file__).resolve().parent
r=json.loads((P/'results.json').read_text()); phi=np.load(P/'relative_state.npy'); angles=np.load(P/'loader_multiplexor_angles.npy')
prep=QuantumCircuit(13); offset=0
for level in range(9):
 target=8-level; n=2**level
 prep.append(UCRYGate(angles[offset:offset+n].tolist()),[target]+list(range(target+1,9)))
 offset+=n
for b in range(3):
 prep.h(9+b)
 for f in range(3): prep.cx(9+b,3*f+b)
state=Statevector.from_instruction(prep).data[:4096]
ix=np.arange(4096); pos=np.array([(ix>>(3*f))&7 for f in range(4)])
rel=pos[:3]^pos[3]; expected=phi[rel[0]+8*rel[1]+64*rel[2]]/np.sqrt(8)
loader_err=float(np.linalg.norm(state-expected)); assert loader_err<1e-10
V=np.zeros(4096)
for a,b in itertools.combinations(range(4),2): V+=r['C0_MeV']*(pos[a]==pos[b])
for a,b,c in itertools.combinations(range(4),3): V+=r['D0_MeV']*((pos[a]==pos[b])&(pos[b]==pos[c]))
coef=V.copy(); h=1
while h<len(coef):
 z=coef.reshape(-1,2*h); a=z[:,:h].copy(); b=z[:,h:].copy(); z[:,:h]=a+b; z[:,h:]=a-b; h*=2
coef/=4096
# Explicit parity networks + controlled rotations. A controlled Rz is two Rz and two CNOTs.
def kinetic(c,dt):
 for q in range(12):
  c.h(q); c.crz(-4*r['t_MeV']*dt,12,q); c.h(q)
def potential(c,dt):
 c.rz(-coef[0]*dt,12) # controlled identity phase, up to global phase
 for mask in np.flatnonzero(np.abs(coef[1:])>1e-10)+1:
  qs=[q for q in range(12) if mask&(1<<q)]; target=qs[-1]
  for q in qs[:-1]: c.cx(q,target)
  c.crz(2*coef[mask]*dt,12,target)
  for q in reversed(qs[:-1]): c.cx(q,target)
tau=.01; steps=2; dt=tau/steps
vals=[]; counts=None
for za,zb in itertools.product([0,3],repeat=2):
 c=prep.copy(); c.h(12); c.cz(12,zb); kinetic(c,dt/2)
 for s in range(steps):
  potential(c,dt); kinetic(c,dt if s<steps-1 else dt/2)
 c.cz(12,za); c.rz(r['energy_without_kinetic_constant_MeV']*tau,12)
 sv=Statevector.from_instruction(c).data
 vals.append(2*np.vdot(sv[:4096],sv[4096:]))
 if counts is None:
  flat=transpile(c,basis_gates=['rx','ry','rz','h','s','sdg','cx','x'],optimization_level=0)
  counts=dict(flat.count_ops()); qasm2.dump(flat,P/'hadamard_00_unsynthesized.qasm')
  # Z readout after H measures Re; sdg,H measures Im. QASM exports the Re variant.
  flat.h(12); flat.measure_all(); qasm2.dump(flat,P/'hadamard_00_real_measurement.qasm')
measured=sum(vals)/4
row=next(x for x in r['rows'] if x['tau_MeV_inv']==tau and x['steps']==steps)
err=abs(measured-complex(*row['approx'])); assert err<1e-10
out=dict(loader_statevector_error=loader_err,full_Hadamard_vs_split_error=float(err),normalized_correlator=[measured.real,measured.imag],decomposed_gate_counts=counts,qubits=13,arbitrary_rotations_not_Clifford_T_synthesized=True)
(P/'circuit_check.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
