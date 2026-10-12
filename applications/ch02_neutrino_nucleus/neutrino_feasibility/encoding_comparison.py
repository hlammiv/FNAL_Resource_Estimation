"""Audited compact counts and explicit occupation-encoding construction comparison.
Does NOT claim optimality or a lower bound on either encoding.
"""
from pathlib import Path
import json,math,sys
from qiskit import qasm2,transpile
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/'scripts'))
from estimates.common import eps_per_rotation,t_per_rotation
rows=[]
for steps in [2,6]:
 r=json.loads((P/f'larger_circuit_counts_steps{steps}.json').read_text())
 qc=qasm2.load(P/f'larger_36q_steps{steps}_real_measurement_unsynthesized.qasm',custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
 flat=transpile(qc,basis_gates=['rx','ry','rz','p','h','s','sdg','t','tdg','cx','x','measure'],optimization_level=0)
 arbitrary=0;fixedT=0
 for entry in flat.data:
  op=entry.operation
  if op.name in ['t','tdg']:fixedT+=1
  elif op.name in ['rx','ry','rz','p']:
   units=float(op.params[0])/(math.pi/4)
   if abs(units-round(units))<1e-10:fixedT+=round(units)%2
   else:arbitrary+=1
 assert arbitrary==r['complete_counts']['rotations'],(arbitrary,r)
 assert fixedT==r['complete_counts']['fixed_T']+7*r['complete_counts']['Toffoli'],(fixedT,r)
 # Same 8-layer unitary ansatz in site occupation registers (four species, 27 sites).
 # V: 14 diagonal Pauli strings/site =>378 rotations/layer (identity phase separately).
 # K: 4 species*3 directions*9 lines=108 three-mode blocks, each F diag Fdag.
 # F is two Givens (one generic angle, one pi/4): F+Fdag gives 4 generic and 4 exact T.
 # Uncontrolled diag adds1 rotation; controlled diag adds3; Jordan-Wigner strings Clifford.
 # Initial product of uniform 27-site one-particle states: 4*(1+3+9)=52 F blocks.
 prep_rot=104+8*(378+540); prep_fixed=104+8*432
 evolution_rot=steps*(2*378+1)+(steps+1)*(756+1)
 # Two controlled current insertions, each 18 nonzero site phases *3 rotations; E0 phase1.
 total_rot=prep_rot+evolution_rot+108+1
 fixed=prep_fixed+(steps+1)*432
 eps=eps_per_rotation(.005,total_rot,'coherent')
 rows.append(dict(steps=steps,compact_qubits=36,compact_T_leading=r['T_RS_leading_estimate'],compact_T_plus10=r['T_RS_plus10_per_rotation'],compact_arbitrary_rotations_checked=arbitrary,compact_exact_T_checked=fixedT,occupation_qubits=109,occupation_rotations=total_rot,occupation_exact_T=fixed,occupation_T_leading=fixed+total_rot*t_per_rotation(eps,'rs'),occupation_T_report_RUS=fixed+total_rot*t_per_rotation(eps_per_rotation(.01,total_rot),'rus')))
(P/'encoding_comparison.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
