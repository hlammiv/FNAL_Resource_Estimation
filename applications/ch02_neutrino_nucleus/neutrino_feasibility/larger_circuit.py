"""Construct a coordinate-wise compact circuit for the 3^3 model.
24 position bits + 11 reusable clean work bits + 1 Hadamard ancilla.
Exports arbitrary-rotation QASM and exact primitive counts, NOT compiled T circuits.
"""
from pathlib import Path
import json,itertools,sys,argparse
import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit,qasm2
from qiskit.quantum_info import Operator
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/'scripts'))
from estimates.common import t_per_rotation,eps_per_rotation
T=10.5794;C0=-98.2265511;D0=127.839693

def ry(c,theta,q,fixed=False):
 c.ry(float(theta),q)
 if hasattr(c,'tally'):c.tally['fixed_T' if fixed else 'rotations']+=1

def phase(c,theta,q,control=None):
 if control is None:
  c.p(float(theta),q)
  if hasattr(c,'tally'):c.tally['rotations']+=1
 else:
  c.cp(float(theta),control,q)
  if hasattr(c,'tally'):c.tally['rotations']+=3

def toff(c,a,b,d):
 c.ccx(a,b,d)
 if hasattr(c,'tally'):c.tally['Toffoli']+=1

def prepare3(c,lo,hi,inverse=False):
 # |00> -> (|00>+|01>+|10>)/sqrt(3). CH uses two exact pi/4 rotations.
 theta=2*np.arcsin(1/np.sqrt(3))
 if not inverse:ry(c,theta,hi)
 c.x(hi);ry(c,np.pi/4,lo,True);c.cx(hi,lo);ry(c,-np.pi/4,lo,True);c.x(hi)
 if inverse:ry(c,-theta,hi)

def kinetic_coord(c,lo,hi,dt,control=None,work=24):
 prepare3(c,lo,hi,True);c.x(lo);c.x(hi)
 if control is None:phase(c,3*T*dt,lo,hi)
 else:
  toff(c,lo,hi,work);phase(c,3*T*dt,work,control);toff(c,lo,hi,work)
 c.x(lo);c.x(hi);prepare3(c,lo,hi)
 # -t I phase on the valid qutrit: needed only when the full evolution is controlled.

def kinetic(c,dt,control=None):
 for j in range(12):kinetic_coord(c,2*j,2*j+1,dt,control)
 if control is not None:phase(c,-12*T*dt,control)

def contact(c,particles,angle,control=None):
 first=particles[0];comparisons=[(6*first+b,6*f+b) for f in particles[1:] for b in range(6)]
 checks=[b for a,b in comparisons]
 for a,b in comparisons:c.cx(a,b);c.x(b)
 work=list(range(24,24+len(checks)-1))
 toff(c,checks[0],checks[1],work[0])
 for i in range(2,len(checks)):toff(c,work[i-2],checks[i],work[i-1])
 phase(c,angle,work[-1],control)
 for i in reversed(range(2,len(checks))):toff(c,work[i-2],checks[i],work[i-1])
 toff(c,checks[0],checks[1],work[0])
 for a,b in reversed(comparisons):c.x(b);c.cx(a,b)

def potential(c,dt,control=None):
 for pair in itertools.combinations(range(4),2):contact(c,pair,-C0*dt,control)
 for triple in itertools.combinations(range(4),3):contact(c,triple,-D0*dt,control)

# Independent local operator validation including controlled branch and clean work return.
checks={}
u=QuantumCircuit(2);prepare3(u,0,1)
checks['uniform3_error']=float(np.linalg.norm(Operator(u).data[:,0]-np.array([1,1,1,0])/np.sqrt(3)))
for controlled in [False,True]:
 c=QuantumCircuit(4 if controlled else 2);kinetic_coord(c,0,1,.003,2 if controlled else None,3)
 if controlled:phase(c,-T*.003,2)
 unit=Operator(c).data
 # un-controlled omits global exp(-it dt); compare accordingly.
 H=-T*(np.ones((3,3))-np.eye(3)); target=expm(-1j*.003*H)
 if not controlled:target*=np.exp(1j*T*.003)
 inds=np.arange(3)+(4 if controlled else 0)
 err=np.linalg.norm(unit[np.ix_(inds,inds)]-target)
 checks['controlled_kinetic_error' if controlled else 'kinetic_error']=float(err)
 if controlled:
  checks['control_off_error']=float(np.linalg.norm(unit[:3,:3]-np.eye(3)))
  checks['work_leakage']=float(np.linalg.norm(unit[8:,inds]))
assert max(checks.values())<1e-12,checks
# Reversible comparison AND network is exhaustively Boolean checked for legal 27-site triples.
for vals in itertools.product(range(27),repeat=3):
 coords=[(v%3,(v//3)%3,v//9) for v in vals]
 encoded=[sum(x<<(2*i) for i,x in enumerate(p)) for p in coords]
 flags=[not((encoded[0]^b)>>i&1) for b in encoded[1:] for i in range(6)]
 assert all(flags)==(vals[0]==vals[1]==vals[2])
parser=argparse.ArgumentParser();parser.add_argument('--steps',type=int,default=2);parser.add_argument('--tau',type=float,default=.01);args=parser.parse_args()
prep=json.loads((P/'unitary_prep_results.json').read_text())[-1];layers=prep['layers']
c=QuantumCircuit(36);c.tally=dict(rotations=0,fixed_T=0,Toffoli=0)
for j in range(12):prepare3(c,2*j,2*j+1)
initial=dict(c.tally)
for a,b in np.array(prep['theta']).reshape(-1,2):potential(c,a/100);kinetic(c,b/100)
afterprep=dict(c.tally)
# rho(q)/2 = (U_p0 + U_p1)/2; each of four pair Hadamard circuits has same count.
anc=35;q=2*np.pi/3;tau=args.tau;steps=args.steps
c.h(anc)
for bit in range(2):phase(c,q*(2**bit),bit,anc)
kinetic(c,tau/(2*steps),anc)
for j in range(steps):
 potential(c,tau/steps,anc);kinetic(c,tau/steps if j<steps-1 else tau/(2*steps),anc)
for bit in range(2):phase(c,-q*(2**bit),bit,anc)
E=json.loads((P/'larger_box_results.json').read_text())['ground_energy_MeV']-24*T
phase(c,E*tau,anc); c.h(anc);c.measure_all()
qasm2.dump(c,P/f'larger_36q_steps{steps}_real_measurement_unsynthesized.qasm')
r=c.tally['rotations'];epsilon=eps_per_rotation(.005,r,'coherent');fixed=7*c.tally['Toffoli']+c.tally['fixed_T']
out=dict(preparation_layers=layers,evolution_steps=steps,tau_MeV_inv=tau,qubits=36,initial_counts=initial,after_preparation=afterprep,complete_counts=c.tally,arbitrary_rotation_tolerance=epsilon,T_RS_leading_estimate=fixed+r*t_per_rotation(epsilon,'rs'),T_RS_plus10_per_rotation=fixed+r*(10+t_per_rotation(epsilon,'rs')),T_RUS_expected=fixed+r*t_per_rotation(eps_per_rotation(.01,r),'rus'),checks=checks,rotation_synthesis_pending=True)
(P/f'larger_circuit_counts_steps{steps}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
