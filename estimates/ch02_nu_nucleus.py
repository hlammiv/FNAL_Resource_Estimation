"""Ch. 2 — Neutrino–nucleus response benchmarks.
Reproduces applications/app01_neutrino_nucleus.tex after the accepted 2026-10-04 revision.

INPUTS AND SOURCES
2028: explicit 36-qubit circuits in editorial_review/neutrino_feasibility.
2033: 908-qubit construction in editorial_review/2033_model_study/TRITON.md.
Both include deterministic preparation and use coherent synthesis allowance .005,
leading Ross–Selinger 3 log2(N_R/.005), seven T per Toffoli.
Historical assumptions are retained only for the codesign scaling comparator.

WHAT IS NOT DERIVED HERE
Physical accuracy, exact synthesis and logical-noise behavior. The 2033 primitive
counts are constructive allowances, not a compiled full-register circuit. The
codesign row remains an evolution scaling estimate, not a complete response shot.

WALL TIME AND FACTORIES
No T-depth or device schedule has been derived for the accepted benchmarks.
Their depth/factory/wall-time exports are explicitly None, with status not_established.
Shot counts and gate volumes are reported instead. The legacy timing model applies
only to the retained codesign scaling comparison.
"""
from dataclasses import dataclass
import math
from . import ch02_legacy as legacy
from .common import Result, Published, Primitive, CircuitStatus, DEPTH_EXPORTS, Tagged, Assumed

@dataclass(frozen=True)
class Assumptions(legacy.Assumptions):
    coherent_synthesis_allowance: Tagged = Assumed(.005, "Coherent per-circuit synthesis allowance", "app01:eq:nu_synthesis")


def counts(era, steps):
    if era == '2028':
        # Independent gate recount: initial + 8 preparation layers + insertion/evolution.
        return 646+91*steps, 1208+172*steps, 456+48*steps
    if era == '2033':
        return (2**20-1+21+432+58752*steps+9504*(steps+1)+1,
                6*63*11*7+3*216*2*17+432*20+256+32, 0)
    raise ValueError(era)


def circuit(a, era, steps):
    nr,nt,nfixed=counts(era,steps)
    if not 0<a.coherent_synthesis_allowance.lo<1:
        raise ValueError('coherent synthesis allowance must be between zero and one')
    rate=3*math.log2(nr/a.coherent_synthesis_allowance.lo)
    return nr*rate+7*nt+nfixed, rate


def model(a,era):
    if era=='codesign':
        return legacy.model(a,era)
    if era not in ('2028','2033'):raise ValueError(era)
    steps=(2,4,6) if era=='2028' else (25,50,100)
    lq=36 if era=='2028' else 908
    costs=[circuit(a,era,n)[0] for n in steps]
    cost,rate=circuit(a,era,steps[-1]);nr,nt,nfixed=counts(era,steps[-1])
    repetitions=math.ceil(2*math.log(12/.05)/.02**2)
    campaign=6*repetitions
    # The 2028 manuscript intentionally quotes a conservative deepest-shot bound.
    volume=campaign*cost if era=='2028' else 2*repetitions*sum(costs)
    inter={key:None for key in DEPTH_EXPORTS}
    inter.update(t_per_shot=(cost,cost),depth_status='not_established',
                 rotations=nr,toffolis=nt,fixed_T=nfixed,eps_rotation=a.coherent_synthesis_allowance.lo/nr,
                 steps=steps,costs_by_time=costs,plus10_T=cost+10*nr,
                 shots_per_estimate=repetitions,campaign_shots=campaign,campaign_T=volume,
                 campaign_T_scope='deepest-shot upper estimate' if era=='2028' else 'sum of time-specific costs',
                 complete_circuit=True,exact_synthesis=False)
    if era=='2033':
        nfirst=2*math.ceil(2*math.log(4/.05)/.05**2)
        inter.update(first_result_shots=nfirst,first_result_T=nfirst*costs[1])
    status=CircuitStatus.COMPILED if era=='2028' else CircuitStatus.SCALING
    src='editorial_review/neutrino_feasibility/LARGER_BOX.md' if era=='2028' else 'editorial_review/2033_model_study/TRITON.md'
    breakdown=(Primitive('rotations',nr,rate,status,src,'leading synthesis estimate, subleading term pending'),
               Primitive('Toffoli',nt,7,status,src),Primitive('fixed_T',nfixed,1,status,src))
    return Result(era=era,lq=(lq,lq),hard_ops=(cost,cost),breakdown=breakdown,
                  intermediates=inter,shots=(campaign,campaign),wall_time_s=None,
                  epsilon_l=(a.fault_budget.lo/cost,)*2,
                  notes=('Complete correlator; classically checkable finite-volume model. No hardware schedule.',))

PUBLISHED={
    '2028':Published((36,36),(8.03e4,8.03e4),'app01:2028 complete four-nucleon correlator box',rel_tol=.005),
    '2033':Published((908,908),(7.23e8,7.23e8),'app01:2033 complete mobile-triton correlator box',rel_tol=.005),
    'codesign':Published((800,1000),(1.7e12,5.5e12),'app01:longer-term carbon EVOLUTION scaling comparison',rel_tol=.1),
}

def INSTANCE_ROWS(a,era,r):
    if era=='codesign':
        return [('12C converged-basis evolution scaling (incomplete circuit)',r.lq,r.hard_ops,
                 {'route':'chiral scaling','nucleus':'12C','t_note':'Evolution scaling, not complete response; no step-error validation'})]
    label='Four-nucleon charge correlator, 3^3, complete' if era=='2028' else 'Mobile triton charge response, 6^3, complete'
    return [(label,r.lq,r.hard_ops,{'nucleus':'4N' if era=='2028' else '3H',
                                 't_note':'Leading synthesis estimate; complete circuit, no hardware schedule'})]
