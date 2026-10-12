"""Ch. 2 — Neutrino–nucleus response benchmarks.
Reproduces applications/app01_neutrino_nucleus.tex after the accepted 2026-10-04 revision.

INPUTS AND SOURCES
2028: explicit 36-qubit circuits in editorial_review/neutrino_feasibility.
2033: 908-qubit construction in editorial_review/2033_model_study/TRITON.md.
Both include deterministic preparation and seven T per Toffoli. Since the 2026-10-12 ruling the
rotations are priced with the report convention (randomized BRS fit at eps_syn = 1e-2 per circuit,
E20 / R-TOL); the study files' coherent leading Ross-Selinger count, 3 log2(N_R/0.005), is kept as
circuit(..., convention='study') and reproduced by the tests (2.06x and 3.49x the headline).
Historical assumptions are retained only for the codesign scaling comparator.

WHAT IS NOT DERIVED HERE
Physical accuracy, exact synthesis and logical-noise behavior. The 2033 primitive
counts are constructive allowances, not a compiled full-register circuit. The
codesign row remains an evolution scaling estimate, not a complete response shot.

WALL TIME AND FACTORIES
No T-depth or device schedule has been derived for the accepted benchmarks, so the depth and
factory exports are None (status serial_t_count_bound). Since the 2026-10-12 ruling the walls are
the report's serial T-count bounds on one machine (1 us per T + 0.1 ms per shot): lower bounds,
since a depth-bound circuit on 44 workspace qubits could run several times longer.
"""
from dataclasses import dataclass
import math
from . import ch02_legacy as legacy
from .common import (Result, Published, Primitive, CircuitStatus, DEPTH_EXPORTS, Tagged, Assumed,
                     t_per_rotation, eps_per_rotation)

@dataclass(frozen=True)
class Assumptions(legacy.Assumptions):
    coherent_synthesis_allowance: Tagged = Assumed(.005, "The chapter authors' study convention (sensitivity only since the "
        "2026-10-12 ruling): coherent leading Ross-Selinger synthesis, delta_syn = 0.005 shared among N_R rotations; "
        "reproduces the study files' T_RS_leading_estimate / leading_T", "app01:eq:nu_synthesis (pre-2026-10-12)")


def counts(era, steps):
    if era == '2028':
        # Independent gate recount: initial + 8 preparation layers + insertion/evolution.
        return 646+91*steps, 1208+172*steps, 456+48*steps
    if era == '2033':
        return (2**20-1+21+432+58752*steps+9504*(steps+1)+1,
                6*63*11*7+3*216*2*17+432*20+256+32, 0)
    raise ValueError(era)


def circuit(a, era, steps, convention='report'):
    """T count of the complete circuit with `steps` evolution steps.

    convention='report'  (headline since the 2026-10-12 ruling) randomized synthesis at the report's R-TOL budget
                         eps_syn = 1e-2 per circuit, eps_R = sqrt(eps_syn/N_R), 1.15 log2(1/eps_R) + 9.2 T per rotation (E20)
    convention='coherent' the same BRS fit with coherently adding errors, eps_R = eps_syn/N_R (sensitivity)
    convention='study'   the chapter authors' leading Ross-Selinger count, 3 log2(N_R/delta) at delta = 0.005
                         (what the study files editorial_review/{neutrino_feasibility,2033_model_study} record)
    Returns (total T, T per rotation)."""
    nr,nt,nfixed=counts(era,steps)
    if convention=='study':
        if not 0<a.coherent_synthesis_allowance.lo<1:
            raise ValueError('coherent synthesis allowance must be between zero and one')
        rate=3*math.log2(nr/a.coherent_synthesis_allowance.lo)
    elif convention in ('report','coherent'):
        eps_r=eps_per_rotation(a.eps_syn.lo,nr,'incoherent' if convention=='report' else 'coherent')
        rate=t_per_rotation(eps_r,'rus')
    else:
        raise ValueError(convention)
    return nr*rate+7*nt+nfixed, rate


def model(a,era):
    if era=='codesign':
        return legacy.model(a,era)
    if era not in ('2028','2033'):raise ValueError(era)
    steps=(2,4,6) if era=='2028' else (25,50,100)
    lq=36 if era=='2028' else 908
    costs=[circuit(a,era,n)[0] for n in steps]
    cost,rate=circuit(a,era,steps[-1]);nr,nt,nfixed=counts(era,steps[-1])
    study_cost=circuit(a,era,steps[-1],'study')[0]
    coherent_cost=circuit(a,era,steps[-1],'coherent')[0]
    repetitions=math.ceil(2*math.log(12/.05)/.02**2)
    campaign=6*repetitions
    # The 2028 box quotes a conservative deepest-shot bound for the campaign volume.
    volume=campaign*cost if era=='2028' else 2*repetitions*sum(costs)
    tg=a.t_gate_s.lo; t0=a.shot_overhead_s.lo
    # Ruling 2026-10-12 (2): walls are serial T-count bounds on one machine (1 us per T + 0.1 ms per shot);
    # no T-depth or factory schedule exists, so the depth keys stay None and the walls are lower bounds.
    wall_campaign=volume*tg+campaign*t0
    inter={key:None for key in DEPTH_EXPORTS}
    inter.update(t_per_shot=(cost,cost),depth_status='serial_t_count_bound',
                 wall_campaign_s=(wall_campaign,wall_campaign),wall_per_shot_s=cost*tg+t0,
                 rotations=nr,toffolis=nt,fixed_T=nfixed,
                 eps_rotation=eps_per_rotation(a.eps_syn.lo,nr,'incoherent'),t_per_rotation=rate,
                 synthesis_convention='randomized BRS fit at eps_syn = 1e-2 per circuit (E20, R-TOL; ruling 2026-10-12)',
                 study_T_per_shot=study_cost,study_T_per_rotation=circuit(a,era,steps[-1],'study')[1],
                 coherent_T_per_shot=coherent_cost,study_over_report=study_cost/cost,
                 steps=steps,costs_by_time=costs,costs_by_time_study=[circuit(a,era,n,'study')[0] for n in steps],
                 shots_per_estimate=repetitions,campaign_shots=campaign,campaign_T=volume,
                 campaign_T_scope='deepest-shot upper estimate' if era=='2028' else 'sum of time-specific costs',
                 complete_circuit=True,exact_synthesis=False)
    if era=='2033':
        nfirst=2*math.ceil(2*math.log(4/.05)/.05**2)
        first_T=nfirst*costs[1]
        inter.update(first_result_shots=nfirst,first_result_T=first_T,
                     wall_first_result_s=(first_T*tg+nfirst*t0,)*2)
    status=CircuitStatus.COMPILED if era=='2028' else CircuitStatus.SCALING
    src='editorial_review/neutrino_feasibility/LARGER_BOX.md' if era=='2028' else 'editorial_review/2033_model_study/TRITON.md'
    breakdown=(Primitive('rotations',nr,rate,status,src,'randomized BRS synthesis at eps_R = sqrt(1e-2/N_R) (ruling 2026-10-12); '
                         'the study files price the same rotations coherently at 3 log2(N_R/0.005)'),
               Primitive('Toffoli',nt,7,status,src),Primitive('fixed_T',nfixed,1,status,src))
    return Result(era=era,lq=(lq,lq),hard_ops=(cost,cost),breakdown=breakdown,
                  intermediates=inter,shots=(campaign,campaign),wall_time_s=(wall_campaign,wall_campaign),
                  epsilon_l=(a.fault_budget.lo/cost,)*2,
                  notes=('Complete correlator; classically checkable finite-volume model. Walls are serial T-count bounds; no depth schedule.',))

PUBLISHED={
    '2028':Published((36,36),(3.89e4,3.89e4),'app01:2028 complete four-nucleon correlator box (report synthesis convention, ruling 2026-10-12; study convention 8.03e4)',rel_tol=.005),
    '2033':Published((908,908),(2.07e8,2.07e8),'app01:2033 complete mobile-triton correlator box (report synthesis convention, ruling 2026-10-12; study convention 7.23e8)',rel_tol=.005),
    'codesign':Published((800,1000),(1.7e12,5.5e12),'app01:longer-term carbon EVOLUTION scaling comparison',rel_tol=.1),
}

def INSTANCE_ROWS(a,era,r):
    if era=='codesign':
        return [('12C converged-basis evolution scaling (incomplete circuit)',r.lq,r.hard_ops,
                 {'route':'chiral scaling','nucleus':'12C','t_note':'Evolution scaling, not complete response; no step-error validation'})]
    label='Four-nucleon charge correlator, 3^3, complete' if era=='2028' else 'Mobile triton charge response, 6^3, complete'
    return [(label,r.lq,r.hard_ops,{'nucleus':'4N' if era=='2028' else '3H',
                                 't_note':'Randomized BRS synthesis at eps_syn 1e-2 (report convention); complete circuit; walls are serial T-count bounds, no depth schedule'})]
