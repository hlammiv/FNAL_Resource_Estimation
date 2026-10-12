"""Accepted Chapter 2 benchmarks versus independently produced study artifacts.
The superseded model/tests are preserved in editorial_review/ch02_integration_20261004/before.
"""
import json
import math
from pathlib import Path
import pytest
from estimates import ch02_nu_nucleus as m

ROOT=Path(__file__).resolve().parents[3]            # the report tree (editorial_review/...)
MIRROR=Path(__file__).resolve().parents[2]/'applications'/'ch02_neutrino_nucleus'   # the resource repository's copy

def read(path):
    """Study files live under editorial_review/ in the report tree and under applications/ch02_neutrino_nucleus/
    in the resource repository (mirrored 2026-10-12); read whichever is present."""
    for root,rel in ((ROOT,path),(MIRROR,path.replace('editorial_review/',''))):
        if (root/rel).is_file():return json.loads((root/rel).read_text())
    raise FileNotFoundError(str(ROOT/path))

@pytest.mark.parametrize('steps',[2,6])
def test_2028_matches_independent_circuit_recount(steps):
    x=read(f'editorial_review/neutrino_feasibility/larger_circuit_counts_steps{steps}.json')
    c=x['complete_counts']
    assert m.counts('2028',steps)==(c['rotations'],c['Toffoli'],c['fixed_T'])
    assert m.circuit(m.Assumptions(),'2028',steps,'study')[0]==pytest.approx(x['T_RS_leading_estimate'])

def test_2033_matches_construction_and_loader_checks():
    x=read('editorial_review/2033_model_study/triton_resources.json')
    for row in x['rows']:
        assert m.counts('2033',row['steps'])[:2]==(row['rotations'],x['total_toffoli_allowance'])
        assert m.circuit(m.Assumptions(),'2033',row['steps'],'study')[0]==pytest.approx(row['leading_T'])
    assert m.model(m.Assumptions(),'2033').lq==(x['total_qubits'],)*2
    assert x['validation']['loader_amplitude_error']<1e-12

def test_fermionic_reference_and_observable_checks():
    x=read('editorial_review/2033_model_study/triton_crosscheck.json')
    assert x['energy_error']<1e-9
    assert max(x['correlator_errors'])<1e-9
    d=read('editorial_review/2033_model_study/triton_results.json')
    e={r['dt']:r['bias'] for r in d['rows'] if r['tau']==.2}
    assert 3.9<e[.002]/e[.001]<4.1
    assert d['L8_energy_same_couplings']-d['ground_energy']>2 # do not erase the volume limitation

@pytest.mark.parametrize('era,limit',[('2028',1e5),('2033',1e9)])
def test_complete_counts_shots_and_no_invented_schedule(era,limit):
    r=m.model(m.Assumptions(),era);i=r.intermediates
    assert r.breakdown_total()==pytest.approx(r.hard_ops[0])
    assert r.hard_ops[0]<limit and i['study_T_per_shot']<limit
    assert r.shots==(164424,164424)
    assert i['shots_per_estimate']==27404
    assert i['complete_circuit'] and not i['exact_synthesis']
    assert i['t_depth_per_shot'] is None and i['f_star'] is None and i['depth_status']=='serial_t_count_bound'
    assert r.wall_time_s==i['wall_campaign_s']
    if era=='2033':
        assert i['first_result_shots']==7012
        assert i['first_result_T']==pytest.approx(7012*i['costs_by_time'][1])
        assert i['campaign_T']==pytest.approx(2*27404*sum(i['costs_by_time']))

def test_tighter_synthesis_budget_increases_count():
    a=m.Assumptions();b=m.Assumptions(eps_syn=m.Assumed(.0025))
    assert m.model(b,'2033').hard_ops[0]>m.model(a,'2033').hard_ops[0]
    c=m.Assumptions(coherent_synthesis_allowance=m.Assumed(.0025))   # study convention only
    assert m.model(c,'2033').hard_ops[0]==m.model(a,'2033').hard_ops[0]
    assert m.model(c,'2033').intermediates['study_T_per_shot']>m.model(a,'2033').intermediates['study_T_per_shot']


def test_ruling_2026_10_12_report_synthesis_convention_and_serial_walls():
    """Ruling 2026-10-12 (Ch. 2 items 1 and 2): rotations at the report's randomized BRS fit with eps_syn = 1e-2
    per circuit; walls are serial T-count bounds (1 us per T + 0.1 ms per shot, one machine)."""
    from estimates.common import t_per_rotation, eps_per_rotation
    a=m.Assumptions()
    r28=m.model(a,'2028');i28=r28.intermediates
    r33=m.model(a,'2033');i33=r33.intermediates
    for r,i in ((r28,i28),(r33,i33)):
        nr=i['rotations']
        assert i['eps_rotation']==pytest.approx(eps_per_rotation(1e-2,nr,'incoherent'))
        assert i['t_per_rotation']==pytest.approx(t_per_rotation(i['eps_rotation'],'rus'))
        assert r.hard_ops[0]==pytest.approx(nr*i['t_per_rotation']+7*i['toffolis']+i['fixed_T'])
        assert i['study_T_per_rotation']==pytest.approx(3*math.log2(nr/.005))
        assert i['coherent_T_per_shot']>r.hard_ops[0]
    # headlines and the study ratio (box: 3.89e4 and 2.07e8; study 8.03e4 and 7.23e8)
    assert r28.hard_ops[0]==pytest.approx(3.8948e4,rel=1e-4) and i28['study_T_per_shot']==pytest.approx(8.0302e4,rel=1e-4)
    assert r33.hard_ops[0]==pytest.approx(2.0694e8,rel=1e-4) and i33['study_T_per_shot']==pytest.approx(7.2310e8,rel=1e-4)
    assert i28['study_over_report']==pytest.approx(2.06,abs=.01) and i33['study_over_report']==pytest.approx(3.49,abs=.01)
    assert i28['coherent_T_per_shot']==pytest.approx(5.05e4,rel=2e-3) and i33['coherent_T_per_shot']==pytest.approx(3.41e8,rel=2e-3)
    # per-tau tables (chapter): 2.68e4, 3.29e4, 3.89e4 and 7.04e7, 1.15e8, 2.07e8
    assert [round(c/1e4,2) for c in i28['costs_by_time']]==[2.68,3.29,3.89]
    assert [round(c/1e7,1) for c in i33['costs_by_time']]==[7.0,11.5,20.7]
    # campaign volumes, first result, eps_l
    assert i28['campaign_T']==pytest.approx(164424*r28.hard_ops[0]) and round(i28['campaign_T']/1e9,2)==6.40
    assert round(i33['first_result_T']/1e11,2)==8.09 and round(i33['campaign_T']/1e13,2)==2.15
    assert round(r28.epsilon_l[0]/1e-6,2)==2.57 and round(r33.epsilon_l[0]/1e-10,2)==4.83
    # serial walls: 2028 campaign ~1.8 h; 2033 per shot ~3.4 min, first result ~9.4 d, campaign ~250 d (0.68 yr)
    assert i28['wall_per_shot_s']==pytest.approx(r28.hard_ops[0]*1e-6+1e-4)
    assert round(r28.wall_time_s[0]/3600,1)==1.8
    assert round(i33['wall_per_shot_s']/60,1)==3.4
    assert round(i33['wall_first_result_s'][0]/86400,1)==9.4
    assert round(r33.wall_time_s[0]/86400)==249 and round(r33.wall_time_s[0]/3.15576e7,2)==0.68
