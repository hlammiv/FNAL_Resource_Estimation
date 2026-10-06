"""Accepted Chapter 2 benchmarks versus independently produced study artifacts.
The superseded model/tests are preserved in editorial_review/ch02_integration_20261004/before.
"""
import json
from pathlib import Path
import pytest
from estimates import ch02_nu_nucleus as m

ROOT=Path(__file__).resolve().parents[3]

def read(path):return json.loads((ROOT/path).read_text())

@pytest.mark.parametrize('steps',[2,6])
def test_2028_matches_independent_circuit_recount(steps):
    x=read(f'editorial_review/neutrino_feasibility/larger_circuit_counts_steps{steps}.json')
    c=x['complete_counts']
    assert m.counts('2028',steps)==(c['rotations'],c['Toffoli'],c['fixed_T'])
    assert m.circuit(m.Assumptions(),'2028',steps)[0]==pytest.approx(x['T_RS_leading_estimate'])

def test_2033_matches_construction_and_loader_checks():
    x=read('editorial_review/2033_model_study/triton_resources.json')
    for row in x['rows']:
        assert m.counts('2033',row['steps'])[:2]==(row['rotations'],x['total_toffoli_allowance'])
        assert m.circuit(m.Assumptions(),'2033',row['steps'])[0]==pytest.approx(row['leading_T'])
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
    assert i['plus10_T']<limit
    assert r.shots==(164424,164424)
    assert i['shots_per_estimate']==27404
    assert i['complete_circuit'] and not i['exact_synthesis']
    assert r.wall_time_s is None and i['t_depth_per_shot'] is None
    if era=='2033':
        assert i['first_result_shots']==7012
        assert i['first_result_T']==pytest.approx(2.79989774790415e12)
        assert i['campaign_T']==pytest.approx(7.474327097304122e13)

def test_tighter_synthesis_allowance_increases_count():
    a=m.Assumptions();b=m.Assumptions(coherent_synthesis_allowance=m.Assumed(.0025))
    assert m.model(b,'2033').hard_ops[0]>m.model(a,'2033').hard_ops[0]
