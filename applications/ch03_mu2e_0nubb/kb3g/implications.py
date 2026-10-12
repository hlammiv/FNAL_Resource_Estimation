import sys, dataclasses, json, glob, os
sys.dont_write_bytecode = True
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
from estimates import ch03_mu2e_0nubb as m
from estimates.common import Assumed
D = str(__import__('pathlib').Path(__file__).resolve().parent)
snt = sys.argv[1] if len(sys.argv) > 1 else 'kb3g'
res = {os.path.basename(f).split('_')[2].split('.')[0]: json.load(open(f)) for f in glob.glob(f'{D}/lam_{snt}_*.json')}
a0 = m.Assumptions()
kw = {}
if 'f7p3' in res:
    kw['lam_over_m_f7p3'] = Assumed(res['f7p3']['lam_over_M'], f"COMPUTED {snt}", "x")
    kw['p_ref_f7p3'] = Assumed(res['f7p3']['p_ref'], f"COMPUTED {snt}", "x")
if 'f7p3p1' in res:
    kw['lam_over_m_f7p3p1'] = Assumed(res['f7p3p1']['lam_over_M'], f"COMPUTED {snt}", "x")
if 'pf' in res:
    r = res['pf']['lam_over_M']; kw['lam_over_m_pf'] = Assumed((r, r), f"COMPUTED {snt}", "x")
elif len(sys.argv) > 2:
    lo, hi = (float(x) for x in sys.argv[2].split(',')); kw['lam_over_m_pf'] = Assumed((lo, hi), "placeholder", "x")
a1 = dataclasses.replace(a0, **kw)
print("inputs:", {k: v.value for k, v in kw.items()})
R0, R1 = m.simplified_0nubb(a0), m.simplified_0nubb(a1)
d = 86400.0; yr = 365.25 * d
def fmt(x):
    if isinstance(x, tuple): return f"{x[0]:.3g}-{x[1]:.3g}"
    return f"{x:.3g}"
for key in R0:
    print(f"--- {key}")
    for f in ('shots', 'n_qdrift', 't', 't_depth', 'wall_s', 'eps_l'):
        if f in R0[key]:
            v0, v1 = R0[key][f], R1[key][f]
            extra = ''
            if f == 'wall_s':
                extra = f"   [{fmt(tuple(x/d for x in v0))} d -> {fmt(tuple(x/d for x in v1))} d; {fmt(tuple(x/yr for x in v1))} yr]"
            if f == 't':
                extra = f"   [x target: {fmt(tuple(x/1e9 for x in v0))} -> {fmt(tuple(x/1e9 for x in v1))}]"
            print(f"  {f:9s} {fmt(v0):>18s} -> {fmt(v1):>18s}{extra}")
# pf projection shots at p_ref = 1 (chapter: 'at least 1.0-1.4e6 projection shots')
for a, lab in ((a0, 'now'), (a1, 'new')):
    lm = a.lam_over_m_pf.value; lm = lm if isinstance(lm, tuple) else (lm, lm)
    sh = tuple(r ** 2 / (4 * a.eps_simplified.value ** 2) for r in lm)
    print(f"pf projection shots ({lab}): {fmt(sh)}")
# model() 2033 era outputs that might carry these (intermediates mentioning rabi/f7)
for era in ('2033',):
    r0, r1 = m.model(a0, era), m.model(a1, era)
    for k in r0.intermediates:
        if r0.intermediates[k] != r1.intermediates[k]:
            print(f"model {era} intermediate {k}: {r0.intermediates[k]} -> {r1.intermediates[k]}")
    print(f"model {era} hard_ops {r0.hard_ops} -> {r1.hard_ops}; lq {r0.lq} -> {r1.lq}; wall {r0.wall_time_s} -> {r1.wall_time_s}")
