"""Discussion figure; no manuscript integration."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
d=json.loads((P/'triton_results.json').read_text());s=np.load(P/'triton_states.npz')
tau=np.linspace(0,.2,301);C=np.sum(s['krylov_weights'][:,None]*np.exp(-1j*(s['krylov_energies'][:,None]-d['ground_energy'])*tau),axis=0)
pole=d['elastic_weight']*np.exp(-1j*(d['moving_ground_energy']-d['ground_energy'])*tau)
rows=[r for r in d['rows'] if r['dt']==.002]
x=np.array([r['tau'] for r in rows]);v=np.array([complex(*r['approximate']) for r in rows])
fig,axes=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
for ax,y,z,title in ((axes[0],C,v,'Full charge correlator'),(axes[1],C-pole,v-d['elastic_weight']*np.exp(-1j*(d['moving_ground_energy']-d['ground_energy'])*x),'Elastic pole subtracted')):
    ax.plot(tau*197.32698,y.real,label='Real, reference',color='#176b99')
    ax.plot(tau*197.32698,y.imag,label='Imaginary, reference',color='#b45309')
    ax.scatter(x*197.32698,z.real,color='#176b99',marker='o',label='Budgeted circuit')
    ax.scatter(x*197.32698,z.imag,color='#b45309',marker='o')
    ax.axhline(0,color='0.75',lw=.6);ax.set(xlabel='Time (fm/c)',ylabel='Normalized correlator',title=title)
axes[0].legend(fontsize=8)
axes[2].plot(np.array([4,6,8])*1.4,[d['L4_energy_same_couplings'],d['ground_energy'],d['L8_energy_same_couplings']],'o-',color='#7251a2')
axes[2].annotate('Fitted here',(8.4,d['ground_energy']),xytext=(7,-12),arrowprops=dict(arrowstyle='->'),fontsize=9)
axes[2].set(xlabel='Box side (fm)',ylabel='Ground energy (MeV)',title='Same couplings: volume dependence')
fig.suptitle('Three mobile nucleons: finite-volume model, not a converged breakup cross section',fontsize=12)
for ext in ('png','pdf'):fig.savefig(P/f'triton_check.{ext}',dpi=180)
