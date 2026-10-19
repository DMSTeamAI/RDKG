exec(open('a10.py').read().split("if sys.argv[1:]==['eval']")[0])
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
M='MONDO:0009265'; sc=scores(M); P=ann[M]; icP=ic(P); known={g for g in drugs if M in treated[g]}
order=sorted(drugs,key=lambda g:-sc[g][1])[:12]
alld=sorted(ann); prng=random.Random(3); rows=[]
for g in order:
    ks=len(treated[g]); nl=[]
    for _ in range(200):
        ds=prng.sample(alld,ks); U=set().union(*[ann[x] for x in ds]); inter=P&U; nl.append(ic(inter)/(icP+ic(U)-ic(inter)))
    nl=np.array(nl); rows.append((nm[g],sc[g][1],(sc[g][1]-nl.mean())/(nl.std()+1e-9),g in known))
ert={'Imiglucerase','Velaglucerase alfa'}
col=lambda n,k:'#2e8b57' if k else ('#e08a1e' if n in ert else '#7a4fa3')
fig,ax=plt.subplots(figsize=(7.2,3.1),dpi=300)
y=range(len(rows))[::-1]
ax.barh(list(y),[r[1] for r in rows],color=[col(r[0],r[3]) for r in rows])
ax.set_yticks(list(y));ax.set_yticklabels([r[0] for r in rows],fontsize=7)
for yi,r in zip(y,rows): ax.text(r[1]+0.002,yi,f'z={r[2]:.1f}',va='center',fontsize=6.5)
ax.set_xlabel('IC-weighted Jaccard score (Gaucher disease type I, MONDO:0009265)',fontsize=7);ax.tick_params(axis='x',labelsize=7)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='#2e8b57',label='Known type I treatment'),Patch(color='#e08a1e',label='Approved Gaucher drug, not linked to type I'),Patch(color='#7a4fa3',label='Candidate (phenotype-matched)')],fontsize=6.5,loc='lower right',frameon=False)
for s in('top','right'): ax.spines[s].set_visible(False)
ax.set_xlim(0,max(r[1] for r in rows)*1.18);plt.tight_layout();plt.savefig('fig3_gaucher.png');print(rows)
