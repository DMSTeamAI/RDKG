import pandas as pd, collections
import os; R=os.path.join(os.environ.get('RDKG_DIR','.'),'')
e=pd.read_csv(R+'edges_no_trials.csv',dtype=str); n=pd.read_csv(R+'nodes_no_trials.csv',dtype=str)
S,O,T=':START_ID',':END_ID',':TYPE'
pre=lambda s:s.str.split(':').str[0]
rel=e[e[T]=='RELATED_TO'].copy()
# disease-gene mirrored
dg=rel[(rel.subject_type=='biolink:Disease')&(rel.object_type=='biolink:Gene')]
gd=rel[(rel.subject_type=='biolink:Gene')&(rel.object_type=='biolink:Disease')]
a=set(zip(dg[S],dg[O])); b=set((x,y) for y,x in zip(gd[S],gd[O]))
print('DG',len(a),'GD',len(b),'overlap',len(a&b),'union unique pairs',len(a|b))
print('gene id prefixes',pre(dg[O]).value_counts().to_dict(), pre(gd[S]).value_counts().to_dict())
# disease-disease
dd=rel[(rel.subject_type=='biolink:Disease')&(rel.object_type=='biolink:Disease')]
print(pd.crosstab(pre(dd[S]),pre(dd[O])))
dd_pairs=set(zip(dd[S],dd[O])); print('dd mirrored',sum((y,x) in dd_pairs for x,y in dd_pairs))
# overlap with subclass edges
sc=e[e[T]=='SUBCLASS_OF']; scp=set(zip(sc[S],sc[O])); print('dd also subclass',sum((x,y) in scp or (y,x) in scp for x,y in dd_pairs))
# related_to disease->gene: also genetic_association/other
ga=e[e[T]=='GENETIC_ASSOCIATION']; print('GA pairs also in related',len(set(zip(ga[S],ga[O]))&a))
# orphanet nodes
nm=n.set_index('id:ID')['name']
orph=[i for i in n['id:ID'] if i.startswith('Orphanet')]
anyedge=set(e[S])|set(e[O])
conn=[o for o in orph if o in anyedge]
import random; random.seed(1)
print('isolated sample',[nm[o] for o in random.sample([o for o in orph if o not in anyedge],25)])
print('connected sample',[nm[o] for o in random.sample(conn,10)])
x=e[e[S].isin(conn)|e[O].isin(conn)]; print(x[T].value_counts().to_dict())
# keyword stratification of isolated names
iso=[nm[o] for o in orph if o not in anyedge]
import re
kw=collections.Counter()
for s in iso:
    s=str(s).lower()
    for k in ['syndrome','deficiency','type','susceptibility','isolated','familial','with','due to','subtype','non-syndromic','autosomal','x-linked','cancer','carcinoma','tumor','hereditary','inherited']:
        if k in s: kw[k]+=1
print(len(iso),kw.most_common())
# grouped nodes: how many members
g=[i for i in n['id:ID'] if i.startswith('MONDO_GROUPED')]
print('grouped',len(g),collections.Counter(len(i.split(':')[1].split('_'))-0 for i in g).most_common(8))
# node counts by prefix for Drug
print(pre(n[n.category=='biolink:Drug']['id:ID']).value_counts().to_dict())
# ClinVar variants: no disease link?
var=set(n[n.category=='biolink:SequenceVariant']['id:ID'])
cc=e[e[T]=='CAUSES_CONDITION']; print('variants with CAUSES',cc[S].nunique(),'targets prefix',pre(cc[O]).value_counts().to_dict())
cb=e[e[T]=='CONTRIBUTES_TO']; print('variants contributes',len(set(cb[S])&var))
print('variants with any disease link',len((set(cc[S])|set(cb[S]))&var),'of',len(var))
