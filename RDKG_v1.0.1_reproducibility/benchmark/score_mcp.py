"""Score MCP-benchmark answers. Usage: python score_mcp.py questions.csv answers.csv nodes.csv
answers.csv columns: id, condition (no_tools | rdkg_mcp), answer (free text containing identifiers).
Reports per condition and category: mean precision, recall, F1 of identifiers vs gold, and the
fabricated-identifier rate (identifiers of the expected type that do not exist in RDKG nodes)."""
import sys, re, json, pandas as pd
q=pd.read_csv(sys.argv[1]); a=pd.read_csv(sys.argv[2]); nodes=set(pd.read_csv(sys.argv[3],dtype=str)['id:ID'])
pat={'HP':r'HP:\d{7}','MONDO':r'MONDO:\d{7}','NCBIGene':r'(?:NCBIGene:)?\b\d{2,9}\b','Drug':r'(?:DrugBank:DB\d{5}|CTDRUG:[A-Za-z0-9_\-]+)'}
rows=[]
for _,r in a.merge(q,on='id').iterrows():
    gold=set(json.loads(r.gold)); t=r.id_type
    found=set(re.findall(pat[t],str(r.answer)))
    if t=='NCBIGene': found={f if f.startswith('NCBIGene:') else 'NCBIGene:'+f for f in found}
    tp=len(found&gold); p=tp/len(found) if found else 0; rec=tp/len(gold)
    fab=sum(f not in nodes for f in found)/len(found) if found else 0
    rows.append(dict(id=r.id,condition=r.condition,category=r.category,precision=p,recall=rec,f1=2*p*rec/(p+rec) if p+rec else 0,fabricated=fab))
d=pd.DataFrame(rows); print(d.groupby('condition')[['precision','recall','f1','fabricated']].mean().round(3)); print(d.groupby(['condition','category'])[['precision','recall','f1','fabricated']].mean().round(3))
