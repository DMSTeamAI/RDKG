import zlib, json, random, math, collections, numpy as np
ann={k:set(v) for k,v in json.load(open('ann.json')).items()}
par={k:set(v) for k,v in json.load(open('par.json')).items()}
memo={}
def anc(t):
    if t in memo: return memo[t]
    s={t}
    for p in par.get(t,()): s|=anc(p)
    memo[t]=s; return s
ROOT='HP:0000001'; PHENO='HP:0000118'
def closure(ts):
    c=set()
    for t in ts: c|=anc(t)
    return c-{ROOT,PHENO}
dis=sorted(ann); idx={d:i for i,d in enumerate(dis)}; N=len(dis)
clo={d:closure(ann[d]) for d in dis}
df=collections.Counter(t for d in dis for t in clo[d])
IC={t:-math.log(c/N) for t,c in df.items()}
inv=collections.defaultdict(list)
for d in dis:
    for t in clo[d]: inv[t].append(idx[d])
inv={t:np.array(v) for t,v in inv.items()}
# exact index
inv_ex=collections.defaultdict(list)
for d in dis:
    for t in ann[d]: inv_ex[t].append(idx[d])
inv_ex={t:np.array(v) for t,v in inv_ex.items()}
def score(q, method, excl=None):
    s=np.zeros(N)
    for t in q:
        if method=='exact':
            if t in inv_ex: s[inv_ex[t]]+=1
        elif method=='closure_count':
            for a in anc(t)-{ROOT,PHENO}:
                pass
            if t in inv: s[inv[t]]+=1
        elif method=='ic_overlap':   # IC of query term if in disease closure
            if t in inv: s[inv[t]]+=IC.get(t,0)
        elif method=='resnik':       # best IC of common ancestor, averaged over query terms
            best=np.zeros(N)
            for a in anc(t)-{ROOT,PHENO}:
                if a in inv: 
                    ic=IC.get(a,0); v=inv[a]; best[v]=np.maximum(best[v],ic)
            s+=best/len(q)
    return s
def rank_of(s,i):
    return 1+(s>s[i]).sum()+0.5*((s==s[i]).sum()-1)
def sim(rng,d,k,imprec,noise,holdout):
    a=sorted(ann[d]); q=rng.sample(a,k)
    q2=[]
    for t in q:
        if rng.random()<imprec:
            ps=sorted(par.get(t,()))
            q2.append(rng.choice(ps) if ps and ps[0]!=ROOT else t)
        else: q2.append(t)
    allt=sorted(inv_ex)
    for _ in range(noise): q2.append(rng.choice(allt))
    return q,q2
rng=random.Random(42)
cand=[d for d in dis if len(ann[d])>=8]
sample=rng.sample(cand,1000)
methods=['exact','closure_count','ic_overlap','resnik']
scen={'A_exact_in_graph':dict(imprec=0,noise=0,holdout=False),
      'B_imprecise_noisy':dict(imprec=0.5,noise=2,holdout=False),
      'C_heldout_imprecise_noisy':dict(imprec=0.5,noise=2,holdout=True)}
res={}
for sname,p in scen.items():
    ranks={m:[] for m in methods}
    for d in sample:
        r=random.Random(zlib.crc32(f"{sname}|{d}".encode()))
        q,q2=sim(r,d,5,p['imprec'],p['noise'],p['holdout'])
        if p['holdout']:
            # remove sampled annotations (and queries' terms) from this disease profile
            saved=(ann[d],clo[d]); newann=ann[d]-set(q)
            # temporarily rebuild this disease's contribution
            old_ex=[(t,) for t in ann[d]]
            # modify indices
            i=idx[d]
            for t in ann[d]:
                inv_ex[t]=inv_ex[t][inv_ex[t]!=i]
            for t in clo[d]:
                inv[t]=inv[t][inv[t]!=i]
            ann[d]=newann; clo[d]=closure(newann)
            for t in ann[d]: inv_ex[t]=np.append(inv_ex.get(t,np.array([],dtype=int)),i)
            for t in clo[d]: inv[t]=np.append(inv.get(t,np.array([],dtype=int)),i)
        for m in methods:
            ranks[m].append(rank_of(score(q2,m),idx[d]))
        if p['holdout']:
            i=idx[d]
            for t in ann[d]: inv_ex[t]=inv_ex[t][inv_ex[t]!=i]
            for t in clo[d]: inv[t]=inv[t][inv[t]!=i]
            ann[d]=saved[0]; clo[d]=saved[1]
            for t in ann[d]: inv_ex[t]=np.append(inv_ex[t],i)
            for t in clo[d]: inv[t]=np.append(inv[t],i)
    res[sname]={}
    for m in methods:
        r=np.array(ranks[m])
        # bootstrap CI for MRR and top10
        b=np.random.default_rng(0)
        def ci(f):
            vals=[f(r[b.integers(0,len(r),len(r))]) for _ in range(500)]
            return float(np.percentile(vals,2.5)),float(np.percentile(vals,97.5))
        top1=float((r<=1.5).mean()); top10=float((r<=10).mean()); mrr=float((1/r).mean())
        res[sname][m]=dict(top1=top1,top10=top10,mrr=mrr,mrr_ci=ci(lambda x:(1/x).mean()),top10_ci=ci(lambda x:(x<=10).mean()),median=float(np.median(r)))
        print(sname,m,{k:(round(v,3) if not isinstance(v,tuple) else tuple(round(z,3) for z in v)) for k,v in res[sname][m].items()})
json.dump(res,open('bench.json','w'),indent=1)
print('pool',N,'sampled',len(sample))
