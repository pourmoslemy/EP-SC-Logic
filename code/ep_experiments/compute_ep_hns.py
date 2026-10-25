#!/usr/bin/env python3
"""Compute EP Heyting-negation and cross-group metrics (stdlib/numpy/pandas/networkx)."""
from __future__ import annotations
import json, math, random, re, time, zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
import networkx as nx
import numpy as np
import pandas as pd

SEED=42; TAUS=np.round(np.arange(.10,.901,.02),2); CAP=100_000
HERE=Path(__file__).resolve().parent; SOURCE=HERE.parent/'VoteWatch-EP-voting-data_2004-2022'
FILES={'EP6':SOURCE/'EP6_RCVs_2022_06_13.xlsx','EP9':SOURCE/'EP9_RCVs_2022_06_22.xlsx'}
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
def log(s): print(f'[INFO] {s}',flush=True)
def col_index(ref):
    z=0
    for c in re.match(r'[A-Z]+',ref).group(0): z=z*26+ord(c)-64
    return z-1
def strings(z):
    if 'xl/sharedStrings.xml' not in z.namelist(): return []
    root=ET.fromstring(z.read('xl/sharedStrings.xml'))
    return [''.join(t.text or '' for t in x.iter(NS+'t')) for x in root]
def value(c,ss):
    if c.get('t')=='inlineStr': return ''.join(t.text or '' for t in c.iter(NS+'t'))
    v=c.find(NS+'v')
    if v is None or v.text is None:return ''
    return ss[int(v.text)] if c.get('t')=='s' else v.text
def group(raw):
    s=re.sub(r'\s+',' ',str(raw or '').strip()).upper()
    rules=[('The Left',('THE LEFT','GUE','NGL','UNITED LEFT')),('Greens/EFA',('GREEN','FREE ALLIANCE')),('EPP',('EPP','EUROPEAN PEOPLE')),('S&D',('S&D','SOCIALIST','PES')),('Renew',('RENEW','REG','ALDE','LIBERAL')),('ECR',('ECR','EUROPEAN CONSERVATIVES','UNION FOR EUROPE OF THE NATIONS')),('ID',('IDG','IDENTITY AND DEMOCRACY','INDEPENDENCE/DEMOCRACY')),('NI',('NON-ATTACHED','NON ATTACHED','NI','NA'))]
    for name,needles in rules:
        if any(x==s or x in s for x in needles): return name
    return 'NI'
def load_data(path,term):
    cache=HERE/f'ep_hns_weights_{term}.npz'
    if cache.exists():
        d=np.load(cache,allow_pickle=False);g=d['groups'].astype(str).tolist();m=json.loads(str(d['metadata']));m['group_counts']=dict(Counter(g));return d['weights'],g,m
    started=time.perf_counter(); rows=[]; nvotes=0; vote_cols=[]
    with zipfile.ZipFile(path) as z:
        ss=strings(z); header=2 if term=='EP6' else 1; start=4 if term=='EP6' else 2
        with z.open('xl/worksheets/sheet1.xml') as stream:
            for _,el in ET.iterparse(stream,events=('end',)):
                if el.tag!=NS+'row':continue
                rn=int(el.get('r','0')); vals={col_index(c.get('r')):value(c,ss) for c in el.findall(NS+'c')}
                if rn==header:
                    vote_cols=[j for j in sorted(vals) if j>=10 and vals[j] not in ('',None)];nvotes=len(vote_cols)
                if rn>=start and vals.get(0,'') and nvotes:
                    try: mid=int(float(vals[0]))
                    except ValueError: el.clear(); continue
                    yes=np.zeros(nvotes,dtype=np.uint8)
                    for j,x in vals.items():
                        if j>=10 and j-10<nvotes and x=='1':yes[j-10]=1
                    rows.append((mid,group(vals.get(7,'')),yes))
                el.clear()
    if not rows:raise RuntimeError(f'No vote rows parsed from {path}')
    matrix=np.vstack([r[2] for r in rows]).astype(np.uint16)
    weights=(matrix@matrix.T).astype(float)/matrix.shape[1];np.fill_diagonal(weights,0)
    groups=[r[1] for r in rows]
    meta={'term':term,'source':path.name,'members':len(rows),'roll_calls':int(matrix.shape[1]),'group_counts':dict(Counter(groups)),'coverage_note':'EP9 source ends in June 2022; it does not cover the full 2019-2024 term.' if term=='EP9' else ''}
    np.savez_compressed(cache,weights=weights,groups=np.asarray(groups),metadata=json.dumps(meta,sort_keys=True))
    log(f'{term}: parsed {len(rows)} MEP rows and {matrix.shape[1]} votes in {time.perf_counter()-started:.1f}s')
    return weights,groups,meta
def direct_defect(g,v):
    faces=[('v',u) for u in g]+[('e',min(u,w),max(u,w)) for u,w in g.edges]
    V={('v',v)}|{('e',min(v,u),max(v,u)) for u in g.neighbors(v)}
    def cof(f):
        if f[0]=='e':return {f}
        u=f[1];return {f}|{('e',min(u,w),max(u,w)) for w in g.neighbors(u)}
    neg={f for f in faces if cof(f).isdisjoint(V)}
    return {f for f in faces if cof(f).isdisjoint(neg)}!=V
def metrics(weights,groups,tau):
    n=len(groups); pairs=np.transpose(np.nonzero(np.triu(weights>=tau,1)))
    g=nx.Graph();g.add_nodes_from(range(n));g.add_edges_from((int(u),int(v)) for u,v in pairs)
    deg=np.fromiter((d for _,d in g.degree()),dtype=np.int64,count=n);active=deg>0;leaves=deg==1;defective=np.zeros(n,dtype=bool)
    for leaf in np.flatnonzero(leaves):defective[next(iter(g.neighbors(int(leaf))))]=True
    ac=int(active.sum());dc=int((defective&active).sum());hns=dc/ac if ac else 0.
    sizes=[];cross=0;cap=False
    for i,c in enumerate(nx.find_cliques(g)):
        if i>=CAP:cap=True;break
        sizes.append(len(c));cross+=len({groups[v] for v in c})>=2
    nc=len(sizes);beta=cross/nc if nc else 0.;counts=Counter(sizes)
    if len(counts)>1:
        p=np.asarray(list(counts.values()),float)/nc;entropy=float(-np.sum(p*np.log(p))/math.log(len(counts)))
    else:entropy=0.
    m=g.number_of_edges();cross_edges=sum(groups[u]!=groups[v] for u,v in g.edges);gc=Counter(groups);allpairs=n*(n-1)//2;same=sum(x*(x-1)//2 for x in gc.values());possible=allpairs-same
    density=m/allpairs if allpairs else 0.;cd=(cross_edges/possible)/(density+1e-9) if density and possible else 0.
    sample=sorted(set(np.flatnonzero(active).tolist()[:3]+np.flatnonzero(defective).tolist()[:3]))
    audits=[{'vertex_index':int(v),'direct_defect':direct_defect(g,int(v)),'closed_form_defect':bool(defective[v]),'match':direct_defect(g,int(v))==bool(defective[v])} for v in sample]
    out={'HNS':hns,'beta_cross':beta,'bool_core':1-hns,'coal_entropy':entropy,'cross_density':cd}
    audit={'tau':tau,'nodes':n,'edges':m,'active_vertices':ac,'isolated_vertices':int((deg==0).sum()),'leaf_vertices':int(leaves.sum()),'defective_active_vertices':dc,'maximal_cliques_enumerated':nc,'cross_group_maximal_cliques':cross,'clique_cap':CAP,'clique_cap_hit':cap,'poset_identity_audit':audits,**out}
    return out,audit
def summary(results,meta):
    lines=['# European Parliament Heyting Negation Spectrum summary','', 'The calculation uses the 1-skeleton face poset with Alexandrov upsets. Here, `HNS > 0` is exact evidence of non-Boolean behavior **within this specified finite Heyting model**: for at least one active vertex, double negation of its principal upset differs from that upset. It is not, by itself, a model-independent causal claim about political intent.','', '| Term | tau | HNS | Boolean core | Cross-group maximal-clique ratio | Clique-size entropy | Relative cross-group density |','|---|---:|---:|---:|---:|---:|---:|']
    for term in ('EP6','EP9'):
        for tau in (.30,.50,.70):
            r=results[term][tau];lines.append(f"| {term} | {tau:.2f} | {r['HNS']:.6f} | {r['bool_core']:.6f} | {r['beta_cross']:.6f} | {r['coal_entropy']:.6f} | {r['cross_density']:.6f} |")
    lines+=['','## Interpretation and data scope','', 'In a graph regarded as a one-dimensional simplicial complex, double negation adds precisely those degree-one neighbors attached to the focal vertex. Consequently, HNS is the share of active vertices adjacent to at least one leaf. The reported cross-group quantities are combinatorial summaries of thresholded co-Yes-voting graphs and should not be read as latent preferences or causal effects.','',f"EP6 uses {meta['EP6']['members']} membership rows and {meta['EP6']['roll_calls']} roll calls. EP9 uses {meta['EP9']['members']} membership rows and {meta['EP9']['roll_calls']} roll calls. The supplied EP9 VoteWatch source ends in June 2022 and therefore does not cover the full 2019-2024 term.",'','Maximal-clique results are flagged in the diagnostics whenever the 100,000-clique cap is reached; capped values describe the deterministic enumeration prefix rather than the full clique family.']
    (HERE/'ep_hns_summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
def main():
    random.seed(SEED);np.random.seed(SEED);started=time.perf_counter();results={x:{} for x in ('EP6','EP9')};diag={'seed':SEED,'thresholds':TAUS.tolist(),'terms':{},'method':{'clique_cap':CAP}};meta={}
    for term in ('EP6','EP9'):
        w,groups,md=load_data(FILES[term],term);meta[term]=md;diag['terms'][term]={'metadata':md,'thresholds':[]}
        for step,tv in enumerate(TAUS,1):
            tau=float(tv);tick=time.perf_counter();r,a=metrics(w,groups,tau);results[term][tau]=r;diag['terms'][term]['thresholds'].append(a)
            log(f"{term} tau={tau:.2f} ({step:02d}/41): edges={a['edges']}, HNS={r['HNS']:.6f}, cliques={a['maximal_cliques_enumerated']}{' [CAP]' if a['clique_cap_hit'] else ''} in {time.perf_counter()-tick:.2f}s")
    base=pd.read_csv(HERE/'ep_filtration_trajectories.csv');wide=None
    for term in ('EP6','EP9'):
        f=base.loc[base.term==term].drop(columns='term').copy();f=f.rename(columns={c:f'{c}_{term}' for c in f if c!='tau'});wide=f if wide is None else wide.merge(f,on='tau',how='outer',validate='one_to_one')
    for name in ('HNS','beta_cross','bool_core','coal_entropy','cross_density'):
        for term in ('EP6','EP9'):wide[f'{name}_{term}']=wide.tau.map(lambda x,t=term,m=name:results[t][round(float(x),2)][m])
    requested=[f'{n}_{t}' for n in ('HNS','beta_cross','bool_core','coal_entropy','cross_density') for t in ('EP6','EP9')];wide=wide[[c for c in wide if c not in requested]+requested].sort_values('tau').reset_index(drop=True);wide.to_csv(HERE/'ep_filtration_trajectories_v2.csv',index=False,float_format='%.12g')
    bounded=[c for c in requested if not c.startswith('cross_density_')];diag['validation']={'trajectory_rows':len(wide),'all_requested_columns_present':all(c in wide for c in requested),'all_metrics_finite':bool(np.isfinite(wide[requested].to_numpy()).all()),'unit_interval_metrics_valid':bool(((wide[bounded]>=0)&(wide[bounded]<=1)).all().all()),'all_poset_audits_match':all(a['match'] for t in diag['terms'].values() for x in t['thresholds'] for a in x['poset_identity_audit'])}
    (HERE/'ep_hns_diagnostics.json').write_text(json.dumps(diag,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8');summary(results,meta);log(f"Completed in {time.perf_counter()-started:.1f}s; validation={diag['validation']}")
if __name__=='__main__':main()
