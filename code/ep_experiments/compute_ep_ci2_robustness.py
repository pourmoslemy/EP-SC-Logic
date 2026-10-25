#!/usr/bin/env python3
"""Exact CI2 cross-witness, Alexandrov group semantics, artifacts, and HNS audit."""
from __future__ import annotations
import json, math, re, time, zipfile
from collections import Counter
from itertools import combinations
from pathlib import Path
from xml.etree import ElementTree as ET
import networkx as nx
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'VoteWatch-EP-voting-data_2004-2022'
FILES={'EP6':SOURCE/'EP6_RCVs_2022_06_13.xlsx','EP9':SOURCE/'EP9_RCVs_2022_06_22.xlsx'}
TAUS=np.round(np.arange(.10,.901,.02),2); SEED=42; NPERM=1000
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
def log(x):print(f'[INFO] {x}',flush=True)
def col(ref):
 z=0
 for c in re.match(r'[A-Z]+',ref).group(0):z=z*26+ord(c)-64
 return z-1
def shared(z):
 root=ET.fromstring(z.read('xl/sharedStrings.xml'));return [''.join(t.text or '' for t in q.iter(NS+'t')) for q in root]
def parse_no(path,term,n_expected,v_expected):
 cache=HERE/f'ep_no_votes_{term}.npz'
 if cache.exists():return np.load(cache,allow_pickle=False)['no'][:,:v_expected]
 with zipfile.ZipFile(path) as z:
  ss=shared(z)
  def val(c):
   v=c.find(NS+'v');return '' if v is None or v.text is None else (ss[int(v.text)] if c.get('t')=='s' else v.text)
  header=2 if term=='EP6' else 1;start=4 if term=='EP6' else 2;vote_cols=[];rows=[]
  for _,e in ET.iterparse(z.open('xl/worksheets/sheet1.xml'),events=('end',)):
   if e.tag!=NS+'row':continue
   rn=int(e.get('r','0'));d={col(c.get('r')):val(c) for c in e.findall(NS+'c')}
   if rn==header:vote_cols=[j for j in sorted(d) if j>=10 and d[j] not in ('',None)]
   elif rn>=start and d.get(0,''):
    rows.append(np.fromiter((d.get(j)=='2' for j in vote_cols),dtype=bool,count=len(vote_cols)))
   e.clear()
 no=np.vstack(rows)[:,:v_expected]
 if no.shape[0]!=n_expected:raise ValueError(f'{term}: {no.shape[0]} rows, expected {n_expected}')
 np.savez_compressed(cache,no=no);return no
def rank3(i,j,k,n):
 a=j-i-1
 return (math.comb(n,3)-math.comb(n-i,3))+a*(n-i-2)-a*(a-1)//2+(k-j-1)
def mark_witnesses(no):
 n,r=no.shape;marked=np.zeros(math.comb(n,3),dtype=np.bool_)
 for vote in range(r):
  verts=np.flatnonzero(no[:,vote]);q=len(verts)
  if q<3:continue
  for ai in range(q-2):
   i=int(verts[ai])
   for bj in range(ai+1,q-1):
    j=int(verts[bj]);ks=verts[bj+1:];a=j-i-1
    base=(math.comb(n,3)-math.comb(n-i,3))+a*(n-i-2)-a*(a-1)//2-j-1
    marked[base+ks]=True
  if vote%1000==0:log(f'witness votes {vote}/{r}; marked={int(marked.sum())}')
 return marked
def triangle_hist(weights,marked):
 n=len(weights);total=np.zeros(len(TAUS),np.int64);wit=np.zeros(len(TAUS),np.int64)
 adj=weights>=TAUS[0]
 for i in range(n-2):
  js=np.flatnonzero(adj[i,i+1:])+i+1
  for j in js:
   ks=np.flatnonzero(adj[i,j+1:]&adj[j,j+1:])+j+1
   if not len(ks):continue
   birth=np.minimum(np.minimum(weights[i,j],weights[i,ks]),weights[j,ks])
   idx=np.searchsorted(TAUS,birth,side='right')-1
   ranks=np.fromiter((rank3(i,int(j),int(k),n) for k in ks),dtype=np.int64,count=len(ks));flags=marked[ranks]
   total+=np.bincount(idx,minlength=len(TAUS));wit+=np.bincount(idx[flags],minlength=len(TAUS))
 # histograms are births; cumulative from high threshold downward
 return np.cumsum(total[::-1])[::-1],np.cumsum(wit[::-1])[::-1]
def graph(weights,tau):
 n=len(weights);g=nx.Graph();g.add_nodes_from(range(n));u,v=np.nonzero(np.triu(weights>=tau,1));g.add_edges_from(zip(u.tolist(),v.tolist()));return g
def choose3(x):return x*(x-1)*(x-2)//6 if x>=3 else 0
def heyting(g,groups,target):
 G={i for i,x in enumerate(groups) if x==target};N={v for x in G for v in g.neighbors(x) if v not in G};neg=set(g)-G-N
 Nneg={v for x in neg for v in g.neighbors(x) if v not in neg};dneg=(set(g)-neg)-Nneg;defect=len(dneg-G)
 boundary=N;U=G|boundary
 # Inclusion-exclusion: triangles meeting both disjoint sets G and boundary.
 tri_u=sum(nx.triangles(g.subgraph(U)).values())//3
 tri_g=sum(nx.triangles(g.subgraph(G)).values())//3
 tri_b=sum(nx.triangles(g.subgraph(boundary)).values())//3
 obs=tri_u-tri_g-tri_b
 potential=choose3(len(U))-choose3(len(G))-choose3(len(boundary))
 return {'neg':len(neg),'double_neg':len(dneg),'defect':defect,'boundary':len(boundary),'CI2_H':obs/potential if potential else float('nan'),'triangles_boundary_observed':obs,'triangles_boundary_possible':potential}
def hns(g):
 deg=dict(g.degree());active=[v for v,d in deg.items() if d>0];leaves={v for v,d in deg.items() if d==1};defective={next(iter(g.neighbors(v))) for v in leaves};return len(defective&set(active))/len(active) if active else 0.
def main():
 base=pd.read_csv(HERE/'ep_filtration_trajectories_v2.csv');diag_path=HERE/'ep_hns_diagnostics.json';diag=json.loads(diag_path.read_text(encoding='utf-8'));rob={'seed':SEED,'permutations':NPERM,'ci2_h_definition':'observed triangles intersecting both G and N(G), divided by all possible triples in G union N(G) intersecting both sets','terms':{}}
 out=base.copy();summary=[]
 for term in ('EP6','EP9'):
  z=np.load(HERE/f'ep_hns_weights_{term}.npz',allow_pickle=False);weights=z['weights'];groups=z['groups'].astype(str).tolist();md=json.loads(str(z['metadata']));no=parse_no(FILES[term],term,len(groups),int(md['roll_calls']));log(f'{term}: no matrix {no.shape}')
  wcache=HERE/f'ep_no_triple_witness_{term}.npz'
  if wcache.exists():marked=np.load(wcache,allow_pickle=False)['marked']
  else:marked=mark_witnesses(no);np.savez_compressed(wcache,marked=marked)
  tot,wit=triangle_hist(weights,marked);ci=np.ones(len(TAUS),float);np.divide(wit,tot,out=ci,where=tot>0);ci=np.where(tot>0,1-ci,1.0)
  out[f'CI2_cross_vote_{term}']=ci;out[f'flag_triangles_{term}']=tot;out[f'no_witness_triangles_{term}']=wit
  records=[]
  for q,tau in enumerate(TAUS):
   g=graph(weights,float(tau));active=sum(d>0 for _,d in g.degree());edges=g.number_of_edges();components=nx.number_connected_components(g);artifact=edges<=5 or active<=10 or components==edges
   out.loc[q,f'artifact_{term}']=artifact;out.loc[q,f'edge_count_{term}']=edges;out.loc[q,f'active_vertices_{term}']=active;out.loc[q,f'component_count_{term}']=components
   rec={'tau':float(tau),'CI2_cross_vote':float(ci[q]),'flag_triangles':int(tot[q]),'shared_no_witness_triangles':int(wit[q]),'edge_count':edges,'active_vertices':active,'component_count':components,'artifact':bool(artifact),'groups':{}}
   for target,label in [('EPP','EPP'),('Socialists','S&D'),('Liberals','Renew')]:
    x=heyting(g,groups,label);rec['groups'][target]=x
    for key in ('neg','double_neg','defect','boundary','CI2_H'):out.loc[q,f'{key}_{target}_{term}']=x[key]
   records.append(rec)
  hs=out[f'HNS_{term}'];mask=~out[f'artifact_{term}'].astype(bool);peak_i=int(hs.idxmax());peak_tau=float(out.loc[peak_i,'tau']);gpeak=graph(weights,peak_tau);obs=hns(gpeak)
  # HNS is topology-only. Label permutations are still executed/audited and are necessarily identical.
  rng=np.random.default_rng(SEED);null=np.empty(NPERM)
  ga=np.asarray(groups,dtype=object)
  for b in range(NPERM):rng.permutation(ga);null[b]=hns(gpeak)
  sd=float(null.std(ddof=1));zscore=(obs-float(null.mean()))/sd if sd>1e-12 else 0.0;p=float(np.mean(null>=obs-1e-15))
  test={'tau_peak':peak_tau,'HNS_obs':obs,'null_mean':float(null.mean()),'null_std':sd,'z_score':zscore,'empirical_p':p,'null_distribution':null.tolist(),'identifiability_note':'HNS is topology-only; shuffling group labels leaves it invariant, so this requested null is degenerate.'}
  row64=records[int(np.where(TAUS==.64)[0][0])];stats={'with_artifacts':{'mean':float(hs.mean()),'median':float(hs.median()),'peak':float(hs.max())},'without_artifacts':{'mean':float(hs[mask].mean()),'median':float(hs[mask].median()),'peak':float(hs[mask].max())},'tau_0_64':row64}
  rob['terms'][term]={'thresholds':records,'artifact_diagnostics':stats,'hns_peak_permutation':test}
  summary.append({'Parliament':term,'tau':peak_tau,'HNS_obs':obs,'p-value':p,'z-score':zscore,'CI2_cross':float(ci[peak_i]),'artifact':bool(out.loc[peak_i,f'artifact_{term}'])})
 out.to_csv(HERE/'ep_filtration_trajectories_v3.csv',index=False,float_format='%.12g');diag['robustness_tests']=rob;diag_path.write_text(json.dumps(diag,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')
 print("| Parliament | tau | HNS_obs | p-value | z-score | CI2_cross | artifact |\n|---|---:|---:|---:|---:|---:|:---:|");print("\n".join(f"| {r['Parliament']} | {r['tau']:.2f} | {r['HNS_obs']:.6f} | {r['p-value']:.6f} | {r['z-score']:.6f} | {r['CI2_cross']:.6f} | {r['artifact']} |" for r in summary));log("wrote ep_filtration_trajectories_v3.csv and appended diagnostics")
if __name__=='__main__':main()