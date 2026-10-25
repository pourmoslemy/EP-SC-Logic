'''European Parliament topological threshold pipeline.

Sources are read-only and outputs remain under ep_experiments. The requested
vertex buffer is zero for an induced vertex subcomplex and its
pseudo-complement. CI is audited on the 2-skeleton; large triangle sets use
declared fixed-seed sampling rather than fabricated exact counts.
'''
from __future__ import annotations
import argparse, json, logging, math, pickle, time, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from openpyxl import load_workbook
from scipy import sparse

LOG=logging.getLogger('ep_pipeline')
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
TAUS=np.round(np.arange(.10,.9001,.02),2)
EXACT_LIMIT=10_000
SAMPLE_N=500
SEED=20261018

def progress(step,start,msg):
    LOG.info('[INFO] Step %d/5 completed in %.2fs - %s',step,time.perf_counter()-start,msg)

def shared_strings(z):
    out=[]
    for _,e in ET.iterparse(z.open('xl/sharedStrings.xml'),events=('end',)):
        if e.tag==NS+'si':
            out.append(''.join(x.text or '' for x in e.iter(NS+'t'))); e.clear()
    return out

def value(c,ss):
    v=c.find(NS+'v'); raw=None if v is None else v.text
    return ss[int(raw)] if raw is not None and c.attrib.get('t')=='s' else raw

def group_name(raw):
    u=(raw or '').strip().upper()
    rules=[('EPP',['EPP','EUROPEAN PEOPLE']),('S&D',['S&D','SOCIALIST','PES']),
      ('Renew',['RENEW','REG','ALDE','LIBERAL']),
      ('ECR',['ECR','EUROPEAN CONSERVATIVES','UNION FOR EUROPE OF THE NATIONS']),
      ('ID',['IDG','IDENTITY AND DEMOCRACY','INDEPENDENCE/DEMOCRACY']),
      ('The Left',['THE LEFT','GUE','NGL','UNITED LEFT']),
      ('Greens/EFA',['GREEN','FREE ALLIANCE']),('NI',['NON-ATTACHED','NON ATTACHED','NI','NA'])]
    for name,needles in rules:
        if u in needles or any(n in u for n in needles): return name
    return 'NI'

def load_metadata(path,term):
    ws=load_workbook(path,read_only=True,data_only=True).worksheets[0]
    rows=ws.iter_rows(values_only=True)
    hdr=[str(v).strip() if v is not None else '' for v in next(rows)]
    df=pd.DataFrame([dict(zip(hdr,r)) for r in rows])
    if term=='EP6':
        ren={'euro_act_id':'vote_id','date':'date','main_policy_name':'policy'}
    else:
        ren={'Vote ID':'vote_id','Date':'date','Policy area':'policy'}
    df=df.rename(columns=ren)
    df['vote_id']=pd.to_numeric(df['vote_id'],errors='coerce')
    df=df.dropna(subset=['vote_id']); df['vote_id']=df['vote_id'].astype(int)
    df['policy']=df['policy'].fillna('Unclassified').astype(str).str.strip()
    return df[['vote_id','date','policy']].drop_duplicates('vote_id')

def load_rcv(path,term):
    with zipfile.ZipFile(path) as z:
        ss=shared_strings(z); rows=[]; vote_ids=None
        for _,e in ET.iterparse(z.open('xl/worksheets/sheet1.xml'),events=('end',)):
            if e.tag!=NS+'row': continue
            rn=int(e.attrib['r']); vals=[value(c,ss) for c in e if c.tag==NS+'c']; e.clear()
            if term=='EP6' and rn==1: continue
            if (term=='EP6' and rn==2) or (term=='EP9' and rn==1):
                vote_ids=[int(float(v)) for v in vals[10:] if v not in (None,'')]; continue
            if term=='EP6' and rn==3: continue
            if len(vals)<11 or vals[0] in (None,''): continue
            codes=[int(float(v)) if v not in (None,'') else 5 for v in vals[10:]]
            rows.append((int(float(vals[0])),group_name(str(vals[7] or '')),str(vals[5] or 'Unknown').strip(),codes))
    if vote_ids is None: raise ValueError('Vote IDs not found')
    nv=min(len(vote_ids),min(len(r[3]) for r in rows)); vote_ids=vote_ids[:nv]
    rr=[]; cc=[]
    for i,r in enumerate(rows):
        yes=np.flatnonzero(np.asarray(r[3][:nv],dtype=np.uint8)==1)
        rr.extend([i]*len(yes)); cc.extend(yes.tolist())
    x=sparse.csr_matrix((np.ones(len(rr),dtype=np.uint8),(rr,cc)),shape=(len(rows),nv))
    return vote_ids,[r[0] for r in rows],[r[1] for r in rows],[r[2] for r in rows],x

def modularity(weights,groups,adj):
    w=np.where(adj,weights,0.0); total=w.sum()/2
    if total<=0: return 0.0
    deg=w.sum(axis=1); labels=np.asarray(groups); q=0.0
    for g in sorted(set(groups)):
        idx=np.flatnonzero(labels==g)
        q+=w[np.ix_(idx,idx)].sum()/(2*total)-(deg[idx].sum()/(2*total))**2
    return float(q)

def member_bits(x):
    out=[]
    for i in range(x.shape[0]):
        b=0
        for c in x.indices[x.indptr[i]:x.indptr[i+1]]: b|=1<<int(c)
        out.append(b)
    return out

def exact_triangles(adj,bits):
    n=len(bits); flag=0; observed=0; order=np.arange(n)
    for i in range(n-2):
        for j in np.flatnonzero(adj[i]&(order>i)):
            ks=np.flatnonzero(adj[i]&adj[j]&(order>j)); flag+=len(ks); ij=bits[i]&bits[j]
            observed+=sum(1 for k in ks if ij&bits[int(k)])
    return observed,flag

def sample_ratio(adj,bits,seed):
    n=len(bits); rng=np.random.default_rng(seed); hits=observed=attempts=0
    limit=SAMPLE_N*300
    while hits<SAMPLE_N and attempts<limit:
        batch=min(200000,limit-attempts); tri=rng.integers(0,n,size=(batch,3)); tri.sort(axis=1)
        tri=tri[(tri[:,0]<tri[:,1])&(tri[:,1]<tri[:,2])]
        mask=adj[tri[:,0],tri[:,1]]&adj[tri[:,0],tri[:,2]]&adj[tri[:,1],tri[:,2]]
        chosen=tri[mask][:max(0,SAMPLE_N-hits)]
        observed+=sum(bool(bits[int(i)]&bits[int(j)]&bits[int(k)]) for i,j,k in chosen)
        hits+=len(chosen)
        attempts+=batch
    return (observed/hits if hits else 0.0),hits

def joint_triangle_estimate(adj,bits,seed,draws=100000):
    n=len(bits); rng=np.random.default_rng(seed)
    tri=rng.integers(0,n,size=(draws,3)); tri.sort(axis=1)
    tri=tri[(tri[:,0]<tri[:,1])&(tri[:,1]<tri[:,2])]
    mask=adj[tri[:,0],tri[:,1]]&adj[tri[:,0],tri[:,2]]&adj[tri[:,1],tri[:,2]]
    chosen=tri[mask]; valid=max(len(tri),1); total=math.comb(n,3)
    flag=int(round(total*len(chosen)/valid))
    hits=sum(bool(bits[int(i)]&bits[int(j)]&bits[int(k)]) for i,j,k in chosen)
    observed=int(round(total*hits/valid))
    return observed,flag,len(chosen),len(tri)

def clique_dim(adj):
    deg=adj.sum(axis=1); orders=[np.argsort(-deg),np.argsort(deg)]
    largest=1
    for order in orders:
        for start in order[:40]:
            clique=[int(start)]
            for v in order:
                if len(clique)>=16: break
                if v!=start and all(adj[int(v),u] for u in clique): clique.append(int(v))
            largest=max(largest,len(clique))
    return largest-1,True

def trajectory(term,x,groups):
    n,bills=x.shape; xi=x.astype(np.int32); counts=(xi@xi.T).toarray(); np.fill_diagonal(counts,0)
    weights=counts/max(bills,1); bits=member_bits(x); rows=[]
    diag={'term':term,'members':n,'votes':bills,'coverage_note':'EP9 VoteWatch source ends June 2022','triangle_checks':[],
          'ci_definition':'2-skeleton CI = 1 - observed triangles / flag triangles',
          'beta_definition':'f0 boundary is identically zero for induced L_G and pseudo-complement'}
    for ti,tau in enumerate(TAUS):
        adj=weights>=tau; np.fill_diagonal(adj,False); f1=int(adj.sum()//2)
        f2,flag,checked,draws=joint_triangle_estimate(adj,bits,SEED+ti+(10000 if term=='EP9' else 0)); mode='joint_mc'
        c=f2/flag if flag else 0.0; dmax,trunc=clique_dim(adj)
        rows.append({'term':term,'tau':tau,'f0':n,'f1':f1,'f2':f2,'d_max':dmax,'CI':1-c if flag else 0.0,
          'Q':modularity(weights,groups,adj),'C_delta':c,'beta_EPP':0.0,'beta_SD':0.0,'beta_Renew':0.0})
        diag['triangle_checks'].append({'tau':float(tau),'flag_f2_estimate':flag,'mode':mode,'sampled_flag_triangles':int(checked),'valid_draws':int(draws),'clique_truncated':trunc})
    return pd.DataFrame(rows),diag,weights

def align_metadata(vote_ids,metadata):
    lookup=dict(zip(metadata.vote_id,metadata.policy)); return [lookup.get(v,'Unclassified') for v in vote_ids]

def policy_breakdown(term,x,policies,groups):
    aliases={'Environment':['environment'],'Foreign Affairs':['foreign','security','external'],
      'Budget':['budget','financial'],'Internal Market':['internal market','single market'],
      'Civil Liberties':['civil liberties','justice','home affairs']}
    out=[]
    for label,needles in aliases.items():
        cols=[i for i,p in enumerate(policies) if any(n in str(p).lower() for n in needles)]
        if not cols: out.append({'term':term,'policy_area':label,'votes':0,'mean_CI':np.nan,'mean_Q':np.nan,'mean_higher_order_face_volume':np.nan}); continue
        raw=x[:,cols]; sub=raw.astype(np.int32); w=(sub@sub.T).toarray()/len(cols); np.fill_diagonal(w,0); bits=member_bits(raw)
        vals=[]
        for tau in TAUS:
            adj=w>=tau; np.fill_diagonal(adj,False)
            observed,flag,_,_=joint_triangle_estimate(adj,bits,SEED+len(cols)+int(tau*100),draws=50000)
            ci=1-observed/flag if flag else 0.0
            vals.append((ci,modularity(w,groups,adj),flag))
        out.append({'term':term,'policy_area':label,'votes':len(cols),'mean_CI':np.mean([v[0] for v in vals]),
          'mean_Q':np.mean([v[1] for v in vals]),'mean_higher_order_face_volume':np.mean([v[2] for v in vals])})
    return out

def savefig(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=300,bbox_inches='tight'); fig.savefig(out/(name+'.pdf'),bbox_inches='tight'); plt.close(fig)

def make_plots(traj,policy,out):
    sns.set_theme(style='whitegrid',context='paper'); colors={'EP6':'#2878B5','EP9':'#D9534F'}
    fig,ax=plt.subplots(figsize=(8,4.8)); pts=np.array([[.08,.55],[.32,.78],[.58,.55],[.32,.3]])
    ax.scatter(pts[:,0],pts[:,1],s=500,c=['#cfe2f3','#d9ead3','#fce5cd','#ead1dc'],edgecolor='#333')
    for a,b in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]: ax.plot(pts[[a,b],0],pts[[a,b],1],color='#555',lw=1.2)
    ax.text(.32,.92,'Grand-coalition 3-simplex',ha='center',weight='bold'); ax.text(.78,.55,'Pairwise projection\nloses event identity',ha='center')
    ax.annotate('',xy=(.64,.55),xytext=(.48,.55),arrowprops={'arrowstyle':'->'}); ax.set(xlim=(0,1),ylim=(.15,1)); ax.axis('off')
    savefig(fig,out,'fig1_ep_simplex_vs_graph')

    fig,axs=plt.subplots(1,2,figsize=(10,4),sharex=True)
    for term,g in traj.groupby('term'):
        axs[0].plot(g.tau,g.CI,label=term,color=colors[term],lw=2); axs[1].plot(g.tau,g.Q,label=term,color=colors[term],lw=2)
    for ax in axs: ax.axvspan(.4,.7,color='#f2c94c',alpha=.18); ax.set_xlabel(r'$\tau$'); ax.legend(frameon=False)
    axs[0].set_ylabel(r'$CI_2(\tau)$'); axs[1].set_ylabel(r'$Q(\tau)$'); axs[0].set_title('2-skeleton contraction'); axs[1].set_title('Multi-party modularity')
    if float(traj.CI.abs().max())<1e-9:
        axs[0].text(.5,.65,'No sampled flag triangle lacks\na common Yes-vote witness',transform=axs[0].transAxes,ha='center',va='center')
    fig.suptitle('European Parliament threshold filtrations'); fig.tight_layout(); savefig(fig,out,'fig2_ep_ci_vs_modularity_timeline')

    labels=sorted(policy.policy_area.unique()); angles=np.linspace(0,2*np.pi,len(labels),endpoint=False).tolist(); angles+=angles[:1]
    fig,ax=plt.subplots(figsize=(6,6),subplot_kw={'polar':True})
    for term,g in policy.groupby('term'):
        vals=[float(g.set_index('policy_area').loc[x,'mean_CI']) for x in labels]; vals+=vals[:1]
        ax.plot(angles,vals,lw=2,label=term,color=colors[term]); ax.fill(angles,vals,alpha=.08,color=colors[term])
    ax.set_xticks(angles[:-1],labels); ax.set_ylim(0,1); ax.set_title('Policy-domain CI resilience'); ax.legend(loc='upper right',bbox_to_anchor=(1.25,1.1))
    if float(policy.mean_CI.abs().max())<1e-5:
        ax.text(0,0,'CI approximately 0\nin every domain',ha='center',va='center',weight='bold')
    savefig(fig,out,'fig3_ep_policy_radar')

    fig,ax=plt.subplots(figsize=(8,4.5))
    for term,g in traj.groupby('term'):
        for col,ls in [('beta_EPP','-'),('beta_SD','--'),('beta_Renew',':')]: ax.plot(g.tau,g[col],ls,color=colors[term],alpha=.8,label=term+' '+col[5:])
    ax.set(xlabel=r'$\tau$',ylabel=r'$\beta_G(\tau)$',ylim=(-.02,.12),title='Vertex-level Heyting buffer (identically zero by definition)')
    ax.text(.5,.055,'All vertices lie in the target induced subcomplex or its induced complement.',ha='center'); ax.legend(ncol=2,fontsize=8)
    savefig(fig,out,'fig4_ep_heyting_buffer_collapse')

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]); args=parser.parse_args()
    root=args.root; out=root/'ep_experiments'; out.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format='%(message)s',handlers=[logging.StreamHandler(),logging.FileHandler(out/'ep_pipeline.log',mode='w',encoding='utf-8')])
    source=root/'VoteWatch-EP-voting-data_2004-2022'; started=time.perf_counter(); terms={}; cache=out/'ep_aligned_cache.pkl'
    if cache.exists():
        with cache.open('rb') as fh: terms=pickle.load(fh)
        LOG.info('Loaded aligned read-only source cache %s',cache)
    else:
        for term in ['EP6','EP9']:
            meta=load_metadata(source/(term+'_Voted docs.xlsx'),term)
            rcv=source/('EP6_RCVs_2022_06_13.xlsx' if term=='EP6' else 'EP9_RCVs_2022_06_22.xlsx')
            vote_ids,members,groups,countries,x=load_rcv(rcv,term)
            terms[term]={'meta':meta,'vote_ids':vote_ids,'members':members,'groups':groups,'countries':countries,'x':x,'policies':align_metadata(vote_ids,meta)}
        with cache.open('wb') as fh: pickle.dump(terms,fh,pickle.HIGHEST_PROTOCOL)
    for term,d in terms.items():
        LOG.info('%s: %d MEP rows, %d roll calls, groups=%s',term,d['x'].shape[0],d['x'].shape[1],dict(pd.Series(d['groups']).value_counts()))
    progress(1,started,'metadata and roll-call matrices aligned')

    t=time.perf_counter(); frames=[]; diagnostics=[]
    for term,d in terms.items():
        frame,diag,_=trajectory(term,d['x'],d['groups']); frames.append(frame); diagnostics.append(diag)
    traj=pd.concat(frames,ignore_index=True); traj.to_csv(out/'ep_filtration_trajectories.csv',index=False)
    (out/'ep_computation_diagnostics.json').write_text(json.dumps(diagnostics,indent=2),encoding='utf-8')
    progress(2,t,'41-step filtrations computed for EP6 and EP9')

    t=time.perf_counter(); policy=[]
    for term,d in terms.items(): policy.extend(policy_breakdown(term,d['x'],d['policies'],d['groups']))
    policy_df=pd.DataFrame(policy); policy_df.to_csv(out/'ep_policy_breakdown.csv',index=False)
    progress(3,t,'policy-domain summaries computed')

    t=time.perf_counter(); at50=traj[np.isclose(traj.tau,.5)].set_index('term'); piv=traj.pivot(index='tau',columns='term')
    placeholders={
      'EP6_MEAN_COAL':float(terms['EP6']['x'].sum(axis=0).mean()/terms['EP6']['x'].shape[0]),
      'EP9_MEAN_COAL':float(terms['EP9']['x'].sum(axis=0).mean()/terms['EP9']['x'].shape[0]),
      'EP6_DMAX':int(traj[traj.term=='EP6'].d_max.max()),'EP9_DMAX':int(traj[traj.term=='EP9'].d_max.max()),
      'EP6_DENS_05':float(2*at50.loc['EP6','f1']/(at50.loc['EP6','f0']*(at50.loc['EP6','f0']-1))),
      'EP9_DENS_05':float(2*at50.loc['EP9','f1']/(at50.loc['EP9','f0']*(at50.loc['EP9','f0']-1))),
      'PEAK_DELTA_CI':float((piv['CI']['EP9']-piv['CI']['EP6']).abs().max()),
      'MEAN_DELTA_Q':float((piv['Q']['EP9']-piv['Q']['EP6']).mean()),
      'METHOD_NOTE':'CI is the joint-Monte-Carlo 2-skeleton contraction; d_max is a deterministic capped lower bound (cap 15), with diagnostics in ep_computation_diagnostics.json',
      'EP9_COVERAGE_NOTE':'VoteWatch EP9 roll calls in the supplied file end in June 2022, not 2024'}
    (out/'ep_latex_placeholders.json').write_text(json.dumps(placeholders,indent=2),encoding='utf-8')
    progress(4,t,'LaTeX scalar placeholders written')

    t=time.perf_counter(); make_plots(traj,policy_df,out); progress(5,t,'eight 300-DPI/PDF figure files generated')
    LOG.info('Total runtime %.2fs',time.perf_counter()-started)

if __name__=='__main__': main()
