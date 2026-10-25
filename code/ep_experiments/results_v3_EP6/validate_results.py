"""Independent saved-output and direct-definition checks; no old writes."""
from common import *
import networkx as nx
def main():
    records=[]
    def check(label,condition):
        records.append(dict(check=label,passed=bool(condition)));assert condition,label
    obs=pd.read_csv(OUT/'observed.csv');mesh=pd.read_csv(OUT/'mesh.csv');steps=pd.read_csv(OUT/'main_steps.csv')
    check('observed 3 minima by 41 rows',len(obs)==123)
    check('mesh 81+41+17 rows',len(mesh)==139)
    for filename,frame,keys in [('observed',obs,['minimum_copresent']),('mesh',mesh,['mesh_step'])]:
        check(filename+' all timeouts have both bounds',bool(frame[['d_lower','d_upper']].notna().all().all()))
        check(filename+' bound ordering',bool((frame.d_lower<=frame.d_upper).all()))
        exact=frame.status.eq('exact')
        check(filename+' timeout is not exact',bool(frame.loc[~exact,'d_max'].isna().all()))
        check(filename+' exact bounds agree',bool((frame.loc[exact,'d_lower']==frame.loc[exact,'d_max']).all() and (frame.loc[exact,'d_upper']==frame.loc[exact,'d_max']).all()))
        for key,g in frame.groupby(keys):
            g=g.sort_values('tau');e=g[g.status=='exact']
            check(f'{filename} {key} monotone certified depth',bool((e.d_max.diff().dropna()<=0).all()))
            check(f'{filename} {key} all deltaQ',np.allclose(g.delta_Q.iloc[1:],g.Q.diff().iloc[1:],atol=1e-12))
            mask=(g.tau>=.38-1e-8)&(g.tau<=.52+1e-8)
            check(f'{filename} {key} requested regional budget',bool((g.loc[mask,'requested_budget_seconds']==10).all()))
    for r in steps.itertuples():
        g=mesh[np.isclose(mesh.mesh_step,r.mesh_step)];g=g.sort_values('tau');sel=g[np.isclose(g.tau,r.tau_after)].iloc[0]
        check(f'mesh {r.mesh_step} selection largest certified step',np.isclose(r.delta_d,g.delta_d.min()))
        if r.selection=='certified_global_drop':
            competitor=g[~np.isclose(g.tau,r.tau_after)].delta_d_lower.dropna()
            check(f'mesh {r.mesh_step} GLOBAL certificate',bool((competitor>=r.delta_d).all()))
    z,a,m=inputs();num,den=m.counts(a)
    rng=np.random.default_rng(SEED)
    for ix in rng.choice(len(m.u),12,replace=False):
        u,v=m.u[ix],m.v[ix]
        check(f'pair {u},{v} co-presence',den[ix]==np.sum(np.isin(a[u],[1,2,3])&np.isin(a[v],[1,2,3])))
        check(f'pair {u},{v} co-Yes',num[ix]==np.sum((a[u]==1)&(a[v]==1)))
    w=m.omega(num,den,20)
    for tau in [.42,.44,.50]:
        keep=w>=tau;g=nx.Graph();g.add_nodes_from(range(m.n));g.add_weighted_edges_from(zip(m.u[keep],m.v[keep],w[keep]))
        partition=[set(np.flatnonzero(z['groups']==label)) for label in np.unique(z['groups'])]
        q=nx.algorithms.community.modularity(g,partition,weight='weight')
        row=obs[(obs.minimum_copresent==20)&np.isclose(obs.tau,tau)].iloc[0]
        check(f'Q {tau} vs independent networkx',np.isclose(q,row.Q,atol=1e-12))
        check(f'edges {tau} vs direct threshold',g.number_of_edges()==row.edges)
    ps=pd.read_csv(OUT/'null_pvalues.csv') if (OUT/'null_pvalues.csv').exists() else pd.DataFrame()
    draws=pd.read_csv(OUT/'null_draws.csv') if not ps.empty else pd.DataFrame()
    for r in ps.itertuples():
        g=draws[np.isclose(draws.tau,r.tau)]
        check(f'null {r.tau} 1000 unique draws',len(g)==1000 and g.replicate.nunique()==1000)
        check(f'null {r.tau} pQ recomputed',np.isclose(r.p_Q,(1+np.sum(g.Q>=r.obs_Q-1e-12))/1001))
        yes=g.d_lower>=r.obs_d_upper;no=g.d_upper<r.obs_d_lower;amb=~(yes|no)
        check(f'null {r.tau} pd bounds recomputed',np.allclose([r.p_d_lower,r.p_d_upper],[(1+yes.sum())/1001,(1+yes.sum()+amb.sum())/1001]))
        check(f'null {r.tau} ambiguous recomputed',int(amb.sum())==r.ambiguous_count)
    check('UNKNOWN group distinct from NI',int(np.sum(z['groups']=='UNKNOWN'))==0 and int(np.sum(z['present'][z['groups']=='UNKNOWN']))==0)
    six=np.load(OLD/'input_EP6.npz')
    check('EP6 Yes disagreement independently rechecked',np.count_nonzero((six['current']==1)!=(six['coordinate']==1))==0)
    check('stored cells NOT presence',np.sum(z['stored_vote_cells']&~z['present'])==2066039)
    if ps.empty:log('verification','null p-value tests pending; current run validates observed outputs only')
    ci_path=OUT/'bootstrap_CI.csv'
    if ci_path.exists():
        ci=pd.read_csv(ci_path)
        check('complete bootstrap CI row count',len(ci)==18)
        for length,g in ci.groupby('block_length'):
            draws=pd.read_csv(OUT/f'bootstrap_draws_EP6_fixed_{length}.csv')
            check(f'L{length} 2000 unique bootstrap draws',len(draws)==2000 and draws.replicate.nunique()==2000)
            for r in g.itertuples():
                vals=draws[r.metric].dropna()
                check(f'L{length} {r.metric} percentile CI',np.allclose([r.CI_025,r.CI_975],vals.quantile([.025,.975]),atol=1e-10) and r.valid_replicates==len(vals))
    fingerprints(verify=True)
    save('output_validation.csv',pd.DataFrame(records),replace=(OUT/'output_validation.csv').exists())
    log('verification',f'PASS {len(records)} independent result checks')
if __name__=='__main__':main()
