"""Run last: fixed EP9, exact full-pair block sums, declared greedy depth."""
from common import *
import pickle
import run_bootstrap as reused
from bootstrap_core import curve,statistics
N=2000
# Keep the v3 log schema; the reused helper otherwise writes legacy 5-column logs.
reused.log=lambda phase,term,mode,detail:log(phase,f'{term} {mode}: {detail}')
def main():
    assert (OUT/'null_pvalues.csv').exists(),'Required order is observed, null, bootstrap'
    if (OUT/'bootstrap_CI.csv').exists() and len(pd.read_csv(OUT/'bootstrap_CI.csv'))==18:
        raise FileExistsError('Completed bootstrap outputs already exist')
    z,a,m=inputs()
    with (ROOT/'ep_aligned_cache.pkl').open('rb') as f:d=pickle.load(f)['EP9']
    d['vote_ids']=z['vote_ids']
    order_time=reused.chronological('EP9',d);a=a[:,order_time]
    num,den=m.counts(a);w=m.omega(num,den,20)
    matrix=np.zeros((m.n,m.n));matrix[m.u,m.v]=w;matrix[m.v,m.u]=w
    order=np.argsort(-matrix.sum(1),kind='stable');del matrix
    step=pd.read_csv(OUT/'main_steps.csv');step=step[np.isclose(step.mesh_step,.02)].iloc[0]
    low=int(np.flatnonzero(np.isclose(TAUS,step.tau_before))[0]);high=int(np.flatnonzero(np.isclose(TAUS,step.tau_after))[0])
    legacy_low=int(np.flatnonzero(np.isclose(TAUS,.44))[0]);legacy_high=int(np.flatnonzero(np.isclose(TAUS,.46))[0])
    def stats(q,d):
        r=statistics(q,d,low,high)
        r['delta_Q_at_legacy_step']=float(q[legacy_high]-q[legacy_low]);r['delta_d_approx_at_legacy_step']=float(d[legacy_high]-d[legacy_low])
        return r
    q,ds=curve(m,w,order);point=stats(q,ds)
    if not (OUT/'bootstrap_point.csv').exists():save('bootstrap_point.csv',pd.DataFrame([dict(seed=SEED,tau_before=step.tau_before,tau_after=step.tau_after,**point)]))
    ci=[]
    for length in [20,50,100]:
        path=OUT/f'bootstrap_draws_EP9_fixed_{length}.csv';rows=pd.read_csv(path).to_dict('records') if path.exists() else []
        started=time.perf_counter()
        if len(rows)<N:
            npth,dpth=reused.counts_for_replicates('EP9','fixed',length,a,m)
            ns=np.load(npth,mmap_mode='r');dens=np.load(dpth,mmap_mode='r')
            for rep in range(len(rows),N):
                ww=m.omega(ns[rep],dens[rep],20);qq,dd=curve(m,ww,order)
                rows.append(dict(replicate=rep+1,**stats(qq,dd)))
                if (rep+1)%50==0:
                    save(path.name,pd.DataFrame(rows),replace=path.exists());log('bootstrap',f'EP9 fixed L={length}, {rep+1}/{N}; seed={SEED}; {time.perf_counter()-started:.1f}s')
            del ns,dens
            # Only task-created, regenerable work arrays inside v3/_cache.
            for pathwork in [npth,dpth,reused.CACHE/f'ready_EP9_fixed_{length}.csv']:
                assert pathwork.resolve().parent==reused.CACHE.resolve();pathwork.unlink()
        frame=pd.DataFrame(rows);assert len(frame)==N
        for metric in point:
            values=frame[metric].dropna()
            ci.append(dict(term='EP9',mode='fixed',block_length=length,seed=SEED,replicates=N,valid_replicates=len(values),metric=metric,point_estimate=point[metric],CI_025=values.quantile(.025),CI_975=values.quantile(.975),tau_before=.44 if metric.endswith('legacy_step') else step.tau_before,tau_after=.46 if metric.endswith('legacy_step') else step.tau_after,depth_approximation='Uncapped fixed observed-strength order, two-start greedy, backward monotone lower-bound envelope; NOT exact d'))
        save('bootstrap_CI.csv',pd.DataFrame(ci),replace=(OUT/'bootstrap_CI.csv').exists())
        log('bootstrap-case',f'L={length} complete with {N} resamples')
    fingerprints(verify=True);log('stage-2','complete')
if __name__=='__main__':main()
