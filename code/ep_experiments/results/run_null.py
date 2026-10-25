import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import time
import numpy as np
import pandas as pd
from analysis_core import OUT,Metrics,Depth
from run_observed import log
N=200;SEED=42
def event_ge(lo,hi,obslo,obshi):
    yes=lo>=obshi;no=hi<obslo
    return yes,~(yes|no)
def probability(yes,unknown):
    low=(1+np.count_nonzero(yes))/(N+1);high=(1+np.count_nonzero(yes|unknown))/(N+1)
    return (low if low==high else np.nan),low,high
def main():
    observed=pd.read_csv(OUT/'observed.csv');steps=pd.read_csv(OUT/'main_steps.csv');depth=Depth();allrows=[];ps=[];joint=[]
    for term in ['EP6','EP9']:
        z=np.load(OUT/f'input_{term}.npz');m=Metrics(z['groups']);prior=None
        for mode in ['missing','current']:
            path=OUT/f'null_draws_{term}_{mode}.csv';a=z[mode];rng=np.random.default_rng(SEED)
            step=steps[(steps.term==term)&(steps['mode']==mode)&np.isclose(steps.mesh_step,.02)].iloc[0]
            low,high=float(step.tau_before),float(step.tau_after)
            ts=np.unique(np.clip(np.round([low-.02,low,high,high+.02],2),.1,.9))
            obs=observed[(observed.term==term)&(observed['mode']==mode)&(observed.minimum_copresent==20)].set_index('tau')
            if path.exists() and pd.read_csv(path,usecols=['replicate']).replicate.nunique()==N:draws=pd.read_csv(path)
            elif prior is not None and np.array_equal(a,z['missing']):
                draws=prior.copy();draws['mode']=mode;draws.to_csv(path,index=False)
            else:
                _,den=m.counts(a);pieces=[];groups=[]
                for group in np.unique(z['groups']):
                    ix=np.flatnonzero(z['groups']==group);aa=a[ix];eligible=aa>0
                    dest=np.argsort(~eligible,axis=0,kind='stable');groups.append((ix,aa,eligible,dest,np.sum(aa==1,axis=0)))
                started=time.perf_counter()
                for rep in range(N):
                    perm=a.copy()
                    for ix,aa,eligible,dest,yescounts in groups:
                        keys=rng.random(aa.shape);keys[~eligible]=np.inf
                        order=np.argsort(keys,axis=0,kind='stable');shuffled=np.take_along_axis(aa,order,axis=0)
                        target=np.empty_like(aa);np.put_along_axis(target,dest,shuffled,axis=0)
                        assert np.array_equal(target>0,eligible)
                        assert np.array_equal(np.sum(target==1,axis=0),yescounts)
                        perm[ix]=target
                    y=(perm==1).astype(np.float32);num=(y@y.T)[m.u,m.v]
                    frame=m.trajectory(m.omega(num,den),ts,depth=depth,budget=.12)
                    frame['replicate']=rep+1;frame['term']=term;frame['mode']=mode;pieces.append(frame)
                    if (rep+1)%20==0:
                        snapshot=pd.concat(pieces,ignore_index=True);pending=path.with_suffix('.pending.csv');snapshot.to_csv(pending,index=False);pending.replace(path);depth.save()
                        log('party_null',term,mode,f'{rep+1}/{N}; seed={SEED}; {time.perf_counter()-started:.1f}s')
                draws=pd.concat(pieces,ignore_index=True)
            # An interrupted partial table is not a completed null ensemble.
            if draws.replicate.nunique()!=N:
                raise RuntimeError(f'Incomplete {path.name}: rerun after preserving/renaming partial output; seed replay required')
            prior=draws.copy() if mode=='missing' else prior
            allrows.append(draws)
            for tau in ts:
                g=draws[np.isclose(draws.tau,tau)];o=obs.loc[tau]
                yes,unknown=event_ge(g.d_lower.to_numpy(),g.d_upper.to_numpy(),o.d_lower,o.d_upper)
                p,plo,phi=probability(yes,unknown)
                ps.append(dict(term=term,mode=mode,tau=tau,N=N,obs_d=o.d_max,obs_d_lower=o.d_lower,obs_d_upper=o.d_upper,obs_Q=o.Q,null_Q_mean=g.Q.mean(),p_Q=(1+sum(g.Q>=o.Q-1e-12))/(N+1),p_d=p,p_d_lower=plo,p_d_upper=phi,ambiguous_d_comparisons=int(unknown.sum()),exact_null_depths=sum(g.status=='exact'),null_d_mean=g.d_max.mean() if g.status.eq('exact').all() else np.nan,null_d_lower_mean=g.d_lower.mean(),null_d_upper_mean=g.d_upper.mean()))
            before=draws[np.isclose(draws.tau,low)].set_index('replicate');after=draws[np.isclose(draws.tau,high)].set_index('replicate')
            dq=after.Q-before.Q;dl=after.d_lower-before.d_upper;du=after.d_upper-before.d_lower
            odq=float(obs.loc[high].Q-obs.loc[low].Q);odl=float(obs.loc[high].d_lower-obs.loc[low].d_upper);odu=float(obs.loc[high].d_upper-obs.loc[low].d_lower)
            qevent=(dq.to_numpy()>=odq-1e-12)
            yes=qevent&(du.to_numpy()<=odl);unknown=qevent&~((du.to_numpy()<=odl)|(dl.to_numpy()>odu))
            pj,pl,pu=probability(yes,unknown)
            joint.append(dict(term=term,mode=mode,tau_before=low,tau_after=high,N=N,obs_delta_Q=odq,obs_delta_d=step.delta_d,obs_delta_d_lower=odl,obs_delta_d_upper=odu,obs_delta_d_greedy=step.delta_d_greedy,null_delta_Q_mean=dq.mean(),null_delta_d_lower_mean=dl.mean(),null_delta_d_upper_mean=du.mean(),p_joint=pj,p_joint_lower=pl,p_joint_upper=pu,ambiguous_joint_comparisons=int(unknown.sum()),main_step_selection=step.selection))
            details=pd.DataFrame(dict(replicate=before.index,delta_Q=dq,delta_d_lower=dl,delta_d_upper=du));details.to_csv(OUT/f'null_step_draws_{term}_{mode}.csv',index=False)
    pd.concat(allrows,ignore_index=True).to_csv(OUT/'null_draws.csv',index=False)
    pd.DataFrame(ps).to_csv(OUT/'null_pvalues.csv',index=False);pd.DataFrame(joint).to_csv(OUT/'null_joint.csv',index=False);depth.save()
if __name__=='__main__':main()
