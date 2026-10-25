"""1000 seeded party-aware null draws; exact Q, certified depth p-bounds."""
from common import *
import json
N=1000
def probability(yes,unknown):
    lo=(1+int(np.sum(yes)))/(N+1);hi=(1+int(np.sum(yes|unknown)))/(N+1)
    return (lo if lo==hi else np.nan),lo,hi
def main():
    assert (OUT/'observed.csv').exists(),'Execute stage 1 first'
    if (OUT/'null_pvalues.csv').exists():raise FileExistsError('Completed null outputs already exist')
    z,a,m=inputs();observed=pd.read_csv(OUT/'observed.csv');observed=observed[observed.minimum_copresent==20].set_index('tau')
    steps=pd.read_csv(OUT/'main_steps.csv');step=steps[np.isclose(steps.mesh_step,.02)].iloc[0]
    taus=np.unique(np.r_[[.40,.44,.46,.50],step.tau_before,step.tau_after]).round(2)
    rng=np.random.default_rng(SEED);depth=Depth();pieces=[];completed=0
    draws_path=OUT/'null_draws.csv';state_path=OUT/'null_rng_state.csv'
    if draws_path.exists():
        frame=pd.read_csv(draws_path);state=pd.read_csv(state_path).iloc[0]
        completed=int(state.completed_draws)
        assert frame.replicate.nunique()==completed and len(frame)==completed*len(taus),'Inconsistent null checkpoint; no silent replay'
        rng.bit_generator.state=json.loads(state.rng_state_json);pieces=[frame]
    _,den=m.counts(a);groups=[];audits=[]
    for group in np.unique(z['groups']):
        ix=np.flatnonzero(z['groups']==group);aa=a[ix];eligible=aa>0
        dest=np.argsort(~eligible,axis=0,kind='stable')
        totals=np.stack([np.sum(aa==code,axis=0) for code in [1,2,3]])
        groups.append((group,ix,aa,eligible,dest,totals))
        audits.append(dict(group=group,members=len(ix),present_cells=int(eligible.sum()),yes_cells=int(totals[0].sum()),absent_missing_fixed=True))
    save('null_party_preservation.csv',pd.DataFrame(audits),replace=(OUT/'null_party_preservation.csv').exists())
    log('party-null',f'N={N}, seed={SEED}, taus={list(taus)}, start={completed}; UNKNOWN is separate group; exact search budget=.12s/null graph, bounds retained')
    started=time.perf_counter()
    for rep in range(completed,N):
        perm=a.copy()
        for group,ix,aa,eligible,dest,totals in groups:
            if not eligible.any():continue
            keys=rng.random(aa.shape);keys[~eligible]=np.inf
            order=np.argsort(keys,axis=0,kind='stable');shuffled=np.take_along_axis(aa,order,axis=0)
            target=np.empty_like(aa);np.put_along_axis(target,dest,shuffled,axis=0)
            assert np.array_equal(target>0,eligible)
            for j,code in enumerate([1,2,3]):assert np.array_equal(np.sum(target==code,axis=0),totals[j])
            perm[ix]=target
        y=(perm==1).astype(np.float32);num=(y@y.T)[m.u,m.v]
        frame=m.trajectory(m.omega(num,den,20),taus,depth=depth,budget=.12)
        frame['replicate']=rep+1;frame['term']='EP9';frame['mode']='fixed';frame['seed']=SEED
        pieces.append(frame)
        if (rep+1)%10==0:
            allrows=pd.concat(pieces,ignore_index=True)
            save('null_draws.csv',allrows,replace=draws_path.exists())
            save('null_rng_state.csv',pd.DataFrame([dict(completed_draws=rep+1,seed=SEED,rng_state_json=json.dumps(rng.bit_generator.state))]),replace=state_path.exists())
            depth.save();log('party-null',f'{rep+1}/{N}; {time.perf_counter()-started:.1f}s since resume')
    draws=pd.concat(pieces,ignore_index=True);ps=[]
    for tau in taus:
        g=draws[np.isclose(draws.tau,tau)];o=observed.loc[tau]
        yes=g.d_lower.to_numpy()>=o.d_upper;no=g.d_upper.to_numpy()<o.d_lower
        unknown=~(yes|no);p,lo,hi=probability(yes,unknown)
        ps.append(dict(term='EP9',tau=tau,N=N,seed=SEED,obs_Q=o.Q,obs_d=o.d_max,obs_d_lower=o.d_lower,obs_d_upper=o.d_upper,p_Q=(1+int(np.sum(g.Q>=o.Q-1e-12)))/(N+1),p_d=p,p_d_lower=lo,p_d_upper=hi,ambiguous_count=int(unknown.sum()),null_exact_depths=int(np.sum(g.status=='exact')),null_depth_timeouts=int(np.sum(g.status=='timeout')),null_Q_mean=g.Q.mean(),null_d_lower_mean=g.d_lower.mean(),null_d_upper_mean=g.d_upper.mean(),Q_comparison_tolerance=1e-12))
    save('null_pvalues.csv',pd.DataFrame(ps))
    # Fixed observed-selected step, not a new step selected in each null draw.
    low=draws[np.isclose(draws.tau,step.tau_before)].set_index('replicate');high=draws[np.isclose(draws.tau,step.tau_after)].set_index('replicate')
    dq=high.Q-low.Q;dl=high.d_lower-low.d_upper;du=high.d_upper-low.d_lower
    save('null_step_draws.csv',pd.DataFrame(dict(replicate=low.index,delta_Q=dq,delta_d=high.d_max-low.d_max,delta_d_lower=dl,delta_d_upper=du)))
    qe=dq.to_numpy()>=step.delta_Q-1e-12
    yes=qe&(du.to_numpy()<=step.delta_d_lower);unknown=qe&~((du.to_numpy()<=step.delta_d_lower)|(dl.to_numpy()>step.delta_d_upper))
    p,lo,hi=probability(yes,unknown)
    save('null_joint.csv',pd.DataFrame([dict(tau_before=step.tau_before,tau_after=step.tau_after,N=N,p_joint=p,p_joint_lower=lo,p_joint_upper=hi,ambiguous_count=int(unknown.sum()),obs_delta_Q=step.delta_Q,obs_delta_d=step.delta_d,obs_delta_d_lower=step.delta_d_lower,obs_delta_d_upper=step.delta_d_upper)]))
    depth.save();fingerprints(verify=True);log('stage-3','complete: 1000 draws at every selected threshold')
if __name__=='__main__':main()
