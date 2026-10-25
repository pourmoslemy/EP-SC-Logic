import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='2'
import sys,time,pickle,csv
from pathlib import Path
import numpy as np
import pandas as pd
from analysis_core import OUT,TAUS,Metrics,Depth
sys.dont_write_bytecode=True
def log(phase,term,mode,detail):
    message=f'{phase} {term} {mode}: {detail}';print('[INFO]',message,flush=True)
    path=OUT/'run_log.csv';exists=path.exists()
    with path.open('a',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        if not exists:w.writerow(['timestamp','phase','term','mode','detail'])
        w.writerow([time.strftime('%Y-%m-%dT%H:%M:%S'),phase,term,mode,detail])
def main():
    if (OUT/'observed.csv').exists():return
    depth=Depth();frames=[];meshes=[];density=[];mains=[];oldcomparison=[]
    with (OUT.parent/'ep_aligned_cache.pkl').open('rb') as f:original=pickle.load(f)
    old=pd.read_csv(OUT.parent/'ep_filtration_trajectories.csv')
    for term in ['EP6','EP9']:
        z=np.load(OUT/f'input_{term}.npz');m=Metrics(z['groups'])
        for mode in ['missing','current']:
            a=z[mode];num,den=m.counts(a);assert np.all(num<=den)
            for minimum in [20,10,50]:
                ts=np.round(np.arange(.1,.90001,.01),2) if minimum==20 else TAUS
                frame=m.trajectory(m.omega(num,den,minimum),ts,depth=depth,budget=.5)
                frame['term']=term;frame['mode']=mode;frame['minimum_copresent']=minimum
                common=frame[frame.tau.isin(TAUS)].copy();common['delta_Q']=common.Q.diff();common['delta_d']=common.d_max.diff();common['delta_d_greedy']=common.d_greedy.diff()
                frames.append(common)
                row=common[np.isclose(common.tau,.5)].iloc[0];oldrow=old[(old.term==term)&np.isclose(old.tau,.5)].iloc[0]
                density.append(dict(term=term,mode=mode,minimum_copresent=minimum,old_density=oldrow.f1/(len(m.u)),new_density=row.density,old_edges=oldrow.f1,new_edges=row.edges))
                for r in common.itertuples():
                    oldr=old[(old.term==term)&np.isclose(old.tau,r.tau)].iloc[0]
                    oldcomparison.append(dict(term=term,mode=mode,minimum_copresent=minimum,tau=r.tau,old_capped_d=oldr.d_max,new_exact_d=r.d_max,new_lower=r.d_lower,new_upper=r.d_upper,status=r.status))
                if minimum==20:
                    for step in [.01,.02,.05]:
                        target=np.round(np.arange(.1,.90001,step),2);mesh=frame[frame.tau.isin(target)].copy()
                        mesh['delta_Q']=mesh.Q.diff();mesh['delta_d']=mesh.d_max.diff();mesh['delta_d_greedy']=mesh.d_greedy.diff();mesh['mesh_step']=step;meshes.append(mesh)
                        certified=mesh.status.eq('exact').all()
                        metric='delta_d' if certified else 'delta_d_greedy';loc=mesh[metric].idxmin();r=mesh.loc[loc]
                        prev=mesh[mesh.tau<r.tau].iloc[-1]
                        mains.append(dict(term=term,mode=mode,mesh_step=step,tau_before=prev.tau,tau_after=r.tau,delta_Q=r.delta_Q,delta_d=r.delta_d,delta_d_greedy=r.delta_d_greedy,selection='exact_global_drop' if certified else 'greedy_selected_global_drop_unverified',before_status=prev.status,after_status=r.status))
                depth.save();log('observed',term,mode,f'minimum={minimum}; exact={(frame.status=="exact").sum()}/{len(frame)}')
    pd.concat(frames,ignore_index=True).to_csv(OUT/'observed.csv',index=False)
    pd.concat(meshes,ignore_index=True).to_csv(OUT/'mesh.csv',index=False)
    pd.DataFrame(density).to_csv(OUT/'density_comparison.csv',index=False)
    pd.DataFrame(oldcomparison).to_csv(OUT/'depth_vs_old.csv',index=False)
    pd.DataFrame(mains).to_csv(OUT/'main_steps.csv',index=False)
    pd.DataFrame([dict(phase='party_null',seed=42,draws=200,minimum_copresent=20,note='Runtime fallback to 200 draws at four thresholds around each term/mode main step; fixed presence, within-party roll-call vote permutations'),dict(phase='block_bootstrap',seed=42,draws=2000,minimum_copresent=20,note='Lengths 20,50,100; resample nonoverlapping contiguous blocks, retaining short terminal block; uncapped fixed-order greedy depth approximation; Q/omega recomputed')]).to_csv(OUT/'seeds_and_methods.csv',index=False)
if __name__=='__main__':main()
