import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,time,pickle
from pathlib import Path
import numpy as np
import pandas as pd
from analysis_core import OUT,TAUS,Metrics
from bootstrap_core import curve,statistics
from run_observed import log
N=2000;SEED=42
CACHE=OUT/'_cache';CACHE.mkdir(exist_ok=True)
def chronological(term,d):
    metadata=d['meta'].drop_duplicates('vote_id').set_index('vote_id')
    values=metadata.reindex(d['vote_ids'])['date']
    def parse_date(v):
        if isinstance(v,str):
            if '.' in v:return pd.to_datetime(v,format='%d.%m.%Y',errors='coerce')
            if '/' in v:return pd.to_datetime(v,format='%m/%d/%Y',errors='coerce')
        return pd.to_datetime(v,errors='coerce')
    dates=values.map(parse_date)
    order=pd.DataFrame(dict(source_index=np.arange(len(values)),date=dates.to_numpy())).sort_values(['date','source_index'],na_position='last').source_index.to_numpy()
    pd.DataFrame(dict(term=term,sorted_index=np.arange(len(order)),source_index=order,vote_id=np.asarray(d['vote_ids'])[order],date=dates.iloc[order].to_numpy(),date_unknown=dates.iloc[order].isna().to_numpy())).to_csv(OUT/f'bootstrap_order_{term}.csv',index=False)
    log('bootstrap_order',term,'both',f'unknown dates={dates.isna().sum()}, known-date inversions before sorting={sum(dates.diff().dropna()<pd.Timedelta(0))}; unknown dates last, stable within date')
    return order
def counts_for_replicates(term,mode,length,a,m):
    key=f'{term}_{mode}_{length}';num_path=CACHE/f'num_{key}.npy';den_path=CACHE/f'den_{key}.npy';marker=CACHE/f'ready_{key}.csv'
    if marker.exists() and num_path.exists() and den_path.exists():return num_path,den_path
    b=a.shape[1];starts=list(range(0,b,length));nb=len(starts);pairs=len(m.u)
    bn=CACHE/f'block_num_{key}.npy';bd=CACHE/f'block_den_{key}.npy'
    blocknum=np.lib.format.open_memmap(bn,mode='w+',dtype=np.uint8,shape=(nb,pairs));blockden=np.lib.format.open_memmap(bd,mode='w+',dtype=np.uint8,shape=(nb,pairs))
    for i,start in enumerate(starts):
        part=a[:,start:min(start+length,b)];yes=(part==1).astype(np.float32);present=(part>0).astype(np.float32)
        blocknum[i]=(yes@yes.T)[m.u,m.v].astype(np.uint8);blockden[i]=(present@present.T)[m.u,m.v].astype(np.uint8)
    blocknum.flush();blockden.flush()
    weights=np.random.default_rng(SEED).multinomial(nb,np.full(nb,1/nb),size=N).astype(np.float32)
    lengths=np.array([min(length,b-start) for start in starts]);sample_sizes=weights@lengths
    pd.DataFrame(dict(replicate=np.arange(1,N+1),roll_calls=sample_sizes.astype(int))).to_csv(OUT/f'bootstrap_sizes_{key}.csv',index=False)
    ns=np.lib.format.open_memmap(num_path,mode='w+',dtype=np.uint16,shape=(N,pairs));ds=np.lib.format.open_memmap(den_path,mode='w+',dtype=np.uint16,shape=(N,pairs))
    for start in range(0,pairs,4096):
        end=min(start+4096,pairs)
        nn=weights@np.asarray(blocknum[:,start:end],dtype=np.float32)
        dd=weights@np.asarray(blockden[:,start:end],dtype=np.float32)
        assert np.all(nn<=dd) and np.max(dd)<65536
        ns[:,start:end]=nn.astype(np.uint16);ds[:,start:end]=dd.astype(np.uint16)
    ns.flush();ds.flush();del ns,ds,blocknum,blockden
    pd.DataFrame([dict(seed=SEED,draws=N,block_length=length,blocks=nb,method='Nonoverlapping contiguous chronological blocks; short final block resampled intact; random total roll-call count')]).to_csv(marker,index=False)
    # Only task-created intermediate block Gram arrays are removed.
    for path in [bn,bd]:
        assert path.resolve().parent==CACHE.resolve();path.unlink()
    return num_path,den_path
def main():
    with (OUT.parent/'ep_aligned_cache.pkl').open('rb') as f:original=pickle.load(f)
    main_steps=pd.read_csv(OUT/'main_steps.csv');intervals=[];points=[]
    for term in ['EP6','EP9']:
        z=np.load(OUT/f'input_{term}.npz');m=Metrics(z['groups']);timeorder=chronological(term,original[term])
        for mode in ['missing','current']:
            a=z[mode][:,timeorder];num,den=m.counts(a);w=m.omega(num,den)
            matrix=np.zeros((m.n,m.n));matrix[m.u,m.v]=w;matrix[m.v,m.u]=w;order=np.argsort(-matrix.sum(1),kind='stable');del matrix
            step=main_steps[(main_steps.term==term)&(main_steps['mode']==mode)&np.isclose(main_steps.mesh_step,.02)].iloc[0]
            low=int(np.flatnonzero(np.isclose(TAUS,step.tau_before))[0]);high=int(np.flatnonzero(np.isclose(TAUS,step.tau_after))[0])
            q,ds=curve(m,w,order);point=statistics(q,ds,low,high)
            points.append(dict(term=term,mode=mode,tau_before=step.tau_before,tau_after=step.tau_after,**point))
            for length in [50,20,100]:
                path=OUT/f'bootstrap_draws_{term}_{mode}_{length}.csv';rows=[]
                if path.exists():rows=pd.read_csv(path).to_dict('records')
                if len(rows)<N and mode=='current' and np.array_equal(z['current'],z['missing']):
                    other=OUT/f'bootstrap_draws_{term}_missing_{length}.csv'
                    if other.exists():
                        rows=pd.read_csv(other).to_dict('records');pd.DataFrame(rows).to_csv(path,index=False)
                        source=OUT/f'bootstrap_sizes_{term}_missing_{length}.csv'
                        pd.read_csv(source).to_csv(OUT/f'bootstrap_sizes_{term}_{mode}_{length}.csv',index=False)
                if len(rows)<N:
                    started=time.perf_counter();npth,dpth=counts_for_replicates(term,mode,length,a,m)
                    ns=np.load(npth,mmap_mode='r');dens=np.load(dpth,mmap_mode='r')
                    for rep in range(len(rows),N):
                        ww=m.omega(ns[rep],dens[rep]);qq,dd=curve(m,ww,order)
                        rows.append(dict(replicate=rep+1,**statistics(qq,dd,low,high)))
                        if (rep+1)%50==0:
                            temp=path.with_suffix('.pending.csv');pd.DataFrame(rows).to_csv(temp,index=False);temp.replace(path)
                            log('bootstrap',term,mode,f'L={length}, {rep+1}/{N}, seed={SEED}; {time.perf_counter()-started:.1f}s')
                    del ns,dens
                    # Discard only regenerable bootstrap count work arrays, not data.
                    for work in [npth,dpth,CACHE/f'ready_{term}_{mode}_{length}.csv']:
                        assert work.resolve().parent==CACHE.resolve();work.unlink()
                frame=pd.DataFrame(rows);assert len(frame)==N
                for metric in point:
                    values=frame[metric].dropna()
                    intervals.append(dict(term=term,mode=mode,block_length=length,seed=SEED,replicates=N,valid_replicates=len(values),metric=metric,point_estimate=point[metric],CI_025=values.quantile(.025),CI_975=values.quantile(.975),tau_before=step.tau_before,tau_after=step.tau_after,depth_approximation='Uncapped fixed observed-strength order, two-start greedy, backward monotone lower-bound envelope; not exact depth'))
                pd.DataFrame(intervals).to_csv(OUT/'bootstrap_CI.csv',index=False)
                pd.DataFrame(points).to_csv(OUT/'bootstrap_points.csv',index=False)
                log('bootstrap_case',term,mode,f'L={length} complete; {N} full-pair reweighted pipelines')
if __name__=='__main__':main()
