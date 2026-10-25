"""Read-only reuse of prior numerical routines; all new writes isolated here."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,time,csv,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;OLD=ROOT/'results'
sys.dont_write_bytecode=True
sys.path.insert(0,str(OLD))
import analysis_core as core
core.OUT=OUT
from analysis_core import Metrics,Depth,TAUS
SEED=42
def log(stage,message):
    print(f'[INFO] {stage}: {message}',flush=True)
    path=OUT/'run_log.csv';exists=path.exists()
    with path.open('a',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        if not exists:writer.writerow(['timestamp','stage','seed','message'])
        writer.writerow([time.strftime('%Y-%m-%dT%H:%M:%S'),stage,SEED,message])
def save(name,frame,replace=False):
    path=OUT/name
    if path.exists() and not replace:raise FileExistsError(path)
    temp=path.with_suffix('.pending.csv');frame.to_csv(temp,index=False);temp.replace(path)
def inputs():
    z=np.load(ROOT/'results_v2/input_EP9_fixed.npz')
    assert np.array_equal(z['present'],np.isin(z['states'],[1,2,3]))
    assert not np.any(z['yes']&z['no'])
    return z,z['states'],Metrics(z['groups'])
def delta(frame):
    f=frame.sort_values('tau').copy()
    f['delta_Q']=f.Q.diff();f['delta_d']=f.d_max.diff();f['delta_d_greedy']=f.d_greedy.diff()
    f['delta_d_lower']=f.d_lower-f.d_upper.shift();f['delta_d_upper']=f.d_upper-f.d_lower.shift()
    return f
def selected(frame):
    f=delta(frame)
    metric='delta_d' if f.delta_d.notna().any() else 'delta_d_greedy'
    i=f[metric].idxmin();r=f.loc[i];prior=f[f.tau<r.tau].iloc[-1]
    # Bounds can sometimes certify a global minimum even with some timeouts.
    certified=bool(pd.notna(r.delta_d) and (f.loc[f.index!=i,'delta_d_lower'].dropna()>=r.delta_d).all())
    selection='certified_global_drop' if certified else 'largest_certified_step_global_drop_unverified' if metric=='delta_d' else 'greedy_candidate_no_certified_steps'
    return dict(tau_before=float(prior.tau),tau_after=float(r.tau),delta_Q=float(r.delta_Q),delta_d=r.delta_d,delta_d_lower=r.delta_d_lower,delta_d_upper=r.delta_d_upper,delta_d_greedy=r.delta_d_greedy,selection=selection,before_status=prior.status,after_status=r.status)
def fingerprints(verify=False):
    paths=[]
    for folder in [OLD,ROOT/'results_v2']:
        paths+=sorted(p for p in folder.rglob('*') if p.is_file())
    paths+=[ROOT/'ep_aligned_cache.pkl',ROOT/'ep_topological_pipeline.py']
    paths+=sorted((ROOT.parents[1]/'elsevier_submission').glob('*.tex'))
    records=[]
    for p in paths:
        with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        records.append(dict(path=str(p),sha256=digest))
    frame=pd.DataFrame(records)
    if verify:
        assert frame.to_dict('records')==pd.read_csv(OUT/'source_fingerprints.csv').to_dict('records'),'An old input/output changed'
        log('preservation','all original input/output/manuscript hashes unchanged')
    else:save('source_fingerprints.csv',frame)
