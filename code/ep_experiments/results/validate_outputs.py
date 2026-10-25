"""Independent checks on saved numerical outputs, without rerunning ensembles."""
from pathlib import Path
import numpy as np
import pandas as pd
O=Path(__file__).resolve().parent
checks=[]
def check(name,condition):
    checks.append(dict(check=name,passed=bool(condition)))
    assert condition,name
obs=pd.read_csv(O/'observed.csv');mesh=pd.read_csv(O/'mesh.csv')
check('observed row count',len(obs)==492)
check('mesh row count',len(mesh)==556)
for name,frame in [('observed',obs),('mesh',mesh)]:
    check(name+' valid depth bounds',bool((frame.d_lower<=frame.d_upper).all()))
    exact=frame.status.eq('exact')
    check(name+' exact values only',bool(frame.loc[~exact,'d_max'].isna().all()))
    check(name+' certified exact bounds',bool((frame.loc[exact,'d_max']==frame.loc[exact,'d_lower']).all() and (frame.loc[exact,'d_max']==frame.loc[exact,'d_upper']).all()))
ci=pd.read_csv(O/'bootstrap_CI.csv')
check('bootstrap 48 CI rows',len(ci)==48)
for (term,mode,length),g in ci.groupby(['term','mode','block_length']):
    draws=pd.read_csv(O/f'bootstrap_draws_{term}_{mode}_{length}.csv')
    check(f'{term} {mode} L{length} 2000 unique draws',len(draws)==2000 and draws.replicate.nunique()==2000)
    for r in g.itertuples():
        vals=draws[r.metric].dropna()
        check(f'{term} {mode} L{length} {r.metric} percentiles',np.allclose([r.CI_025,r.CI_975],vals.quantile([.025,.975]),atol=1e-10) and r.valid_replicates==len(vals))
null=pd.read_csv(O/'null_draws.csv');ps=pd.read_csv(O/'null_pvalues.csv')
for r in ps.itertuples():
    g=null[(null.term==r.term)&(null['mode']==r.mode)&np.isclose(null.tau,r.tau)]
    yes=g.d_lower>=r.obs_d_upper
    unknown=~yes & ~(g.d_upper<r.obs_d_lower)
    check(f'{r.term} {r.mode} {r.tau} null N',len(g)==200 and g.replicate.nunique()==200)
    check(f'{r.term} {r.mode} {r.tau} pQ',np.isclose(r.p_Q,(1+(g.Q>=r.obs_Q-1e-12).sum())/201))
    check(f'{r.term} {r.mode} {r.tau} pd bounds',np.allclose([r.p_d_lower,r.p_d_upper],[(1+yes.sum())/201,(1+yes.sum()+unknown.sum())/201]))
for r in pd.read_csv(O/'null_joint.csv').itertuples():
    g=pd.read_csv(O/f'null_step_draws_{r.term}_{r.mode}.csv')
    qe=g.delta_Q>=r.obs_delta_Q-1e-12
    yes=qe & (g.delta_d_upper<=r.obs_delta_d_lower)
    unknown=qe & ~((g.delta_d_upper<=r.obs_delta_d_lower)|(g.delta_d_lower>r.obs_delta_d_upper))
    check(f'{r.term} {r.mode} joint bounds',np.allclose([r.p_joint_lower,r.p_joint_upper],[(1+yes.sum())/201,(1+yes.sum()+unknown.sum())/201]))
pd.DataFrame(checks).to_csv(O/'output_validation.csv',index=False)
print(f'PASS: {len(checks)} independent output checks',flush=True)
