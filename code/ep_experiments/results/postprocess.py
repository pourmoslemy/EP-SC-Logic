from pathlib import Path
import hashlib,sys,platform
import numpy as np,pandas as pd
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
def main():
    for filename,keys in [('observed.csv',['term','mode','minimum_copresent']),('mesh.csv',['term','mode','mesh_step'])]:
        path=OUT/filename;frame=pd.read_csv(path);pieces=[]
        for _,g in frame.groupby(keys,sort=False):
            g=g.sort_values('tau').copy();g['delta_d_lower']=g.d_lower-g.d_upper.shift();g['delta_d_upper']=g.d_upper-g.d_lower.shift();pieces.append(g)
        pd.concat(pieces,ignore_index=True).to_csv(path,index=False)
    if (OUT/'null_draws.csv').exists():
        draws=pd.read_csv(OUT/'null_draws.csv');joint=pd.read_csv(OUT/'null_joint.csv')
        for row in joint.itertuples():
            g=draws[(draws.term==row.term)&(draws['mode']==row.mode)]
            lo=g[np.isclose(g.tau,row.tau_before)].set_index('replicate');hi=g[np.isclose(g.tau,row.tau_after)].set_index('replicate')
            frame=pd.DataFrame(dict(replicate=lo.index,delta_Q=hi.Q-lo.Q,delta_d=hi.d_max-lo.d_max,delta_d_lower=hi.d_lower-lo.d_upper,delta_d_upper=hi.d_upper-lo.d_lower))
            frame.to_csv(OUT/f'null_step_draws_{row.term}_{row.mode}.csv',index=False)
    originals=[ROOT/'ep_topological_pipeline.py',ROOT/'compute_ep_ci2_robustness.py',ROOT/'ep_aligned_cache.pkl',ROOT/'ep_no_votes_EP6.npz',ROOT/'ep_no_votes_EP9.npz',ROOT/'ep_filtration_trajectories.csv']
    originals+=list((ROOT.parent/'VoteWatch-EP-voting-data_2004-2022').glob('EP[69]_RCVs*.xlsx'))
    records=[]
    for path in originals:
        with path.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        records.append(dict(path=str(path),sha256=digest))
    dest=OUT/'input_fingerprints.csv'
    if dest.exists():assert pd.read_csv(dest).to_dict('records')==records,'Original input fingerprint changed'
    else:pd.DataFrame(records).to_csv(dest,index=False)
    pd.DataFrame([dict(key='seed',value=42),dict(key='python',value=sys.version.split()[0]),dict(key='numpy',value=np.__version__),dict(key='pandas',value=pd.__version__),dict(key='observed_exact_search_seconds',value=.5),dict(key='null_exact_search_seconds',value=.12),dict(key='minimum_copresence',value='20; sensitivity 10,50'),dict(key='null_draws',value=200),dict(key='null_draws_reason',value='Runtime fallback authorized: 1000 full-grid draws require 41 searches per draw versus 4 in selected local test'),dict(key='null_fixed_features',value='Presence/missing masks and per-party per-roll-call Yes,No,Abstain counts'),dict(key='null_Q_comparison_tolerance',value=1e-12),dict(key='bootstrap_reps',value=2000),dict(key='bootstrap_block_lengths',value='50; sensitivity 20,100'),dict(key='bootstrap_depth',value='Approximate: uncapped fixed observed-strength order, two-start greedy, monotone backward lower-bound envelope'),dict(key='bootstrap_counts',value='Full pair counts: exact integer sums of precomputed contiguous-block Gram matrices; equivalent to expanded resampling; no vote subsampling'),dict(key='bootstrap_total_rollcalls',value='Variable: full blocks incl short terminal block resampled intact'),dict(key='missing_dates',value='Known dates sorted stably, unknown appended last; audit in bootstrap_order files'),dict(key='overlap_treatment_basis',value='Original packed Yes/presence basis; primary removes flagged historical cross-cache overlaps, sensitivity leaves current packed states; full coordinate repair NOT performed')]).to_csv(OUT/'runtime_config.csv',index=False)
    print('Postprocessing and input fingerprint check passed',flush=True)
    config=pd.read_csv(OUT/'runtime_config.csv')
    config.loc[len(config)]=dict(key='party_metadata_audit',value='EP9: 20 shifted member rows; 12 cached groups incorrectly NI. Legacy groups retained in Q/null; not coordinate-repaired. See metadata_alignment_audit.csv')
    config.to_csv(OUT/'runtime_config.csv',index=False)
if __name__=='__main__':main()
