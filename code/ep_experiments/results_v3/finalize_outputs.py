"""Metadata alignment and final manifest; writes only inside results_v3."""
from common import *
def main():
    methods=pd.read_csv(OUT/'methods.csv')
    methods.loc[methods.key=='selection','value']='Preliminary greedy budget placement; final largest certified exact step, globally certified only if competing delta lower bounds exclude a larger drop; see selection_methods.csv'
    additions={'density_denominator':'811 choose 2, all cohort members including isolates; active nodes reported separately','Q_when_edgeless':'Q=0 by explicit computational convention','p_value_formula':'(1+count(null >= observed))/(1000+1); Q tolerance 1e-12, d bounds and ambiguous comparisons explicit','stored_vote_cells_definition':'True iff physical raw vote value is nonempty numeric; not a presence or missingness mask; codes 0,4,5,6 are NOT present','EP6_Yes_alignment_disagreements':'0 on all shared cached member/vote IDs','bootstrap_CI_method':'Percentile 0.025 and 0.975 over 2000 full-pair block resamples; approximate-depth intervals are NOT exact-depth intervals','timeout_reproducibility':'Exact results deterministic; timeout bounds may depend on wall-clock scheduling/hardware. Per-graph budgets, elapsed seconds and certified bounds retained.'}
    for key,value in additions.items():
        if key not in set(methods.key):methods.loc[len(methods)]=dict(key=key,value=value)
    save('methods.csv',methods,replace=True)
    fingerprints(verify=True)
    rec=[]
    for path in sorted(OUT.iterdir()):
        if not path.is_file() or path.name=='output_manifest.csv':continue
        with path.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        rec.append(dict(file=path.name,bytes=path.stat().st_size,sha256=digest))
    save('output_manifest.csv',pd.DataFrame(rec),replace=(OUT/'output_manifest.csv').exists())
    print('Final metadata and manifest saved',flush=True)
if __name__=='__main__':main()
