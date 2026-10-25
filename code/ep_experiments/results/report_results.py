"""Concise descriptive report; no manuscript fitting or silently completed depths."""
from pathlib import Path
import sys,hashlib,platform
import numpy as np,pandas as pd
OUT=Path(__file__).resolve().parent
def read(name):
    try:return pd.read_csv(OUT/name)
    except (FileNotFoundError,pd.errors.EmptyDataError):return pd.DataFrame()
def f(x):return 'NA' if pd.isna(x) else f'{float(x):.4g}'
def interval(row,prefix='p_joint'):
    if pd.notna(row[prefix]):return f(row[prefix])
    return f'[{f(row[prefix+"_lower"])},{f(row[prefix+"_upper"])}]'
def main():
    obs=read('observed.csv');audit=read('overlap_audit.csv');density=read('density_comparison.csv');joint=read('null_joint.csv');pvals=read('null_pvalues.csv');steps=read('main_steps.csv');ci=read('bootstrap_CI.csv')
    lines=['# EP6/EP9 revision - seed 42','P=overlap-missing (primary); C=current packed handling; default co-presence minimum=20.']
    lines+=['','**0. Overlap audit**','','|Term|Overlap|Duplicate member/vote IDs|Yes alignment disagreements|','|---|---:|---|---:|']
    for r in audit.itertuples():lines.append(f'|{r.term}|{r.overlap_cells:,}|{r.duplicate_member_rows}/{r.duplicate_cached_vote_ids}|{r.packed_coordinate_yes_disagreements:,}|')
    lines+=['Cause: coordinate-free Yes packing versus coordinate-indexed No parsing; not duplicates/merge. EP9: 20 shifted rows, including 12 cached group labels incorrectly NI (retained, not repaired).','', '**1. Density at tau=0.50**','','|Term|Old|P/C (min20)|P/C (min10)|P/C (min50)|','|---|---:|---|---|---|']
    for term in ['EP6','EP9']:
        g=density[density.term==term];vals=[]
        for cutoff in [20,10,50]:
            sub=g[g.minimum_copresent==cutoff].set_index('mode');vals.append(f'{f(sub.loc["missing","new_density"])}/{f(sub.loc["current","new_density"])}')
        lines.append(f'|{term}|{f(g.old_density.iloc[0])}|'+'|'.join(vals)+'|')
    lines+=['All pair counts and 41-threshold trajectories are in observed.csv.','','**2. Uncapped exact depth**','','|Term|Exact P/C (of 41)|Timeout P/C|Old maximum|','|---|---|---|---:|']
    for term in ['EP6','EP9']:
        counts=[];timeouts=[]
        for mode in ['missing','current']:
            g=obs[(obs.term==term)&(obs['mode']==mode)&(obs.minimum_copresent==20)]
            counts.append(str(sum(g.status=='exact')));timeouts.append(str(sum(g.status=='timeout')))
        lines.append(f'|{term}|{"/".join(counts)}|{"/".join(timeouts)}|15 (capped)|')
    lines+=['Pivoted Bron-Kerbosch has a 0.5s search budget; timed-out d is NA, with certified bounds retained.','','**3. Party-aware null**','','|Term|Main step*|pQ after P/C|pd after P/C|p joint P/C|','|---|---|---|---|---|']
    for term in ['EP6','EP9']:
        if joint.empty or pvals.empty:lines.append(f'|{term}|pending|pending|pending|pending|');continue
        j=joint[joint.term==term].set_index('mode');q=[];d=[];pj=[]
        for mode in ['missing','current']:
            row=j.loc[mode];p=pvals[(pvals.term==term)&(pvals['mode']==mode)&np.isclose(pvals.tau,row.tau_after)].iloc[0]
            q.append(f(p.p_Q));d.append(interval(p,'p_d'));pj.append(interval(row))
        row=j.loc['missing'];lines.append(f'|{term}|{f(row.tau_before)}->{f(row.tau_after)}|{"/".join(q)}|{"/".join(d)}|{"/".join(pj)}|')
    lines+=['Runtime fallback: 200 draws at four nearby thresholds, presence fixed, within-party/roll-call states shuffled; unresolved exact-depth comparisons are p-intervals.','','**4. Mesh (selected drop endpoints*, P/C)**','','|Term|0.01|0.02|0.05|','|---|---|---|---|']
    for term in ['EP6','EP9']:
        vals=[]
        for step in [.01,.02,.05]:
            g=steps[(steps.term==term)&np.isclose(steps.mesh_step,step)].set_index('mode');vals.append(f'{f(g.loc["missing","tau_after"])}/{f(g.loc["current","tau_after"])}')
        lines.append(f'|{term}|'+'|'.join(vals)+'|')
    lines+=['Every consecutive delta_Q, exact delta_d (or NA), and delta_d bounds are in mesh.csv; *global largest exact drop is unverified where searches time out.','','**5. Block bootstrap (primary, L=50, 95% intervals)**','','|Term|Pearson Q vs d*|Spearman Q vs d*|delta_Q|delta_d*|','|---|---|---|---|---|']
    for term in ['EP6','EP9']:
        g=ci[(ci.term==term)&(ci['mode']=='missing')&(ci.block_length==50)].set_index('metric') if not ci.empty else pd.DataFrame()
        vals=[]
        for key in ['pearson_approx_d','spearman_approx_d','delta_Q','delta_d_approx']:
            vals.append(f'[{f(g.loc[key,"CI_025"])},{f(g.loc[key,"CI_975"])}]' if key in g.index else 'pending')
        lines.append(f'|{term}|'+'|'.join(vals)+'|')
    lines+=['2,000 resamples/case; L=20/50/100 and C intervals are in bootstrap_CI.csv; full-pair omega/Q are recomputed, d* is an uncapped fixed-order greedy monotone lower-bound approximation.']
    complete=not ci.empty and len(ci)==48 and not joint.empty
    lines+=['','**Limits/failures:** residual EP9 alignment errors are NOT repaired by masking 27,992 flags; exact-depth timeouts/p-intervals; 1,000-draw full-grid null replaced by permitted 200-draw local test; bootstrap d* CIs are NOT exact-d CIs; chronological unknown dates placed last (see bootstrap_order CSVs); resampled short blocks make total roll-call counts variable.','Status: '+('complete with the listed limitations.' if complete else 'IN PROGRESS; pending entries are not results.')]
    (OUT/'summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    failures=[]
    for term in ['EP6','EP9']:
        for mode in ['missing','current']:
            g=obs[(obs.term==term)&(obs['mode']==mode)]
            failures.append(dict(term=term,mode=mode,item='observed_exact_depth',status='timeouts',count=sum(g.status=='timeout'),detail='Counts include all 3 copresence minima on the 41-point grid; full bounds in observed.csv and depth_searches.csv'))
            failures.append(dict(term=term,mode=mode,item='null',status='budget_fallback',count=200,detail='Four thresholds near heuristic-selected main drop, not 1000 full-grid draws; 0.12s exact depth budget per graph; cached certified results reused'))
            failures.append(dict(term=term,mode=mode,item='bootstrap_depth',status='approximation',count=2000,detail='Exact Q/omega; uncapped fixed-order two-start greedy depth with backward monotone envelope; approximate-depth correlations/deltas, not exact-depth CIs'))
    failures.append(dict(term='EP9',mode='both',item='input_alignment',status='unrepaired_beyond_requested_mask',count=66001,detail='Packed-vs-physical Yes disagreements; coordinate-correct matrices retained in input_EP9.npz for audit, NOT substituted silently'))
    failures.append(dict(term='EP9',mode='both',item='party_metadata',status='legacy_labels_retained',count=12,detail='20 rows omit stored metadata cells; 12 cached groups incorrectly NI; party-aware null and Q retain legacy labels for controlled cache-based sensitivity, not fully coordinate-repaired analysis; see metadata_alignment_audit.csv'))
    pd.DataFrame(failures).to_csv(OUT/'failures_and_approximations.csv',index=False)
    print('SUMMARY', 'complete' if complete else 'pending',flush=True)
if __name__=='__main__':main()
