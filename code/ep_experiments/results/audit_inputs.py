"""Read-only audit of the two historical Excel parsing conventions."""
from pathlib import Path
import sys,importlib.util,pickle,zipfile,re,json
from collections import Counter
from xml.etree import ElementTree as ET
import numpy as np
import pandas as pd
sys.dont_write_bytecode=True
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
s=importlib.util.spec_from_file_location('original',ROOT/'ep_topological_pipeline.py')
p=importlib.util.module_from_spec(s);s.loader.exec_module(p)
def column(ref):
    n=0
    for c in re.match('[A-Z]+',ref)[0]:n=26*n+ord(c)-64
    return n-1
def main():
    with (ROOT/'ep_aligned_cache.pkl').open('rb') as f:data=pickle.load(f)
    audit=[];examples=[]
    for term,d in data.items():
        dest=OUT/f'input_{term}.npz'
        if dest.exists():continue
        raw=ROOT.parent/'VoteWatch-EP-voting-data_2004-2022'/('EP6_RCVs_2022_06_13.xlsx' if term=='EP6' else 'EP9_RCVs_2022_06_22.xlsx')
        yes=d['x'].toarray().astype(bool);no=np.load(ROOT/f'ep_no_votes_{term}.npz')['no'][:,:yes.shape[1]]
        flags=yes&no;nv=yes.shape[1];header=2 if term=='EP6' else 1
        packed=[];physical=[];ids=[];shifted=0;physical_refs=None;voteids=[];raw_width=[]
        with zipfile.ZipFile(raw) as z:
            strings=p.shared_strings(z)
            for _,e in ET.iterparse(z.open('xl/worksheets/sheet1.xml'),events=('end',)):
                if e.tag!=p.NS+'row':continue
                rn=int(e.attrib['r']);cells=list(e.findall(p.NS+'c'))
                vals=[p.value(c,strings) for c in cells];coords=[column(c.attrib['r']) for c in cells];refs=[c.attrib['r'] for c in cells];e.clear()
                row=dict(zip(coords,vals))
                if rn==header:
                    physical_refs=[c for c in sorted(row) if c>=10 and row[c] not in (None,'')][:nv]
                    voteids=[int(float(row[c])) for c in physical_refs];continue
                if rn<=(3 if term=='EP6' else 1) or not vals or vals[0] in (None,''):continue
                ids.append(int(float(vals[0])));raw_width.append(len(vals)-10)
                rawcodes=np.array([int(float(v)) if v not in (None,'') else 0 for v in vals[10:10+nv]],np.int8)
                coordcodes=np.array([int(float(row[c])) if row.get(c) not in (None,'') else 0 for c in physical_refs],np.int8)
                packed.append(rawcodes);physical.append(coordcodes)
                diff=np.array(coords[10:10+nv])!=physical_refs;shifted+=int(diff.sum())
                i=len(ids)-1
                for j in np.flatnonzero(flags[i])[:max(0,30-sum(x['term']==term for x in examples))]:
                    examples.append(dict(term=term,member_id=ids[-1],cached_vote_id=d['vote_ids'][j],yes_source_cell=refs[10+j],no_source_column_index=physical_refs[j],packed_code=int(rawcodes[j]),coordinate_code=int(coordcodes[j])))
        a=np.stack(packed);c=np.stack(physical)
        assert ids==list(d['members']);assert np.array_equal(a==1,yes)
        assert np.array_equal(c==2,no),'Historical No cache has another mismatch'
        modes=dict(current=np.where(np.isin(a,[1,2,3]),a,0).astype(np.int8))
        modes['missing']=modes['current'].copy();modes['missing'][flags]=-1
        np.savez_compressed(dest,**modes,coordinate=np.where(np.isin(c,[1,2,3]),c,0).astype(np.int8),overlap=flags,groups=np.asarray(d['groups'],str),members=np.asarray(ids),vote_ids=np.asarray(d['vote_ids']))
        audit.append(dict(term=term,members=len(ids),votes=nv,overlap_cells=int(flags.sum()),duplicate_member_rows=len(ids)-len(set(ids)),duplicate_cached_vote_ids=nv-len(set(d['vote_ids'])),duplicate_physical_vote_ids=len(voteids)-len(set(voteids)),header_id_sequence_matches=voteids==d['vote_ids'],packed_coordinate_yes_disagreements=int(np.count_nonzero((a==1)!=(c==1))),packed_coordinate_state_disagreements=int(np.count_nonzero(a!=c)),shifted_cell_positions=shifted,physical_yes_no_overlap=0,minimum_packed_row_width=min(raw_width),maximum_packed_row_width=max(raw_width),cause='Yes parser discards XML cell coordinates and packs stored cells; No parser indexes physical header columns. No merging or duplicate member/vote IDs is needed to produce the overlap.',analysis_basis='Existing packed Yes/presence coordinates retained for controlled missing-vs-current comparison; masking does not repair residual coordinate shifts.'))
        pd.DataFrame(audit).to_csv(OUT/'overlap_audit.csv',index=False);pd.DataFrame(examples).to_csv(OUT/'overlap_examples.csv',index=False)
        print('[INFO]',term,audit[-1],flush=True)
if __name__=='__main__':main()
