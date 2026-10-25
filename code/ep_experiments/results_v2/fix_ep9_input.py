"""Coordinate-indexed EP9 reconstruction; immutable old outputs, seed 42.

Excel physical columns are selected by their vote-ID header, never by their
position in the XML cell list. Missing stored cells never shift other cells.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='2'
from pathlib import Path
import csv,hashlib,io,sys,time,zipfile,importlib.util
from xml.etree import ElementTree as ET
import numpy as np
import pandas as pd
from openpyxl import load_workbook
sys.dont_write_bytecode=True
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;OLD=ROOT/'results'
RAW=ROOT.parent/'VoteWatch-EP-voting-data_2004-2022/EP9_RCVs_2022_06_22.xlsx'
EXPORT=ROOT.parent/'export.zip'
SEED=42
spec=importlib.util.spec_from_file_location('legacy_readonly',ROOT/'ep_topological_pipeline.py')
legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
def col(ref):
    n=0
    for ch in ref:
        if not ch.isalpha():break
        n=n*26+ord(ch)-64
    return n-1
def save(name,frame):
    with (OUT/name).open('x',newline='',encoding='utf-8') as f:frame.to_csv(f,index=False)
def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    start=time.perf_counter();np.random.seed(SEED)
    assert not (OUT/'input_EP9_fixed.npz').exists(),'Refusing to overwrite existing fixed input'
    paths=sorted(p for p in OLD.iterdir() if p.is_file())+[RAW,EXPORT,ROOT/'ep_aligned_cache.pkl',ROOT/'ep_topological_pipeline.py']
    hashes={str(p):digest(p) for p in paths}
    prior=np.load(OLD/'input_EP9.npz');members=prior['members'];oldnv=len(prior['vote_ids'])
    with zipfile.ZipFile(RAW) as z:
        strings=legacy.shared_strings(z)
        first=next(e for _,e in ET.iterparse(z.open('xl/worksheets/sheet1.xml'),events=('end',)) if e.tag==legacy.NS+'row')
        header_pairs=[(col(c.attrib['r']),int(float(legacy.value(c,strings)))) for c in first if col(c.attrib['r'])>=10 and legacy.value(c,strings) not in (None,'')]
    vote_ids=np.asarray([v for _,v in sorted(header_pairs)])
    assert np.array_equal(vote_ids[:oldnv],prior['vote_ids']),'Existing vote IDs not a physical prefix'
    member_index={int(v):i for i,v in enumerate(members)};nv=len(vote_ids)
    states=np.zeros((len(members),nv),np.int8);observed=np.zeros_like(states,bool)
    groups=np.full(len(members),'UNKNOWN',dtype='<U40');rawgroups=np.full(len(members),'',dtype='<U100')
    countries=np.full(len(members),'',dtype='<U80');names=np.full(len(members),'',dtype='<U120')
    group_source=np.full(len(members),'physical_H',dtype='<U100')
    metadata=[];missing_cells=[];limitations=[];seen=set();code_counts={};columns=None
    with zipfile.ZipFile(RAW) as z:
        strings=legacy.shared_strings(z)
        for _,element in ET.iterparse(z.open('xl/worksheets/sheet1.xml'),events=('end',)):
            if element.tag!=legacy.NS+'row':continue
            rn=int(element.attrib['r']);cells=list(element.findall(legacy.NS+'c'))
            row={col(c.attrib['r']):legacy.value(c,strings) for c in cells};element.clear()
            if rn==1:
                pairs=[(c,int(float(v))) for c,v in sorted(row.items()) if c>=10 and v not in (None,'')]
                header=dict((v,c) for c,v in pairs)
                assert len(header)==len(pairs),'Duplicate physical vote IDs'
                assert len(pairs)==nv and [v for _,v in pairs]==list(vote_ids),'Vote header/cache ID mismatch'
                columns=np.array([header[int(v)] for v in vote_ids]);continue
            if row.get(0) in (None,''):continue
            mid=int(float(row[0]));assert mid in member_index and mid not in seen
            seen.add(mid);i=member_index[mid];name=str(row.get(3) or '')
            rawgroup=str(row.get(7) or '').strip();rawgroups[i]=rawgroup
            names[i]=name;countries[i]=str(row.get(5) or '')
            if rawgroup:
                groups[i]=legacy.group_name(rawgroup)
                if groups[i]=='NI' and rawgroup.upper() not in ['NI','NA','NON-ATTACHED','NON ATTACHED']:
                    limitations.append(dict(item='unrecognized_group',member_id=mid,status='unresolved',detail=rawgroup))
                    groups[i]='UNKNOWN';group_source[i]='unrecognized_physical_H'
            else:group_source[i]='missing_physical_H'
            missing=sorted(set(range(10))-set(row))
            if missing:missing_cells.append(dict(member_id=mid,worksheet_row=rn,missing_columns_zero_based=','.join(map(str,missing)),raw_group=rawgroup))
            for j,c in enumerate(columns):
                value=row.get(int(c))
                if value in (None,''):continue
                code=int(float(value));code_counts[code]=code_counts.get(code,0)+1
                observed[i,j]=True
                if code in (1,2,3):states[i,j]=code
                elif code not in (0,4,5,6):
                    limitations.append(dict(item='unrecognized_vote_code',member_id=mid,status='unresolved',detail=f'vote_id={vote_ids[j]}; raw_code={code}'))
            metadata.append(dict(member_id=mid,worksheet_row=rn,name=name,raw_group=rawgroup,national_party=str(row.get(6) or ''),missing_metadata_count=len(missing)))
            if len(seen)%100==0:print(f'[INFO] physical EP9 rows {len(seen)}/{len(members)}; seed={SEED}',flush=True)
    assert seen==set(map(int,members))
    # Only missing H cells use independent, term-specific local memberships.
    with zipfile.ZipFile(EXPORT) as z:
        memberships=list(csv.DictReader(io.TextIOWrapper(z.open('group_memberships.csv'),encoding='utf-8')))
    for i in np.flatnonzero(groups=='UNKNOWN'):
        mid=int(members[i]);matches=[r for r in memberships if int(r['member_id'])==mid and r['term']=='9']
        labels={legacy.group_name(r['group_code']) for r in matches}
        if len(labels)==1 and 'NI' not in labels:
            groups[i]=labels.pop();group_source[i]='export.zip/group_memberships.csv;term=9'
        else:
            limitations.append(dict(item='missing_EP9_group',member_id=mid,status='explicit_UNKNOWN',detail=f'Physical H blank; no unambiguous term-9 membership; present_votes={int(np.sum(states[i]>0))}; no NI imputation'))
    # Independent old physical reconstruction plus openpyxl coordinate checks.
    checks=[]
    def check(label,condition,count=0):
        checks.append(dict(check=label,passed=bool(condition),disagreements=int(count)))
        assert condition,label
    dif=int(np.count_nonzero(states[:,:oldnv]!=prior['coordinate']))
    check('shared vote IDs: fixed states vs independent previous physical reconstruction',dif==0,dif)
    yes=states==1;no=states==2;abstain=states==3;present=states>0
    check('shared vote IDs: fixed Yes vs independent physical Yes',np.array_equal(yes[:,:oldnv],prior['coordinate']==1),np.count_nonzero(yes[:,:oldnv]!=(prior['coordinate']==1)))
    check('mutually exclusive Yes No Abstain',not np.any((yes&no)|(yes&abstain)|(no&abstain)))
    wb=load_workbook(RAW,read_only=True,data_only=True);ws=wb.worksheets[0]
    rng=np.random.default_rng(SEED);samples=np.unique(np.r_[0,nv-1,rng.integers(0,nv,20)])
    byrow={r['worksheet_row']:member_index[r['member_id']] for r in metadata if r['worksheet_row']<=53}
    bad=0;tested=0
    for rn,values in enumerate(ws.iter_rows(min_row=2,max_row=53,values_only=True),2):
        if rn not in byrow:continue
        i=byrow[rn]
        for j in samples:
            v=values[int(columns[j])] if int(columns[j])<len(values) else None
            code=int(v) if v not in (None,'') else 0
            expected=code if code in (1,2,3) else 0
            bad+=int(states[i,j]!=expected);tested+=1
        bad+=int(rawgroups[i]!=str(values[7] or '').strip());tested+=1
    wb.close();check(f'independent openpyxl coordinates ({tested} comparisons)',bad==0,bad)
    changed=(prior['groups']!=groups)
    rows=pd.DataFrame(metadata);rows['old_group']=prior['groups'];rows['fixed_group']=groups;rows['group_changed']=changed;rows['group_source']=group_source;rows['present_votes']=present.sum(1)
    save('member_metadata_audit.csv',rows)
    save('missing_metadata_cells.csv',pd.DataFrame(missing_cells))
    save('vote_codes.csv',pd.DataFrame([dict(raw_code=k,count=v,state='Yes' if k==1 else 'No' if k==2 else 'Abstain' if k==3 else 'not_present') for k,v in sorted(code_counts.items())]))
    records=[]
    for label,a,gg in [('legacy_current',prior['current'],prior['groups']),('legacy_missing',prior['missing'],prior['groups']),('fixed_physical',states,groups)]:
        ref=states[:,:a.shape[1]]
        records.append(dict(version=label,members=len(members),votes=a.shape[1],yes_cells=int(np.sum(a==1)),no_cells=int(np.sum(a==2)),abstain_cells=int(np.sum(a==3)),present_cells=int(np.sum(a>0)),yes_vs_physical_disagreements=int(np.sum((a==1)!=(ref==1))),all_state_vs_physical_disagreements=int(np.sum(a!=ref)),group_vs_fixed_disagreements=int(np.sum(gg!=groups)),yes_no_overlap=int(np.sum(prior['overlap'])) if label=='legacy_current' else 0,unknown_groups=int(np.sum(gg=='UNKNOWN')),missing_physical_vote_columns=nv-a.shape[1]))
    save('input_before_after.csv',pd.DataFrame(records));save('input_validation.csv',pd.DataFrame(checks))
    save('restored_vote_columns.csv',pd.DataFrame([dict(vote_id=int(vote_ids[j]),physical_column_zero_based=int(columns[j]),yes_cells=int(yes[:,j].sum()),no_cells=int(no[:,j].sum()),abstain_cells=int(abstain[:,j].sum()),present_cells=int(present[:,j].sum())) for j in range(oldnv,nv)]))
    with (OUT/'input_EP9_fixed.npz').open('xb') as f:
        np.savez_compressed(f,states=states,yes=yes,no=no,abstain=abstain,present=present,stored_vote_cells=observed,groups=groups,raw_groups=rawgroups,group_source=group_source,members=members,vote_ids=vote_ids,physical_columns_zero_based=columns,countries=countries,names=names,seed=np.array(SEED),missing_group=groups=='UNKNOWN')
    save('limitations.csv',pd.DataFrame(limitations,columns=['item','member_id','status','detail']))
    fingerprints=[]
    for p in paths:
        after=digest(p);check('unchanged '+str(p),after==hashes[str(p)])
        fingerprints.append(dict(path=str(p),before_sha256=hashes[str(p)],after_sha256=after,unchanged=after==hashes[str(p)]))
    save('source_preservation.csv',pd.DataFrame(fingerprints))
    save('seeds_and_methods.csv',pd.DataFrame([dict(stage='input_fix',seed=SEED,method='Physical XML column references; vote-ID header mapping; raw H group labels; missing groups use term-9 local membership only; independent openpyxl samples',status='completed')]))
    print(pd.DataFrame(records).to_string(index=False),flush=True)
    print(rows[rows.group_changed][['member_id','old_group','fixed_group','group_source','present_votes']].to_string(index=False),flush=True)
    print(f'[INFO] input fix completed {time.perf_counter()-start:.1f}s; limitations={len(limitations)}',flush=True)
if __name__=='__main__':main()
