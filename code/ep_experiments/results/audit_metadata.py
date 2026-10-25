import sys,zipfile,importlib.util
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np,pandas as pd
sys.dont_write_bytecode=True
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
sp=importlib.util.spec_from_file_location('original',ROOT/'ep_topological_pipeline.py');p=importlib.util.module_from_spec(sp);sp.loader.exec_module(p)
def col(ref):
    n=0
    for c in ref:
        if not c.isalpha():break
        n=n*26+ord(c)-64
    return n-1
def main():
    records=[]
    for term,filename in [('EP6','EP6_RCVs_2022_06_13.xlsx'),('EP9','EP9_RCVs_2022_06_22.xlsx')]:
        cached=np.load(OUT/f'input_{term}.npz');index={int(v):i for i,v in enumerate(cached['members'])}
        with zipfile.ZipFile(ROOT.parent/'VoteWatch-EP-voting-data_2004-2022'/filename) as z:
            strings=p.shared_strings(z)
            for _,e in ET.iterparse(z.open('xl/worksheets/sheet1.xml'),events=('end',)):
                if e.tag!=p.NS+'row':continue
                rn=int(e.attrib['r']);cells=list(e.findall(p.NS+'c'));e.clear()
                if rn<=(3 if term=='EP6' else 1):continue
                first={col(c.attrib['r']):p.value(c,strings) for c in cells[:12]}
                if first.get(0) in (None,''):continue
                mid=int(float(first[0]));i=index[mid];missing=sorted(set(range(10))-set(first))
                physical=p.group_name(str(first.get(7) or ''));legacy=str(cached['groups'][i])
                if missing or physical!=legacy:
                    records.append(dict(term=term,member_id=mid,worksheet_row=rn,missing_metadata_columns_zero_based=','.join(map(str,missing)),cached_group=legacy,coordinate_group=physical,group_disagreement=physical!=legacy))
    frame=pd.DataFrame(records);frame.to_csv(OUT/'metadata_alignment_audit.csv',index=False)
    audit=pd.read_csv(OUT/'overlap_audit.csv')
    for term in ['EP6','EP9']:
        mask=audit.term==term;g=frame[frame.term==term]
        audit.loc[mask,'rows_with_missing_metadata_cells']=len(g)
        audit.loc[mask,'group_label_disagreements']=sum(g.group_disagreement)
    audit.to_csv(OUT/'overlap_audit.csv',index=False)
    print(frame.groupby('term').group_disagreement.agg(['size','sum']).to_string(),flush=True)
if __name__=='__main__':main()
