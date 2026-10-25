import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import numpy as np,pandas as pd,networkx as nx
from analysis_core import OUT,TAUS,Metrics,Depth
from bootstrap_core import curve
def main():
    rng=np.random.default_rng(42);records=[];solver=Depth()
    for n in [1,5,15,28]:
        for probability in [0,.1,.4,.8,1]:
            a=np.triu(rng.random((n,n))<probability,1);a|=a.T
            r=solver.solve(a,5);truth=max(map(len,nx.find_cliques(nx.from_numpy_array(a))))-1
            assert r['status']=='exact' and r['d_max']==truth
            records.append(dict(check='Exact pivoted BK vs networkx',n=n,value=probability,passed=True))
    a=rng.integers(-1,4,size=(17,123),dtype=np.int8);m=Metrics(np.arange(17)%3)
    blocks=[np.arange(start,min(start+20,123)) for start in range(0,123,20)]
    nums=[];dens=[]
    for block in blocks:
        num,den=m.counts(a[:,block]);nums.append(num);dens.append(den)
    for rep in range(5):
        chosen=rng.integers(0,len(blocks),len(blocks));weights=np.bincount(chosen,minlength=len(blocks))
        num,den=m.counts(a[:,np.concatenate([blocks[i] for i in chosen])])
        assert np.array_equal(weights@np.array(nums),num) and np.array_equal(weights@np.array(dens),den)
        records.append(dict(check='Weighted block Grams vs full expanded resample',n=17,value=rep,passed=True))
    num,den=m.counts(a);w=m.omega(num,den,20);q,d=curve(m,w,np.arange(m.n));matrix=np.zeros((m.n,m.n));matrix[m.u,m.v]=w;matrix[m.v,m.u]=w
    for i,tau in enumerate(TAUS):
        weighted=np.where(matrix>=tau,matrix,0.);np.fill_diagonal(weighted,0);volume=weighted.sum()/2
        truth=0.
        if volume:
            degrees=weighted.sum(1)
            for group in np.unique(m.labels):
                ix=np.flatnonzero(m.labels==group);truth+=weighted[np.ix_(ix,ix)].sum()/(2*volume)-(degrees[ix].sum()/(2*volume))**2
        exact=max(map(len,nx.find_cliques(nx.from_numpy_array(matrix>=tau))))-1
        assert np.isclose(q[i],truth,atol=1e-12) and d[i]<=exact
    assert np.all(np.diff(d)<=0)
    records.append(dict(check='Q histogram vs direct formula and bootstrap depth lower-bound/envelope',n=17,value=41,passed=True))
    pd.DataFrame(records).to_csv(OUT/'unit_validation.csv',index=False);print('ALL METHOD UNIT CHECKS PASSED',flush=True)
if __name__=='__main__':main()
