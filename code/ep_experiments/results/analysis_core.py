"""Exact uncapped pivoted Bron--Kerbosch, timeouts, and common metrics."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='2'
import time,sys,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
sys.setrecursionlimit(6000)
TAUS=np.round(np.arange(.10,.9001,.02),2)
OUT=Path(__file__).resolve().parent
def bits_for(adj):return [int.from_bytes(r.tobytes(),'little') for r in np.packbits(adj,axis=1,bitorder='little')]
def greedy_depth(adj,starts=1):
    degree=adj.sum(1);best=1
    for order in [np.argsort(-degree),np.argsort(degree)]:
        bits=[int.from_bytes(r.tobytes(),'little') for r in np.packbits(adj[:,order],axis=1,bitorder='little')]
        for start in order[:starts]:
            size=1;p=bits[int(start)]
            while p:
                v=(p&-p).bit_length()-1;size+=1;p &= bits[int(order[v])]
            best=max(best,size)
    return best-1
class Expired(Exception):pass
class Depth:
    def __init__(self):
        self.cache={};self.path=OUT/'depth_searches.csv'
        if self.path.exists():
            for r in pd.read_csv(self.path).to_dict('records'):self.cache[r['graph_hash']]=r
    def save(self):pd.DataFrame(self.cache.values()).to_csv(self.path,index=False)
    def solve(self,adj,budget=.5):
        key=hashlib.blake2b(np.packbits(adj).tobytes(),digest_size=16).hexdigest()
        if key in self.cache:
            r=self.cache[key]
            if r['status']=='exact' or r['budget_seconds']>=budget:return r.copy()
        begin=time.perf_counter();deadline=begin+budget
        bits=bits_for(adj);n=len(bits);best=greedy_depth(adj,starts=3)+1
        nodes=0
        def coloring(p):
            colors=0
            while p:
                colors+=1;available=p
                while available:
                    vbit=available&-available;v=vbit.bit_length()-1;p^=vbit
                    available&=~(vbit|bits[v])
            return colors
        root_upper=coloring((1<<n)-1)
        def bk(size,p,x):
            nonlocal best,nodes
            nodes+=1
            if time.perf_counter()>deadline:raise Expired
            if size+p.bit_count()<=best:return
            if size+coloring(p)<=best:return
            if not p:
                if not x:best=max(best,size)
                return
            union=p|x;u=-1;score=-1
            while union:
                bit=union&-union;union^=bit;v=bit.bit_length()-1
                val=(bits[v]&p).bit_count()
                if val>score:u=v;score=val
            candidates=p & ~bits[u] if u>=0 else p
            while candidates:
                bit=candidates&-candidates;candidates^=bit;v=bit.bit_length()-1
                npart=p&bits[v]
                if not npart and not (x&bits[v]):best=max(best,size+1)
                else:bk(size+1,npart,x&bits[v])
                p^=bit;x|=bit
                if size+p.bit_count()<=best:return
        status='exact'
        try:
            if best<root_upper:bk(0,(1<<n)-1,0)
        except Expired:status='timeout'
        result=dict(graph_hash=key,status=status,d_max=best-1 if status=='exact' else np.nan,d_lower=best-1,d_upper=best-1 if status=='exact' else root_upper-1,seconds=time.perf_counter()-begin,budget_seconds=budget,search_nodes=nodes)
        self.cache[key]=result.copy();return result
class Metrics:
    def __init__(self,groups):
        self.groups=np.asarray(groups);_,self.labels=np.unique(self.groups,return_inverse=True)
        self.n=len(groups);self.ng=int(max(self.labels))+1
        self.u,self.v=np.triu_indices(self.n,1);self.gu=self.labels[self.u];self.gv=self.labels[self.v]
    def counts(self,a):
        y=(a==1).astype(np.float32);present=(a>0).astype(np.float32)
        return (y@y.T)[self.u,self.v],(present@present.T)[self.u,self.v]
    def omega(self,num,den,minimum=20):
        w=np.zeros(len(num),np.float64);np.divide(num,den,out=w,where=den>=minimum);return w
    def trajectory(self,w,taus=TAUS,depth=None,budget=.5,approximate=False):
        taus=np.asarray(taus);nt=len(taus);binidx=np.searchsorted(taus,w,side='right')-1
        keep=binidx>=0;b=binidx[keep];ww=w[keep];gu=self.gu[keep];gv=self.gv[keep]
        reverse=lambda a:np.cumsum(a[...,::-1],axis=-1)[...,::-1]
        volume=reverse(np.bincount(b,weights=ww,minlength=nt))
        count=reverse(np.bincount(b,minlength=nt))
        inside=(gu==gv)
        iw=reverse(np.bincount(gu[inside]*nt+b[inside],weights=ww[inside],minlength=self.ng*nt).reshape(self.ng,nt))
        strengths=reverse((np.bincount(gu*nt+b,weights=ww,minlength=self.ng*nt)+np.bincount(gv*nt+b,weights=ww,minlength=self.ng*nt)).reshape(self.ng,nt))
        q=np.zeros(nt);valid=volume>0;q[valid]=(iw[:,valid]/volume[valid]-(strengths[:,valid]/(2*volume[valid]))**2).sum(0)
        degrees=reverse((np.bincount(self.u[keep]*nt+b,minlength=self.n*nt)+np.bincount(self.v[keep]*nt+b,minlength=self.n*nt)).reshape(self.n,nt))
        matrix=np.zeros((self.n,self.n));matrix[self.u,self.v]=w;matrix[self.v,self.u]=w
        rows=[]
        for i,tau in enumerate(taus):
            adj=matrix>=tau;np.fill_diagonal(adj,False)
            g=greedy_depth(adj)
            r=dict(tau=tau,edges=int(count[i]),density=count[i]/len(w),active_nodes=int(np.count_nonzero(degrees[:,i])),Q=float(q[i]),d_greedy=g)
            if approximate:r.update(d_max=np.nan,status='approx_uncapped_greedy',d_lower=g,d_upper=self.n-1)
            else:r.update(depth.solve(adj,budget))
            rows.append(r)
        frame=pd.DataFrame(rows);frame['delta_Q']=frame.Q.diff();frame['delta_d']=frame.d_max.diff();frame['delta_d_greedy']=frame.d_greedy.diff()
        return frame
