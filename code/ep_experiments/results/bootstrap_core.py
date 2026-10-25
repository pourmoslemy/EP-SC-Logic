"""Full-pair omega/Q per resample; declared uncapped greedy depth surrogate."""
import numpy as np
from analysis_core import TAUS,bits_for
def curve(m,w,order):
    nt=len(TAUS);idx=np.searchsorted(TAUS,w,side='right')-1;keep=idx>=0
    b=idx[keep];ww=w[keep];gu=m.gu[keep];gv=m.gv[keep]
    reverse=lambda a:np.cumsum(a[...,::-1],axis=-1)[...,::-1]
    volume=reverse(np.bincount(b,weights=ww,minlength=nt));same=gu==gv
    inside=reverse(np.bincount(gu[same]*nt+b[same],weights=ww[same],minlength=m.ng*nt).reshape(m.ng,nt))
    degree=reverse((np.bincount(gu*nt+b,weights=ww,minlength=m.ng*nt)+np.bincount(gv*nt+b,weights=ww,minlength=m.ng*nt)).reshape(m.ng,nt))
    q=np.zeros(nt);ok=volume>0;q[ok]=(inside[:,ok]/volume[ok]-(degree[:,ok]/(2*volume[ok]))**2).sum(0)
    matrix=np.zeros((m.n,m.n));matrix[m.u,m.v]=w;matrix[m.v,m.u]=w;matrix=matrix[np.ix_(order,order)]
    ds=[]
    for tau in TAUS:
        adj=matrix>=tau;np.fill_diagonal(adj,False);degrees=adj.sum(1)
        if not degrees.any():ds.append(0);continue
        bits=bits_for(adj);best=1
        # Fixed observed-strength vertex order, two starts, no dimension cap.
        for start in [int(np.argmax(degrees)),int(np.flatnonzero(degrees)[0])]:
            size=1;p=bits[start]
            while p:
                v=(p&-p).bit_length()-1;size+=1;p &= bits[v]
            best=max(best,size)
        ds.append(best-1)
    # A clique found at a larger tau is also valid at every smaller tau.
    # Backward propagation strengthens the lower bound without inventing faces.
    return q,np.maximum.accumulate(np.array(ds)[::-1])[::-1]
def statistics(q,d,low_index,high_index):
    if np.std(q)==0 or np.std(d)==0:pearson=spearman=np.nan
    else:
        import pandas as pd
        pearson=float(np.corrcoef(q,d)[0,1]);spearman=float(np.corrcoef(pd.Series(q).rank(),pd.Series(d).rank())[0,1])
    return dict(pearson_approx_d=pearson,spearman_approx_d=spearman,delta_Q=float(q[high_index]-q[low_index]),delta_d_approx=float(d[high_index]-d[low_index]))
