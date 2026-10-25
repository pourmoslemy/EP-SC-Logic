"""Quantify approximation error on the observed input, never replace exact d."""
from common import *
from bootstrap_core import curve
def main():
    z,a,m=inputs();num,den=m.counts(a);w=m.omega(num,den,20)
    matrix=np.zeros((m.n,m.n));matrix[m.u,m.v]=w;matrix[m.v,m.u]=w
    order=np.argsort(-matrix.sum(1),kind='stable');q,d=curve(m,w,order)
    observed=pd.read_csv(OUT/'observed.csv');observed=observed[observed.minimum_copresent==20].sort_values('tau')
    assert np.allclose(q,observed.Q,atol=1e-12)
    assert np.all(d<=observed.d_upper.to_numpy())
    frame=observed[['tau','d_max','d_lower','d_upper','status']].copy()
    frame['bootstrap_depth_greedy']=d
    frame['exact_minus_greedy']=frame.d_max-d
    frame['delta_bootstrap_depth_greedy']=np.r_[np.nan,np.diff(d)]
    frame['greedy_is_approximation']=True
    save('bootstrap_depth_approximation_audit.csv',frame)
    step=pd.read_csv(OUT/'main_steps.csv');step=step[np.isclose(step.mesh_step,.02)].iloc[0]
    ix=int(np.flatnonzero(np.isclose(TAUS,step.tau_after))[0])
    point=pd.read_csv(OUT/'bootstrap_point.csv').iloc[0]
    assert point.delta_d_approx==int(d[ix]-d[ix-1])
    log('bootstrap-approximation-audit',f'exact selected delta_d={step.delta_d}; greedy point delta_d={point.delta_d_approx}; exact-depth CIs NOT inferred')
if __name__=='__main__':main()
