from common import *
class ScheduledDepth:
    def __init__(self,solver,taus,window):self.solver=solver;self.taus=list(taus);self.window=window;self.i=0
    def solve(self,adj,budget=.5):
        tau=self.taus[self.i];self.i+=1
        use=10.0 if .38-1e-8<=tau<=.52+1e-8 or self.window[0]-1e-8<=tau<=self.window[1]+1e-8 else .5
        r=self.solver.solve(adj,use)
        if use==10:log('observed-depth',f'tau={tau:.2f}, budget={use}s, {r["status"]}, bounds=[{r["d_lower"]},{r["d_upper"]}]')
        return r
def main():
    assert not (OUT/'observed.csv').exists(),'Refusing to overwrite observed output'
    fingerprints();z,a,m=inputs();num,den=m.counts(a);solver=Depth()
    w=m.omega(num,den,20)
    preliminary=m.trajectory(w,TAUS,approximate=True);choice=selected(preliminary)
    window=(choice['tau_before']-.02,choice['tau_after']+.02)
    log('observed',f'seed={SEED}; selected candidate {choice["tau_before"]:.2f}->{choice["tau_after"]:.2f}; extended budget window={window}')
    frames=[];meshes=[];steps=[];density=[]
    for minimum in [20,10,50]:
        taus=np.round(np.arange(.1,.9001,.01),2) if minimum==20 else TAUS
        f=m.trajectory(m.omega(num,den,minimum),taus,depth=ScheduledDepth(solver,taus,window))
        f['term']='EP6';f['mode']='fixed';f['minimum_copresent']=minimum;f=delta(f)
        f['requested_budget_seconds']=[10.0 if .38-1e-8<=t<=.52+1e-8 or window[0]-1e-8<=t<=window[1]+1e-8 else .5 for t in f.tau]
        common=delta(f[f.tau.isin(TAUS)]);frames.append(common)
        r=common[np.isclose(common.tau,.5)].iloc[0]
        density.append(dict(term='EP6',minimum_copresent=minimum,tau=.5,edges=r.edges,density=r.density,active_nodes=r.active_nodes,Q=r.Q))
        if minimum==20:
            for step in [.01,.02,.05]:
                target=np.round(np.arange(.1,.9001,step),2);mesh=delta(f[f.tau.isin(target)]);mesh['mesh_step']=step;meshes.append(mesh)
                chosen=selected(mesh);chosen.update(term='EP6',mode='fixed',mesh_step=step);steps.append(chosen)
        solver.save();log('observed',f'min co-presence {minimum}: exact={sum(f.status=="exact")}/{len(f)}')
    obs=pd.concat(frames,ignore_index=True);mesh=pd.concat(meshes,ignore_index=True)
    save('observed.csv',obs);save('mesh.csv',mesh);save('main_steps.csv',pd.DataFrame(steps));save('density_at_050.csv',pd.DataFrame(density))
    six=np.load(OLD/'input_EP6.npz');ep6=int(np.count_nonzero((six['current']==1)!=(six['coordinate']==1)))
    stored=z['stored_vote_cells'];present=z['present']
    audit=[dict(metric='EP6 packed_vs_physical_Yes_disagreements',value=ep6,definition='Current packed Yes vs independently reconstructed physical coordinates; same member/vote IDs'),dict(metric='EP6_stored_vote_cells',value=int(stored.sum()),definition='Physical vote cell has a nonempty raw value; may encode absent/not-voting, NOT presence'),dict(metric='EP6_present_cells',value=int(present.sum()),definition='Raw physical code 1 Yes, 2 No or 3 Abstain only'),dict(metric='EP6_stored_but_not_present',value=int(np.sum(stored&~present)),definition='Stored raw code outside 1,2,3; not in co-presence denominator'),dict(metric='EP6_present_without_stored_cell',value=int(np.sum(present&~stored)),definition='Must be zero'),dict(metric='UNKNOWN_members',value=int(np.sum(z['groups']=='UNKNOWN')),definition='Single separate UNKNOWN group, not NI; member-specific counts are recorded in party_groups.csv')]
    save('input_presence_audit_stage1.csv',pd.DataFrame(audit))
    save('party_groups.csv',pd.DataFrame(dict(member_id=z['members'],group=z['groups'],group_source=z['group_source'],present_votes=present.sum(1))))
    save('methods.csv',pd.DataFrame([dict(key=k,value=v) for k,v in dict(seed=42,omega='co-Yes / co-present (Yes,No,Abstain), minimum=20; observed sensitivities 10,50',observed_depth_budget='10s at 0.38-0.52 and +/-0.02 of selected .02-grid step; 0.5s elsewhere',selection='Preliminary greedy budget placement; final largest certified exact step, globally certified only if competing delta lower bounds exclude a larger drop; see selection_methods.csv',null_draws=1000,null_seed=42,null_depth_budget_seconds=.12,null_depth_budget_scope='Observed enhanced budget is step 1; null retains previous 0.12s bounded exact search per graph, no clique cap; p_d bounds handle ambiguity',null_unknown='0 UNKNOWN members retained together as own group, distinct from NI',presence='states in {1,2,3}; NEVER stored_vote_cells',bootstrap_reps=2000,bootstrap_lengths='20,50,100',bootstrap_depth='Approximate uncapped fixed-order greedy lower-bound envelope; not exact depth').items()]))
    fingerprints(verify=True);log('stage-1','complete')
if __name__=='__main__':main()
