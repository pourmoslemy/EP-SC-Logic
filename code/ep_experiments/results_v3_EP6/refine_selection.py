"""Prefer certified exact step drops to greedy differences after stage 1."""
from common import *
def main():
    meshes=pd.read_csv(OUT/'mesh.csv');steps=[]
    for step,g in meshes.groupby('mesh_step',sort=False):
        chosen=selected(g);chosen.update(term='EP6',mode='fixed',mesh_step=step);steps.append(chosen)
    old=pd.read_csv(OUT/'main_steps.csv')
    save('main_steps_initial_greedy.csv',old)
    save('main_steps.csv',pd.DataFrame(steps),replace=True)
    save('selection_methods.csv',pd.DataFrame([dict(stage='initial_budget_placement',method='Preliminary greedy drop before exact searches; expanded window, plus fixed 0.38-0.52'),dict(stage='final_selection',method='Largest exact consecutive depth drop among certified steps; if competitor lower bounds permit a larger drop, flag global maximum unverified; greedy fallback only if no certified step')]))
    # The new main window must already have received the requested 10s budget.
    chosen=pd.DataFrame(steps);r=chosen[np.isclose(chosen.mesh_step,.02)].iloc[0]
    obs=pd.read_csv(OUT/'observed.csv');mask=(obs.tau>=r.tau_before-.02000001)&(obs.tau<=r.tau_after+.02000001)
    assert (obs.loc[mask,'requested_budget_seconds']==10).all()
    log('selection-refinement',f'final main step {r.tau_before:.2f}->{r.tau_after:.2f}; delta_d={r.delta_d}, delta_Q={r.delta_Q:.8f}; {r.selection}')
if __name__=='__main__':main()
