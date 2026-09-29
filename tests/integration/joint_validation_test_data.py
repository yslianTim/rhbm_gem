"""Small report values shared by historical validation unit tests."""


def fixed_result():
    return dict(valid=True,reason='full-profile-operator',mode='composed',normal_q_actions=6,
        state_control=dict(eta=[.2],beta=[1.,.1],free_columns=[0,1],scale=2.,rank_rows=10,
                           residual=[.1,.2],gradient=[.03],objective=.025),
        rank=[dict(valid=True,rank=2,rows=2,columns=2,relative_threshold=1e-12,absolute_override=-1,
                   threshold=1e-12,singular_values=[1.,.5])],
        apply=[.1,.2],adjoint=[.3],normal=[.4],operator_gradient=[.03],
        work=dict(reference_solves=0,derivative_preparations=0),
        steps=[dict(kind=k,valid=True,step=[.2],predicted=.001,true_residual=1e-11) for k in ('identity','diagonal','schwarz')])


def search_result():
    return dict(stage='complete',search_completed=True,search=dict(residual_scale=2),
        returned_assessment=dict(runtime_convergence='passed',runtime_failure='none',runtime_checks=dict(inner=True,local_correction=True),design_spectrum=dict(rank=2),width_spectrum=dict(rank=1),profile_jacobian_spectrum=dict(rank=1)),
        returned_state=dict(objective=.04,beta=[2.,.2],b=[.5],active_atoms=[]))
