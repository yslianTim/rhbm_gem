"""Shared numerical comparisons used by independent validation campaigns."""
import math
import numpy as np


def sparse_parity(a, b):
    checks = {}
    for role in ('primary', 'reference'):
        if role not in a or role not in b:
            checks[role] = dict(available=False)
            continue
        x, y = a[role], b[role]
        same = all(x.get(k) == y.get(k) for k in ('valid', 'reason', 'feasible', 'free_rank', 'kkt_passed', 'active_atoms'))
        available = x['valid'] and y['valid'] and len(x['beta']) == len(y['beta']) and all(i is not None for i in x['beta']+y['beta'])
        error = max((abs(i-j)/(1+max(abs(i),abs(j))) for i,j in zip(x['beta'],y['beta'])), default=0.) if available else None
        objective = abs(x['objective']-y['objective']) if available else None
        normalized = .5*abs(x['relative_residual']**2-y['relative_residual']**2) if available else None
        gradient = available and len(x['b_gradient'])==len(y['b_gradient']) and all(
            i is not None and z is not None and abs(i-z)<=1e-13+2e-9*abs(i) for i,z in zip(x['b_gradient'],y['b_gradient']))
        checks[role] = dict(available=available, same_status_rank_face=same, coefficient_error=error,
                            objective_difference=objective, normalized_objective_difference=normalized, gradient_passed=gradient,
                            passed=same and available and error<=1e-10 and normalized<=1e-12 and gradient)
    return checks


def scaled(a, b):
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        return None
    return float(np.max(np.abs(x-y)/(1+np.maximum(np.abs(x), np.abs(y))), initial=0))


def svd_parity(a, b):
    x, y = np.asarray(a['singular_values'], dtype=float), np.asarray(b['singular_values'], dtype=float)
    same = bool(a['valid'] == b['valid'] and a['rank'] == b['rank'] and x.shape == y.shape
                and np.isfinite(x).all() and np.isfinite(y).all())
    maximum = float(x[0]) if x.size else 0.
    error = float(np.max(np.abs(x-y), initial=0)/max(maximum, np.finfo(float).tiny)) if same else None
    solution = scaled(a['solution'], b['solution'])
    threshold = abs(a['threshold']-b['threshold']) <= 1e-10*max(abs(a['threshold']), np.finfo(float).tiny)
    # Retain the actual weak-end values; a global norm alone does not verify rank.
    near = [dict(index=k, baseline=float(s), candidate=float(y[k]), threshold=a['threshold'])
            for k,s in enumerate(x) if same and (s <= 2*a['threshold'] or k >= len(x)-3)]
    return dict(same_rank_status=same, spectrum_error=error, solution_error=solution,
                threshold_agrees=threshold, weak_values=near,
                passed=same and a['valid'] and error <= 1e-10 and threshold and solution is not None and solution <= 1e-10)


def audit_parity(a, b):
    endpoints = sparse_parity(a, b)
    checks = [svd_parity(x, y) for x,y in zip(a.get('svd_records', []), b.get('svd_records', []))]
    same_calls = (len(a.get('svd_records', [])) == len(b.get('svd_records', [])) > 0
                  and all(all(x.get(k)==y.get(k) for k in ('role','rows','columns','relative_threshold','absolute_override'))
                          for x,y in zip(a.get('svd_records', []), b.get('svd_records', []))))
    left, right = a.get('derivative_audit', {}), b.get('derivative_audit', {})
    arrays = {k: scaled(left.get(k, []), right.get(k, []))
              for k in ('coefficients', 'correction', 'projected', 'jacobian', 'response')}
    passed = (all(c['passed'] for c in endpoints.values()) and a['trust']['passed'] and b['trust']['passed']
              and a['initial_b'] == b['initial_b'] and same_calls and all(c['passed'] for c in checks)
              and left.get('valid') is True and right.get('valid') is True
              and a['derivative_reason'] == b['derivative_reason']
              and all(k in left and k in right and len(left[k]) > 0 for k in arrays)
              and all(e is not None and e <= 1e-10 for e in arrays.values()))
    return dict(endpoints=endpoints, spectra=checks, derivative_errors=arrays, passed=passed)


def finite(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


def valid_state(report):
    state=report.get('returned_state')
    if not isinstance(state,dict): return False
    beta,b,active=(state.get(k) for k in ('beta','b','active_atoms'))
    if not isinstance(b,list) or not b or not all(finite(x) and x>0 for x in b): return False
    if not isinstance(beta,list) or len(beta)!=2*len(b) or not all(finite(x) for x in beta): return False
    if any(x<0 for x in beta[::2]): return False
    if not isinstance(active,list) or any(type(x) is not int or not 0<=x<len(b) for x in active): return False
    if len(set(active))!=len(active) or sorted(active)!=[i for i,x in enumerate(beta[::2]) if x==0]: return False
    scale=report.get('search',{}).get('residual_scale')
    return finite(state.get('objective')) and state['objective']>=0 and finite(scale) and scale>0


def scientific_parity(left,right):
    """Trajectories may differ; returned evidence and qualified endpoints may not."""
    a,b=left.get('returned_assessment'),right.get('returned_assessment')
    if (a is None)!=(b is None): return dict(passed=False,reason='state-availability')
    if a is None:
        if left.get('returned_state') is not None or right.get('returned_state') is not None: return dict(passed=False,reason='uncertified-returned-state')
        a,b=left.get('search',{}).get('initial',{}),right.get('search',{}).get('initial',{})
        same=a.get('valid') is False and b.get('valid') is False and isinstance(a.get('reason'),str) and bool(a['reason']) and a.get('reason')==b.get('reason')
        return dict(passed=same,reason='equivalent-explicit-unavailability' if same else 'uncertified-states')
    if not valid_state(left) or not valid_state(right): return dict(passed=False,reason='invalid-returned-state')
    failures=[]
    if a.get('runtime_convergence') not in ('passed','failed','unavailable','not-run') or not a.get('runtime_failure'): failures.append('missing-runtime-evidence')
    for key in ('runtime_convergence','runtime_checks','runtime_failure'):
        if a.get(key)!=b.get(key): failures.append(key)
    if a.get('runtime_convergence')=='passed' and not isinstance(a.get('runtime_checks'),dict): failures.append('runtime_checks/missing')
    for key in ('design_spectrum','width_spectrum','profile_jacobian_spectrum'):
        if (key in a)!=(key in b): failures.append(key+'/availability')
        if a.get('runtime_convergence')=='passed' and key not in a: failures.append(key+'/missing')
        for record in (a,b):
            if key in record and (not isinstance(record[key],dict) or type(record[key].get('rank')) is not int or record[key]['rank']<0): failures.append(key+'/invalid')
        if isinstance(a.get(key,{}),dict) and isinstance(b.get(key,{}),dict) and a.get(key,{}).get('rank')!=b.get(key,{}).get('rank'): failures.append(key+'/rank')
    sa,sb=left['returned_state'],right['returned_state']
    if len(sa['b'])!=len(sb['b']): failures.append('state-dimensions')
    if sa.get('active_atoms')!=sb.get('active_atoms'): failures.append('active-face')
    if left.get('search_completed') and not right.get('search_completed'): failures.append('search-incomplete')
    oa=sa['objective']/left['search']['residual_scale']/left['search']['residual_scale']
    ob=sb['objective']/right['search']['residual_scale']/right['search']['residual_scale']
    if not finite(oa) or not finite(ob): return dict(passed=False,reason='invalid-normalized-objective')
    if a.get('runtime_convergence')=='passed':
        for key in ('beta','b'):
            x,y=sa.get(key,[]),sb.get(key,[])
            if len(x)!=len(y) or any(abs(u-w)>1e-10*(1+max(abs(u),abs(w))) for u,w in zip(x,y)): failures.append(key)
        if abs(oa-ob)>1e-12: failures.append('normalized-objective')
    elif ob>oa+1e-12: failures.append('normalized-objective-regression')
    return dict(passed=not failures,reason='same-evidence-and-qualified-endpoint' if a.get('runtime_convergence')=='passed' else 'unconverged-not-worse',differences=failures)


def vector_close(a,b,relative=1e-8,absolute=1e-12):
    if not isinstance(a,list) or not isinstance(b,list) or not a or len(a)!=len(b): return False
    if not all(finite(x) for x in a+b): return False
    return math.hypot(*(x-y for x,y in zip(a,b)))<=absolute+relative*math.hypot(*a)


def fixed_action_parity(a,b,kinds=('identity','diagonal','schwarz')):
    failures=[]
    for key in ('valid','reason'):
        if key not in a or a.get(key)!=b.get(key): failures.append(key)
    sa,sb=a.get('state_control',{}),b.get('state_control',{})
    for key in ('eta','beta','free_columns','scale','rank_rows'):
        if key not in sa or sa.get(key)!=sb.get(key): failures.append('state/'+key)
    for key in ('residual','gradient'):
        if not vector_close(sa.get(key),sb.get(key),relative=0,absolute=1e-12): failures.append('state/'+key)
    oa,ob=sa.get('objective'),sb.get('objective')
    if not finite(oa) or not finite(ob) or abs(oa-ob)>1e-12: failures.append('state/objective')
    ra,rb=a.get('rank',[]),b.get('rank',[])
    if len(ra)!=1 or len(rb)!=1: failures.append('rank/missing')
    else:
        for key in ('valid','rank','rows','columns','relative_threshold','absolute_override'):
            if key not in ra[0] or ra[0].get(key)!=rb[0].get(key): failures.append('rank/'+key)
        x,y=ra[0].get('singular_values'),rb[0].get('singular_values')
        if not isinstance(x,list) or not isinstance(y,list) or not x or len(x)!=len(y) or not all(finite(z) for z in x+y) or max(abs(u-w) for u,w in zip(x,y))>1e-10*max(x): failures.append('rank/spectrum')
        x,y=ra[0].get('threshold'),rb[0].get('threshold')
        if not finite(x) or not finite(y) or abs(x-y)>1e-10*max(abs(x),abs(y)): failures.append('rank/threshold')
    if a.get('valid') is not True or b.get('valid') is not True: failures.append('unavailable-operator')
    for key in ('apply','adjoint','normal'):
        if not vector_close(a.get(key),b.get(key)): failures.append(key)
    x,y=a.get('operator_gradient',[]),b.get('operator_gradient',[])
    if not x or len(x)!=len(y) or not all(finite(z) for z in x+y) or any(abs(u-w)>1e-13+2e-9*abs(u) for u,w in zip(x,y)): failures.append('operator-gradient')
    for result in (a,b):
        if result.get('normal_q_actions')!=(2 if result.get('mode')=='normal' else 6): failures.append('normal-q-actions')
        if any(result.get('work',{}).get(k)!=0 for k in ('reference_solves','derivative_preparations')): failures.append('forbidden-work')
    aa,bb=a.get('steps',[]),b.get('steps',[])
    if len(aa)!=len(kinds) or len(bb)!=len(kinds): failures.append('steps/missing')
    for index,kind in enumerate(kinds):
        if len(aa)<=index or len(bb)<=index: continue
        x,y=aa[index],bb[index]
        if x.get('kind')!=kind or y.get('kind')!=kind or x.get('valid') is not True or y.get('valid') is not True: failures.append(kind+'/unavailable'); continue
        if not vector_close(x.get('step'),y.get('step'),1e-10,1e-10): failures.append(kind+'/step')
        if not finite(x.get('predicted')) or not finite(y.get('predicted')) or abs(x['predicted']-y['predicted'])>1e-12: failures.append(kind+'/predicted')
        if any(not finite(r.get('true_residual')) or r['true_residual']>1e-10 for r in (x,y)): failures.append(kind+'/true-residual')
    return dict(passed=not failures,differences=failures)
