"""Independent small-matrix trajectory diagnostics on immutable contributor rows.

Usage: python3 trajectory_replay.py RUN_DIRECTORY OUTPUT.json
This is report reproduction code, not a search or certificate implementation.
"""
import json
import math
from pathlib import Path
import sys

import numpy as np


def read(p):
    return json.loads(p.read_text())


def spectrum(x):
    _, s, v = np.linalg.svd(x, full_matrices=False)
    threshold = np.finfo(float).eps*max(x.shape)*s[0]
    return {'rank': int(np.sum(s > threshold)), 'singular_values': s.tolist(),
            'threshold': float(threshold), 'column_norms': np.linalg.norm(x, axis=0).tolist(),
            'weakest_direction': v[-1].tolist()}


def diagnose(support, y, point):
    beta, b = np.array(point['beta']), np.array(point['b'])
    n, m = len(y), len(b)
    x, dx = np.zeros((n, 2*m)), np.zeros((n, 2*m))
    for a in range(m):
        cells = support[support['atom'] == a]; ids = cells['row'].astype(int); r2 = cells['square']
        r = np.sqrt(r2); exponential = np.exp(-r2/(2*b[a]**2)); center = math.sqrt(2/math.pi)/b[a]
        g = (2*math.pi*b[a]**2)**(-1.5)*exponential
        x[ids, 2*a] = g
        x[ids, 2*a+1] = [center if radius < 1e-5 else math.erf(radius/b[a]/math.sqrt(2))/radius for radius in r]
        dx[ids, 2*a] = g*(r2/b[a]**2-3)
        dx[ids, 2*a+1] = -center*np.where(r2 < 1e-10, 1., exponential)
    norms = np.linalg.norm(x, axis=0); free = np.flatnonzero((np.arange(2*m) % 2 == 1) | (beta > 0))
    if np.any(norms == 0): return {'available': False, 'reason': 'zero design column'}
    z = x[:, free]/norms[free]; q, r = np.linalg.qr(z, mode='reduced')
    residual = x@beta-y; scale = max(1., np.linalg.norm(y))
    raw = (dx*beta).reshape(n, m, 2).sum(axis=2)
    projected = (raw-q@(q.T@raw))/scale
    t = np.zeros((len(free), m))
    for j, col in enumerate(free): t[j, col//2] = dx[:, col]@residual/norms[col]
    jacobian = projected-q@np.linalg.solve(r.T, t)/scale
    gradient = (x.T@residual)/norms/scale; u = norms*beta/scale
    tangent = u-gradient; tangent[::2] = np.maximum(0., tangent[::2])
    return {'available': True, 'design': spectrum(x/norms), 'projected_width': spectrum(projected),
            'profile_jacobian': spectrum(jacobian), 'active_atoms': np.flatnonzero(beta[::2] == 0).tolist(),
            'projected_kkt': float(np.max(np.abs(u-tangent))), 'b_gradient': (raw.T@residual/scale/scale).tolist(),
            'maximum_cancellation_ratio': float(np.max((np.abs(x)@np.abs(beta))/np.maximum(1., np.abs(x@beta)))),
            'residual_norm': float(np.linalg.norm(residual))}


def main():
    root = Path(sys.argv[1])/'datasets/near-0.02'
    support = np.genfromtxt(root/'contributors.csv', delimiter=',', names=True)
    voxels = np.genfromtxt(root/'voxels.csv', delimiter=',', names=True)
    out = {'method': 'independent dense QR/SVD and scalar kernels from frozen contributor CSR; no truth inputs', 'cases': []}
    for variant in ('legacy', 'guarded', 'guarded-log'):
        for precision in ('double', 'float32'):
            p = root/variant/'fits'/f'narrower-{precision}.json'
            if not p.exists(): continue
            fit = read(p); y = voxels['reference_double' if precision == 'double' else 'observed']
            trials = []
            for trial in fit['trials']:
                if len(trial.get('beta', [])) != 24: continue
                result = diagnose(support, y, trial)
                trials.append({'evaluation': trial['evaluation'], 'accepted': trial['accepted'], 'primary_valid': trial['valid'],
                               'b1': trial['b'][0], 'replay': result, 'trust': trial.get('trust')})
            out['cases'].append({'variant': variant, 'precision': precision, 'trials': trials})
    Path(sys.argv[2]).write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__': main()
