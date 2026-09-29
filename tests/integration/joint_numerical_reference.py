"""Independent scalar replay of frozen Joint Component observations."""
import math
import numpy as np

from joint_runtime_support import require


def predict(basis, beta, n):
    out = np.zeros(n)
    for (ids, values), coefficient in zip(basis, beta): out[ids] += values*coefficient
    return out


def snapshot_replay(table, atom_count, y, endpoint):
    """Independent scalar kernels on immutable CSR memberships and distances."""
    beta, widths = np.asarray(endpoint["beta"]), np.asarray(endpoint["b"])
    basis, derivative = [], []
    for a, b in enumerate(widths):
        cells = table[table["atom"] == a]; ids = cells["row"].astype(int); square = cells["square"]
        radius = np.sqrt(square); exponent = np.exp(-square/(2*b*b))
        gaussian = (2*math.pi*b*b)**(-1.5)*exponent; center = math.sqrt(2/math.pi)/b
        charge = np.fromiter((center if r < 1e-5 else math.erf(r/b/math.sqrt(2))/r for r in radius), float)
        basis.extend(((ids, gaussian), (ids, charge)))
        derivative.append((ids, beta[2*a]*gaussian*(square/(b*b)-3)-beta[2*a+1]*center*np.where(square < 1e-10, 1., exponent)))
    require(len(widths) == atom_count, "Wrong endpoint atom count.")
    prediction = predict(basis, beta, len(y)); residual = prediction-y
    scale = max(1., np.linalg.norm(y)); norms = np.array([np.linalg.norm(v) for _, v in basis])
    require(np.all(norms > 0), "Zero design column.")
    u = norms*beta/scale; gradient = np.array([v@residual[ids] for ids, v in basis])/norms/scale
    projected = u-gradient; projected[::2] = np.maximum(0, projected[::2])
    return {"prediction": prediction, "residual": residual, "projected_kkt": float(np.max(np.abs(u-projected))),
            "b_gradient": np.array([v@residual[ids]/scale/scale for ids, v in derivative])}


def replay(data, y, endpoint, scale):
    result = snapshot_replay(data["table"], len(data["ids"]), y, endpoint)
    ratio = max(1., np.linalg.norm(y))/scale
    result["projected_kkt"] *= ratio; result["b_gradient"] *= ratio*ratio
    return result


def replay_passed(data, y, endpoint, scale):
    if not endpoint or not endpoint.get("valid"): return None
    raw = replay(data, y, endpoint, scale)
    return bool(abs(raw["projected_kkt"]-endpoint["projected_kkt"]) <= 1e-13 and
                np.allclose(raw["b_gradient"], endpoint["b_gradient"], rtol=2e-9, atol=1e-13) and
                abs(np.linalg.norm(raw["residual"])-np.sqrt(endpoint["rss"])) <= 512*np.finfo(float).eps*scale)
