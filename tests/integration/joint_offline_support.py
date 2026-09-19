"""Independent endpoint replay and the retained certification-v2 policy."""
from joint_runtime_support import require
import math
import numpy as np

def predict(basis, beta, n):
    out = np.zeros(n)
    for (ids, values), coefficient in zip(basis, beta): out[ids] += values*coefficient
    return out


def certificate(fit, audit, replay):
    raw = fit.get("qualification_checks", {})
    checks = {k: raw.get(k) for k in ("inner", "b_gradient", "local_correction", "identified")}
    local_source = "legacy-double-QR"
    precision = audit.get("precision")
    if precision and precision.get("agreement_passed") and precision.get("fixed_face_feasible"):
        require(precision["local_correction_passed"] == (precision["local_correction_inf"] <= 1e-10), "Incorrect precision local correction flag.")
        checks["local_correction"] = precision["local_correction_passed"]; local_source = "independent-100-digit-QR"
    checks["derivative"] = None if audit["derivative_status"] == "unavailable" else audit["derivative_verified"]
    checks["endpoint_replay"] = replay
    checks["endpoint_trust"] = audit.get("endpoint_trust", {}).get("passed")
    checks["search_completed"] = fit.get("execution_complete", False) and not fit.get("search_stopped_without_convergence", False)
    if precision:
        from decimal import Decimal
        checks["precision_agreement"] = precision["agreement_passed"]
        checks["precision_feasibility"] = precision.get("fixed_face_feasible")
        checks["precision_kkt"] = Decimal(precision["projected_kkt"]) <= Decimal("1e-10") if "projected_kkt" in precision else None
        checks["precision_b_gradient"] = max(map(abs, precision["b_gradient"])) <= 1e-12 if "b_gradient" in precision else None
    # Strict complementarity is required for a regular boundary certificate.
    constraints = audit.get("constraint_evidence", [])
    weak_boundary = any(c["exactly_active"] and (c["dual"] is None or c["dual"] <= 1e-10) for c in constraints)
    checks["regular_active_face"] = not weak_boundary if constraints else None
    failures = [k for k, v in checks.items() if v is False]
    unavailable = [k for k, v in checks.items() if v is None]
    regular = all(v is True for v in checks.values())
    return {"schema_version": 2, "checks": checks, "failures": failures, "unavailable": unavailable,
            "regular_qualified": regular, "derivative_status": audit["derivative_status"],
            "local_correction_source": local_source,
            "identifiability_scope": "endpoint local numerical rank; generating-model evidence is separate",
            "truth_used": False, "boundary_status": "weakly-active-nonregular" if weak_boundary else "interior-or-strict-face", "evidence": evidence(fit, audit),
            "unavailable_reasons": {k: audit.get("reason", "not executed") for k in unavailable}}


def evidence(fit, audit):
    entries = []
    def add(name, value, threshold, source, comparison="<="):
        entries.append({"check": name, "value": value, "threshold": threshold,
                        "comparison": comparison, "source": source,
                        "status": "unavailable" if value is None else "pass" if
                        (value == threshold if comparison == "==" else value <= threshold) else "fail"})
    for label in ("primary", "reference"):
        endpoint = fit.get(label, {})
        for name, threshold, comparison in (("feasible", True, "=="), ("projected_kkt", 1e-10, "<="), ("b_gradient_inf", 1e-12, "<=")):
            add(label+"_"+name, endpoint.get(name), threshold, "fit."+label+"."+name, comparison)
    add("coefficient_agreement", fit.get("scaled_reference_difference"), 1e-10, "fit.scaled_reference_difference")
    for name, expected in (("design_spectrum", len(fit.get("primary", {}).get("beta", []))),
                           ("width_spectrum", len(fit.get("primary", {}).get("b", []))),
                           ("profile_jacobian_spectrum", len(fit.get("primary", {}).get("b", [])))):
        add(name+"_rank", fit.get(name, {}).get("rank"), expected, "fit."+name+".rank", "==")
    add("legacy_local_correction", fit.get("local_correction_inf"), 1e-10, "fit.local_correction_inf")
    if "precision" in audit:
        add("precision_local_correction", audit["precision"].get("local_correction_inf"), 1e-10, "audit.precision.local_correction_inf")
        for name in ("valid50", "valid100", "agreement_passed", "fixed_face_feasible", "derivative_passed"):
            add("fixed_face_reference_"+name, audit["precision"].get(name), True, "audit.precision."+name, "==")
        p = audit["precision"]
        add("precision_projected_kkt", float(p["projected_kkt"]) if "projected_kkt" in p else None, 1e-10, "audit.precision.projected_kkt")
        add("precision_b_gradient_inf", max(map(abs, p["b_gradient"])) if "b_gradient" in p else None, 1e-12, "audit.precision.b_gradient")
    for ladder in audit.get("ladders", []):
        for key, threshold in (("estimated_relative_error", 1e-7), ("first_relative_error", 1e-6), ("second_relative_error", 1e-6)):
            add("direction_"+str(ladder["direction"])+"_"+key, ladder[key], threshold, "audit.ladders."+str(ladder["direction"])+"."+key)
    return entries


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


def validate_audit(audit, fit):
    require(audit["schema_version"] == 2 and audit["case"] == fit["case"] and audit["dataset"] == fit["dataset"] and
            audit["legacy_joint_qualified"] == fit["joint_qualified"], "Wrong audit identity.")
    ladders = audit.get("ladders", [])
    require(not audit.get("precision") or len(ladders) == 3, "Precision target is missing its full derivative ladder.")
    if ladders:
        require(len(ladders) == 3, "Missing derivative direction.")
        for ladder in ladders:
            require(len(ladder["samples"]) == 17 and [s["h"] for s in ladder["samples"]] == [.01*2**-k for k in range(17)], "Changed step ladder.")
            candidates = [c for c in ladder["candidates"] if c["estimated_relative_error"] is not None]
            selected = min(candidates, key=lambda c: (c["estimated_relative_error"], c["index"])) if candidates else None
            require(ladder["selected"] == (selected["index"] if selected else None), "Step selection is not independent minimum estimated error.")
            passed = selected is not None and selected["estimated_relative_error"] <= 1e-7 and ladder["first_relative_error"] <= 1e-6 and ladder["second_relative_error"] <= 1e-6
            require(ladder["passed"] == passed, "Incorrect Richardson gate.")
        precision = audit.get("precision", {})
        passed = all(d["passed"] for d in ladders) and precision.get("agreement_passed", False) and precision.get("derivative_passed", False)
        require(audit["derivative_verified"] == passed, "False derivative certificate.")
    elif audit["derivative_verified"]:
        require(fit.get("derivative_verified") is True, "Missing derivative evidence.")
    if audit.get("precision", {}).get("valid50") and audit["precision"]["valid100"]:
        from decimal import Decimal
        p = audit["precision"]
        require(p["agreement_passed"] == (Decimal(p["maximum_scaled_precision_difference"]) <= Decimal("1e-20")), "False precision agreement.")
        require(p["derivative_passed"] == (p["fixed_face_feasible"] and max(p["derivative_relative_errors"]) <= 1e-6), "False high precision derivative flag.")


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
