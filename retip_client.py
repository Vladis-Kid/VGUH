"""
Bridge to Retip / pyRetip (open-source ML retention-time prediction).
Repo: https://github.com/oloBion/Retip  (R + Python 'pyRetip')

This module is OPTIONAL: if pyRetip is not installed, `predict_rt` raises
ImportError with install instructions instead of crashing the simulator.

Typical use in this project: predict a seed retention time for a new
compound from its structure (SMILES), then convert to k_w/S by fitting
two isocratic simulator runs to match the predicted tR (see `calibrate_k_w_S`).
"""
from typing import Optional


def predict_rt(smiles: str, model_name: str = "default") -> float:
    try:
        import pyretip  # type: ignore
    except ImportError as e:
        raise ImportError(
            "pyRetip is not installed. Install with: pip install pyretip\n"
            "Repo: https://github.com/oloBion/Retip"
        ) from e
    model = pyretip.load_model(model_name)
    return float(model.predict([smiles])[0])


def calibrate_k_w_S(predicted_tR: float, t0: float, phi_low: float = 0.1, phi_high: float = 0.9) -> dict:
    """
    Given a single predicted tR (e.g. from Retip/RT-Transformer) at a *reference*
    isocratic composition, back-calculate an approximate k_w assuming a typical
    S (since one data point cannot resolve two unknowns). This is a seed value,
    intended to be refined once >=2 real/predicted retention times are available
    (then k_w and S can be solved exactly from the LSS linear system).
    """
    import math
    k_ref = max((predicted_tR - t0) / max(t0, 1e-9), 1e-3)
    S_assumed = 4.0  # typical small-molecule RP default
    phi_ref = (phi_low + phi_high) / 2.0
    log_k_w = math.log10(k_ref) + S_assumed * phi_ref
    return {"k_w": round(10 ** log_k_w, 3), "S": S_assumed}


def solve_k_w_S_from_two_points(tR1: float, phi1: float, tR2: float, phi2: float, t0: float) -> dict:
    """Exact solution of log10(k) = log10(k_w) - S*phi from two (phi, tR) pairs."""
    import math
    k1 = max((tR1 - t0) / t0, 1e-6)
    k2 = max((tR2 - t0) / t0, 1e-6)
    S = (math.log10(k1) - math.log10(k2)) / (phi2 - phi1)
    log_k_w = math.log10(k1) + S * phi1
    return {"k_w": round(10 ** log_k_w, 4), "S": round(S, 4)}
