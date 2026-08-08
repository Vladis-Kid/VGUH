"""
hplc_sim.metrics
=================
Derived chromatographic figures of merit.
"""
from __future__ import annotations
import math
from typing import List, Dict


def retention_factor(tR: float, t0: float) -> float:
    """k' = (tR - t0) / t0"""
    return (tR - t0) / max(t0, 1e-9)


def selectivity(k1: float, k2: float) -> float:
    """alpha = k2 / k1 (k2 > k1, adjacent peaks)"""
    lo, hi = sorted([k1, k2])
    return hi / max(lo, 1e-9)


def plate_number(tR: float, width_base: float = None, width_half: float = None) -> float:
    """
    N = 16*(tR/Wb)^2   (baseline width, tangent method)
    or N = 5.54*(tR/W0.5)^2  (width at half height, USP/EP standard)
    """
    if width_half is not None:
        return 5.54 * (tR / width_half) ** 2
    if width_base is not None:
        return 16.0 * (tR / width_base) ** 2
    raise ValueError("Provide width_base or width_half")


def hetp(column_length_mm: float, N: float) -> float:
    """Height Equivalent to a Theoretical Plate (micrometers)."""
    L_um = column_length_mm * 1000.0
    return L_um / max(N, 1e-9)


def resolution(tR1: float, tR2: float, w_half1: float, w_half2: float) -> float:
    """
    Rs = 1.18 * (tR2 - tR1) / (w_h1 + w_h2)   [USP, half-height widths]
    """
    return 1.18 * abs(tR2 - tR1) / max((w_half1 + w_half2), 1e-9)


def resolution_baseline(tR1: float, tR2: float, w_b1: float, w_b2: float) -> float:
    """Rs = 2*(tR2 - tR1) / (Wb1 + Wb2)  [classical baseline-width form]"""
    return 2.0 * abs(tR2 - tR1) / max((w_b1 + w_b2), 1e-9)


def asymmetry_factor(peak_left_half_width: float, peak_right_half_width: float) -> float:
    """
    As (USP tailing factor at 10% height uses similar geometry):
    As = b / a, where a,b are left/right widths at 10% (or 5%) peak height
    from the perpendicular dropped from the apex.
    As = 1 -> symmetric; As > 1 -> tailing; As < 1 -> fronting.
    """
    return peak_right_half_width / max(peak_left_half_width, 1e-9)


def confidence_interval(values: List[float], confidence: float = 0.95) -> Dict[str, float]:
    """Simple normal-approximation CI for repeated-run statistics (n replicate injections)."""
    n = len(values)
    if n == 0:
        return {"mean": float("nan"), "std": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    mean = sum(values) / n
    if n > 1:
        var = sum((v - mean) ** 2 for v in values) / (n - 1)
        std = math.sqrt(var)
    else:
        std = 0.0
    # z ~ 1.96 for 95% (normal approx; for small n a t-table would be more rigorous)
    z = 1.96 if confidence == 0.95 else 2.576 if confidence == 0.99 else 1.645
    se = std / math.sqrt(n) if n > 0 else 0.0
    return {
        "mean": mean, "std": std, "n": n,
        "ci_low": mean - z * se, "ci_high": mean + z * se, "rsd_percent": (std / mean * 100.0) if mean else 0.0
    }


def system_suitability(peaks: List[Dict]) -> Dict:
    """
    peaks: list of dicts with keys tR, N, As, w_half, k
    Returns aggregate system-suitability style summary.
    """
    if not peaks:
        return {}
    Ns = [p["N"] for p in peaks if "N" in p]
    Rs_values = []
    sorted_peaks = sorted(peaks, key=lambda p: p["tR"])
    for i in range(len(sorted_peaks) - 1):
        p1, p2 = sorted_peaks[i], sorted_peaks[i + 1]
        w1 = p1.get("w_half_min", p1.get("w_half"))
        w2 = p2.get("w_half_min", p2.get("w_half"))
        Rs_values.append(resolution(p1["tR"], p2["tR"], w1, w2))
    return {
        "N_min": min(Ns) if Ns else None,
        "N_mean": sum(Ns) / len(Ns) if Ns else None,
        "Rs_min": min(Rs_values) if Rs_values else None,
        "Rs_list": Rs_values,
    }
