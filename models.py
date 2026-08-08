"""
hplc_sim.models
================
Analytical math models for HPLC retention & peak shape.

References (open literature, no proprietary formulas):
- Snyder, L.R.; Dolan, J.W. "High-Performance Gradient Elution" (LSS theory)
- Neue, U.D. "HPLC Columns: Theory, Technology, and Practice"
- Foley, J.P.; Dorsey, J.G. Anal. Chem. 1983, 55, 730 (EMG peak model, asymmetry)
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import List, Callable, Optional

R_GAS = 8.314462618  # J/(mol*K)


# ----------------------------------------------------------------------
# 1. ISOCRATIC RETENTION MODEL (Linear Solvent Strength, LSS)
# ----------------------------------------------------------------------
def k_isocratic(k_w: float, S: float, phi: float) -> float:
    """
    Reverse-phase LSS retention model:
        log10(k) = log10(k_w) - S * phi
    k_w  : retention factor in pure weak solvent (phi=0), compound-specific
    S    : solvent-strength sensitivity, compound-specific (~2-6 for small molecules)
    phi  : volume fraction of strong (organic) solvent [0..1]
    """
    log_k = math.log10(max(k_w, 1e-6)) - S * phi
    return 10 ** log_k


def van_t_hoff_correction(k_ref: float, T_ref: float, T: float, dH: float = -20000.0) -> float:
    """
    Temperature dependence of k' via van't Hoff equation:
        ln(k_T) = ln(k_ref) - (dH/R) * (1/T - 1/T_ref)
    dH   : apparent transfer enthalpy, J/mol (typ. -10 to -30 kJ/mol, negative = adsorption exothermic)
    T, T_ref in Kelvin
    """
    ln_k = math.log(max(k_ref, 1e-6)) - (dH / R_GAS) * (1.0 / T - 1.0 / T_ref)
    return math.exp(ln_k)


# ----------------------------------------------------------------------
# 2. GRADIENT ELUTION MODEL (Snyder-Dolan Linear Solvent Strength gradient)
# ----------------------------------------------------------------------
def gradient_retention_time(
    k_w: float,
    S: float,
    t0: float,
    t_dwell: float,
    flow_ml_min: float,
    gradient_profile: List[tuple],
) -> float:
    """
    Numerically integrate the fundamental gradient-elution equation:

        integral_{0}^{tR - t0 - tD} [1 / (t0 * k(phi(t)))] dt = 1

    gradient_profile: list of (time_min, phi) breakpoints defining phi(t)
                       (phi = fraction of strong solvent B), piecewise-linear,
                       supports isocratic holds and multi-step/step gradients.
    t0      : column dead time (min) = V0 / F
    t_dwell : instrument dwell time (min), delays gradient reaching the column
    Returns retention time tR (min).
    """
    def phi_at(t: float) -> float:
        if t <= gradient_profile[0][0]:
            return gradient_profile[0][1]
        for i in range(len(gradient_profile) - 1):
            t1, p1 = gradient_profile[i]
            t2, p2 = gradient_profile[i + 1]
            if t1 <= t <= t2:
                if t2 == t1:
                    return p2
                frac = (t - t1) / (t2 - t1)
                return p1 + frac * (p2 - p1)
        return gradient_profile[-1][1]

    dt = 0.001  # min, integration step
    integral = 0.0
    t = 0.0
    max_t = 200.0  # safety cap (min)
    while t < max_t:
        t_eff = max(0.0, t - t_dwell)
        phi = phi_at(t_eff)
        k = k_isocratic(k_w, S, phi)
        integral += dt / (t0 * max(k, 1e-6))
        if integral >= 1.0:
            return t + t0  # elution = migration time + one more column volume
        t += dt
    return max_t  # compound never elutes within window -> report cap


# ----------------------------------------------------------------------
# 3. PEAK-SHAPE MODEL: Exponentially Modified Gaussian (EMG)
#    Captures realistic HPLC tailing/fronting (Foley-Dorsey model)
# ----------------------------------------------------------------------
def emg_peak(t: "np.ndarray", area: float, tR: float, sigma: float, tau: float):
    """
    EMG peak profile (convolution of Gaussian + exponential decay -> tailing).
    tau > 0  -> tailing (typical HPLC, adsorption/extra-column effects)
    tau -> 0 -> approaches pure Gaussian
    A small negative-tau approximation is used for fronting by mirroring t.

    Numerically-stable implementation: naive exp(exponent)*erfc(z) overflows
    for large exponent while erfc(z) underflows to ~0 (inf*0 = nan). We use
    the scaled complementary error function erfcx(z) = exp(z^2)*erfc(z) and
    fold z^2 into the exponent before exponentiating, so the product stays
    bounded. For z<0 we use the erfc reflection identity to keep erfcx's
    argument non-negative (erfcx grows very fast for negative arguments).
    """
    import numpy as np
    from scipy import special

    if abs(tau) < 1e-6:
        tau = 1e-6 if tau >= 0 else -1e-6

    sign = 1.0 if tau > 0 else -1.0
    tt = t if sign > 0 else (2 * tR - t)
    tau_abs = abs(tau)

    z = (tR + sign * (sigma ** 2) / tau_abs - tt) / (sigma * math.sqrt(2))
    prefactor = area / (2 * tau_abs)
    exponent = (sigma ** 2) / (2 * tau_abs ** 2) + sign * (tR - tt) / tau_abs

    y = np.zeros_like(np.asarray(t, dtype=float))
    z = np.asarray(z, dtype=float)
    exponent = np.asarray(exponent, dtype=float)

    with np.errstate(over="ignore", invalid="ignore"):
        pos = z >= 0
        # z >= 0: erfc(z) = exp(-z^2) * erfcx(z)  -> combine exponents first
        exp_arg_pos = np.clip(exponent[pos] - z[pos] ** 2, -700, 700)
        y[pos] = prefactor * np.exp(exp_arg_pos) * special.erfcx(z[pos])

        # z < 0: erfc(z) = 2 - erfc(-z) = 2 - exp(-z^2)*erfcx(-z)
        zn = -z[~pos]
        exp_arg_full = np.clip(exponent[~pos], -700, 700)
        exp_arg_neg = np.clip(exponent[~pos] - zn ** 2, -700, 700)
        y[~pos] = prefactor * (2.0 * np.exp(exp_arg_full) - np.exp(exp_arg_neg) * special.erfcx(zn))

    return np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)


def gaussian_peak(t, area: float, tR: float, sigma: float):
    import numpy as np
    return (area / (sigma * math.sqrt(2 * math.pi))) * np.exp(-((t - tR) ** 2) / (2 * sigma ** 2))
