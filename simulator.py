"""
hplc_sim.simulator
===================
HPLCSimulator: central orchestrator.

Design goal: allow "live edits" (flow rate, mobile-phase composition,
temperature, added components, gradient shape) WITHOUT rebuilding the
whole object graph. Column/MobilePhase/Method/Detector are held as
mutable references; `patch_*` methods mutate in place and the simulator
recomputes only what changed the next time `run()` is called.
"""
from __future__ import annotations
from dataclasses import asdict
from typing import List, Dict, Optional
import numpy as np

from .components import Column, MobilePhase, Compound, Method, Detector, GradientStep
from .models import k_isocratic, van_t_hoff_correction, gradient_retention_time, emg_peak
from .metrics import (
    retention_factor, selectivity, plate_number, hetp, resolution,
    asymmetry_factor, confidence_interval, system_suitability,
)


class HPLCSimulator:
    def __init__(
        self,
        column: Optional[Column] = None,
        mobile_phase: Optional[MobilePhase] = None,
        method: Optional[Method] = None,
        detector: Optional[Detector] = None,
        compounds: Optional[List[Compound]] = None,
    ):
        self.column = column or Column()
        self.mobile_phase = mobile_phase or MobilePhase()
        self.method = method or Method()
        self.detector = detector or Detector()
        self.compounds: List[Compound] = compounds or []
        self._t_ref_K = 298.15

    # ------------------------------------------------------------------
    # LIVE-EDIT API (no restart required)
    # ------------------------------------------------------------------
    def patch_flow(self, flow_ml_min: float):
        self.mobile_phase.flow_ml_min = flow_ml_min

    def patch_temperature(self, temperature_C: float):
        self.mobile_phase.temperature_C = temperature_C

    def patch_mobile_phase(self, **kwargs):
        for k, v in kwargs.items():
            if hasattr(self.mobile_phase, k):
                setattr(self.mobile_phase, k, v)

    def patch_method(self, **kwargs):
        for k, v in kwargs.items():
            if hasattr(self.method, k):
                setattr(self.method, k, v)

    def set_gradient(self, profile: List[Dict]):
        """profile: [{'time_min':0,'phi_B':0.05}, {'time_min':10,'phi_B':0.95}, ...]"""
        self.method.mode = "gradient"
        self.method.gradient_profile = [GradientStep(**p) for p in profile]

    def set_step_gradient(self, steps: List[Dict]):
        """Step elution: list of {'time_min':t,'phi_B':phi} held constant between jumps."""
        self.method.mode = "step"
        expanded = []
        for i, s in enumerate(steps):
            expanded.append(GradientStep(s["time_min"], s["phi_B"]))
            if i + 1 < len(steps):
                # hold value until just before next step (creates a sharp step, not a ramp)
                next_t = steps[i + 1]["time_min"]
                expanded.append(GradientStep(next_t - 1e-3, s["phi_B"]))
        self.method.gradient_profile = expanded

    def add_component(self, compound: Compound):
        self.compounds.append(compound)

    def remove_component(self, name: str):
        self.compounds = [c for c in self.compounds if c.name != name]

    def update_component(self, name: str, **kwargs):
        for c in self.compounds:
            if c.name == name:
                for k, v in kwargs.items():
                    if hasattr(c, k):
                        setattr(c, k, v)

    # ------------------------------------------------------------------
    # CORE COMPUTATION
    # ------------------------------------------------------------------
    def _dead_time(self) -> float:
        return self.column.dead_time_min(self.mobile_phase.flow_ml_min)

    def _apply_temperature(self, compound: Compound, k_ref: float) -> float:
        T = self.mobile_phase.temperature_K
        return van_t_hoff_correction(k_ref, self._t_ref_K, T, compound.dH_J_mol)

    def compute_retention(self, compound: Compound) -> float:
        t0 = self._dead_time()
        # temperature-adjust k_w at phi=0 reference, then propagate through model
        k_w_T = self._apply_temperature(compound, compound.k_w)

        if self.method.mode == "isocratic":
            k = k_isocratic(k_w_T, compound.S, self.method.isocratic_phi)
            tR = t0 * (1 + k)
        else:
            profile = self.method.profile_tuples()
            tR = gradient_retention_time(
                k_w=k_w_T,
                S=compound.S,
                t0=t0,
                t_dwell=self.method.dwell_time_min,
                flow_ml_min=self.mobile_phase.flow_ml_min,
                gradient_profile=profile,
            )
        return tR

    def efficiency(self, tR: float) -> float:
        """Plate number scaled from nominal column efficiency, degraded slightly by
        retention factor (band broadening) and flow-rate deviation from optimum."""
        t0 = self._dead_time()
        k = retention_factor(tR, t0)
        N_nominal = self.column.N_per_m_nominal * (self.column.length_mm / 1000.0)
        # mild extra-column / kinetic broadening correction (empirical, bounded)
        broadening_factor = 1.0 / (1.0 + 0.02 * max(k, 0))
        return max(N_nominal * broadening_factor, 100.0)

    def run(self, n_points: int = 6000, add_noise: bool = True, seed: Optional[int] = None) -> Dict:
        """
        Execute the simulation for the current state (column/mobile phase/
        method/detector/compounds) and return a full result payload:
        peak table + time/signal arrays for the chromatogram.
        """
        if seed is not None:
            np.random.seed(seed)

        t0 = self._dead_time()
        t_end = self.method.run_time_min
        t = np.linspace(0, t_end, n_points)
        signal = np.zeros_like(t)

        peaks = []
        for compound in self.compounds:
            tR = self.compute_retention(compound)
            N = self.efficiency(tR)
            sigma = tR / np.sqrt(max(N, 1.0))
            tau = compound.tau_rel * sigma
            area = 1000.0 * compound.response_factor
            y = emg_peak(t, area=area, tR=tR, sigma=sigma, tau=tau)
            signal += y

            k = retention_factor(tR, t0)
            w_half = 2.355 * sigma  # FWHM for a near-Gaussian approximation
            w_base = 4.0 * sigma
            height = float(np.max(y))
            peaks.append({
                "name": compound.name,
                "tR": round(float(tR), 3),
                "k": round(float(k), 3),
                "N": round(float(N), 0),
                "HETP_um": round(hetp(self.column.length_mm, N), 2),
                "area": round(float(area), 2),
                "height_mAU": round(height, 3),
                "w_half_min": round(float(w_half), 4),
                "w_base_min": round(float(w_base), 4),
                "sigma": float(sigma),
                "tau": float(tau),
            })

        # adjacent-peak selectivity & resolution
        peaks_sorted = sorted(peaks, key=lambda p: p["tR"])
        for i in range(len(peaks_sorted) - 1):
            p1, p2 = peaks_sorted[i], peaks_sorted[i + 1]
            k1 = p1["k"] if p1["k"] > 0 else 1e-6
            k2 = p2["k"] if p2["k"] > 0 else 1e-6
            p2["alpha_vs_prev"] = round(selectivity(k1, k2), 4)
            p2["Rs_vs_prev"] = round(
                resolution(p1["tR"], p2["tR"], p1["w_half_min"], p2["w_half_min"]), 3
            )

        if add_noise:
            noise = np.random.normal(0, self.detector.noise_std, size=t.shape)
            drift = self.detector.baseline_drift_per_min * t
            signal = signal + noise + drift

        suitability = system_suitability(peaks_sorted)

        return {
            "time_min": t.tolist(),
            "signal_mAU": signal.tolist(),
            "t0_min": round(t0, 4),
            "peaks": peaks_sorted,
            "system_suitability": suitability,
            "column": asdict(self.column),
            "mobile_phase": asdict(self.mobile_phase),
            "method": {
                "mode": self.method.mode,
                "profile": self.method.profile_tuples(),
                "run_time_min": self.method.run_time_min,
            },
            "detector": asdict(self.detector),
        }

    def run_replicates(self, n: int = 6, **kwargs) -> Dict:
        """Repeat injections with different noise seeds -> statistics per compound (RSD, CI)."""
        tR_by_name: Dict[str, List[float]] = {}
        area_by_name: Dict[str, List[float]] = {}
        for i in range(n):
            result = self.run(seed=i, **kwargs)
            for p in result["peaks"]:
                tR_by_name.setdefault(p["name"], []).append(p["tR"])
                area_by_name.setdefault(p["name"], []).append(p["area"])
        stats = {}
        for name in tR_by_name:
            stats[name] = {
                "tR": confidence_interval(tR_by_name[name]),
                "area": confidence_interval(area_by_name[name]),
            }
        return stats
