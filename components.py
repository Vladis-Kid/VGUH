"""
hplc_sim.components
====================
Domain model: Column, MobilePhase, Compound, Method (isocratic/gradient/step),
Detector. All are mutable dataclasses so the Simulator can patch fields
in real time (no full restart needed).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Literal
import math


@dataclass
class Column:
    name: str = "Generic C18"
    ctype: Literal["RP", "NP", "IEX", "HILIC", "SEC"] = "RP"
    length_mm: float = 150.0
    id_mm: float = 4.6
    particle_um: float = 3.5
    porosity: float = 0.65          # total porosity (interparticle+intraparticle)
    N_per_m_nominal: float = 80000  # manufacturer efficiency at ref. conditions

    @property
    def volume_ml(self) -> float:
        r_cm = (self.id_mm / 10.0) / 2.0
        l_cm = self.length_mm / 10.0
        return math.pi * r_cm ** 2 * l_cm

    def dead_volume_ml(self) -> float:
        return self.volume_ml * self.porosity

    def dead_time_min(self, flow_ml_min: float) -> float:
        return self.dead_volume_ml() / max(flow_ml_min, 1e-6)


@dataclass
class MobilePhase:
    solvent_A: str = "Water + 0.1% formic acid"
    solvent_B: str = "Acetonitrile"
    pH: float = 2.6
    buffer_conc_mM: float = 0.0
    flow_ml_min: float = 1.0
    temperature_C: float = 25.0

    @property
    def temperature_K(self) -> float:
        return self.temperature_C + 273.15


@dataclass
class Compound:
    name: str
    k_w: float = 5.0        # retention factor in 100% weak solvent
    S: float = 4.0          # LSS solvent-strength parameter
    dH_J_mol: float = -20000.0  # van't Hoff enthalpy for temperature dependence
    tau_rel: float = 0.15   # relative EMG tailing factor (tau/sigma), >0 tailing
    response_factor: float = 1.0  # detector-specific relative response
    logP: float = None      # optional, filled by PubChem integration


@dataclass
class GradientStep:
    time_min: float
    phi_B: float  # fraction of strong solvent (0..1)


@dataclass
class Method:
    mode: Literal["isocratic", "gradient", "step"] = "isocratic"
    isocratic_phi: float = 0.4
    gradient_profile: List[GradientStep] = field(default_factory=list)
    dwell_time_min: float = 0.5
    run_time_min: float = 20.0

    def profile_tuples(self) -> List[Tuple[float, float]]:
        if self.mode == "isocratic":
            return [(0.0, self.isocratic_phi), (self.run_time_min, self.isocratic_phi)]
        return [(g.time_min, g.phi_B) for g in self.gradient_profile]


@dataclass
class Detector:
    dtype: Literal["UV", "MS", "RI", "DAD"] = "UV"
    wavelength_nm: float = 254.0
    noise_std: float = 0.15       # mAU baseline noise
    baseline_drift_per_min: float = 0.01
    threshold_mAU: float = 1.0    # integration cutoff
