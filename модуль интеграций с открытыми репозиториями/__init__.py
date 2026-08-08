from .components import Column, MobilePhase, Compound, Method, Detector, GradientStep
from .simulator import HPLCSimulator
from . import metrics, models, optimization

__all__ = [
    "Column", "MobilePhase", "Compound", "Method", "Detector", "GradientStep",
    "HPLCSimulator", "metrics", "models", "optimization",
]
__version__ = "0.1.0"


"""
hplc_sim.api
============
FastAPI application exposing the simulator over REST, so the web dashboard
(dashboard/index.html) or any external LIMS/client can:
  - create a simulation session,
  - PATCH parameters live (flow, mobile phase, temperature, gradient, components),
  - trigger a run and get back the full chromatogram + peak table + metrics.

Run:  uvicorn hplc_sim.api:app --reload --port 8000
"""
from __future__ import annotations
from typing import Dict, List, Optional
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from .components import Column, MobilePhase, Method, Detector, Compound
from .simulator import HPLCSimulator
from .optimization import monte_carlo_optimize, genetic_optimize

app = FastAPI(title="HPLC Simulation & Modeling Module", version="0.1.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

# in-memory session store (swap for Redis/DB in production)
SESSIONS: Dict[str, HPLCSimulator] = {}


class ColumnIn(BaseModel):
    name: str = "Generic C18"
    ctype: str = "RP"
    length_mm: float = 150.0
    id_mm: float = 4.6
    particle_um: float = 3.5
    porosity: float = 0.65
    N_per_m_nominal: float = 80000


class MobilePhaseIn(BaseModel):
    solvent_A: str = "Water + 0.1% formic acid"
    solvent_B: str = "Acetonitrile"
    pH: float = 2.6
    buffer_conc_mM: float = 0.0
    flow_ml_min: float = 1.0
    temperature_C: float = 25.0


class CompoundIn(BaseModel):
    name: str
    k_w: float = 5.0
    S: float = 4.0
    dH_J_mol: float = -20000.0
    tau_rel: float = 0.15
    response_factor: float = 1.0


class CreateSessionIn(BaseModel):
    session_id: str
    column: ColumnIn = ColumnIn()
    mobile_phase: MobilePhaseIn = MobilePhaseIn()
    compounds: List[CompoundIn] = []
    mode: str = "isocratic"
    isocratic_phi: float = 0.4
    gradient_profile: Optional[List[Dict]] = None
    run_time_min: float = 20.0


class PatchIn(BaseModel):
    flow_ml_min: Optional[float] = None
    temperature_C: Optional[float] = None
    mobile_phase: Optional[Dict] = None
    method: Optional[Dict] = None
    gradient_profile: Optional[List[Dict]] = None
    add_compound: Optional[CompoundIn] = None
    remove_compound: Optional[str] = None
    update_compound: Optional[Dict] = None  # {"name": "...", **fields}


def _get_session(session_id: str) -> HPLCSimulator:
    sim = SESSIONS.get(session_id)
    if sim is None:
        raise HTTPException(404, f"Session '{session_id}' not found. Create it first.")
    return sim


@app.post("/sessions")
def create_session(payload: CreateSessionIn):
    column = Column(**payload.column.model_dump())
    mobile_phase = MobilePhase(**payload.mobile_phase.model_dump())
    method = Method(mode=payload.mode, isocratic_phi=payload.isocratic_phi,
                     run_time_min=payload.run_time_min)
    detector = Detector()
    compounds = [Compound(**c.model_dump()) for c in payload.compounds]
    sim = HPLCSimulator(column, mobile_phase, method, detector, compounds)
    if payload.gradient_profile:
        sim.set_gradient(payload.gradient_profile)
    SESSIONS[payload.session_id] = sim
    return {"status": "created", "session_id": payload.session_id}


@app.patch("/sessions/{session_id}")
def patch_session(session_id: str, patch: PatchIn):
    """Live edit without restarting the simulation session."""
    sim = _get_session(session_id)
    if patch.flow_ml_min is not None:
        sim.patch_flow(patch.flow_ml_min)
    if patch.temperature_C is not None:
        sim.patch_temperature(patch.temperature_C)
    if patch.mobile_phase:
        sim.patch_mobile_phase(**patch.mobile_phase)
    if patch.method:
        sim.patch_method(**patch.method)
    if patch.gradient_profile:
        sim.set_gradient(patch.gradient_profile)
    if patch.add_compound:
        sim.add_component(Compound(**patch.add_compound.model_dump()))
    if patch.remove_compound:
        sim.remove_component(patch.remove_compound)
    if patch.update_compound:
        name = patch.update_compound.pop("name")
        sim.update_component(name, **patch.update_compound)
    return {"status": "patched"}


@app.post("/sessions/{session_id}/run")
def run_session(session_id: str, n_points: int = 4000, add_noise: bool = True, seed: Optional[int] = None):
    sim = _get_session(session_id)
    return sim.run(n_points=n_points, add_noise=add_noise, seed=seed)


@app.post("/sessions/{session_id}/replicates")
def run_replicates(session_id: str, n: int = 6):
    sim = _get_session(session_id)
    return sim.run_replicates(n=n)


@app.post("/sessions/{session_id}/optimize")
def optimize(session_id: str, algorithm: str = "genetic", iterations: int = 15):
    sim = _get_session(session_id)
    if algorithm == "monte_carlo":
        return monte_carlo_optimize(sim, n_iter=max(50, iterations * 10))
    return genetic_optimize(sim, generations=iterations)


@app.delete("/sessions/{session_id}")
def delete_session(session_id: str):
    SESSIONS.pop(session_id, None)
    return {"status": "deleted"}


@app.get("/health")
def health():
    return {"status": "ok", "active_sessions": len(SESSIONS)}



import math
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hplc_sim import Column, MobilePhase, Method, Compound, Detector, HPLCSimulator
from hplc_sim.metrics import retention_factor, selectivity, plate_number, resolution


def test_isocratic_monotonic_with_phi():
    """Higher %organic (phi) must reduce retention (weaker retention on RP)."""
    from hplc_sim.models import k_isocratic
    k_low_phi = k_isocratic(k_w=10, S=4, phi=0.2)
    k_high_phi = k_isocratic(k_w=10, S=4, phi=0.6)
    assert k_high_phi < k_low_phi


def test_retention_factor_formula():
    assert math.isclose(retention_factor(tR=5.0, t0=1.0), 4.0)


def test_selectivity_ge_one():
    alpha = selectivity(k1=2.0, k2=3.0)
    assert alpha >= 1.0
    assert math.isclose(alpha, 1.5)


def test_plate_number_positive():
    N = plate_number(tR=5.0, width_half=0.2)
    assert N > 0


def test_resolution_zero_for_identical_peaks():
    rs = resolution(tR1=5.0, tR2=5.0, w_half1=0.2, w_half2=0.2)
    assert rs == 0


def test_simulator_end_to_end_isocratic():
    col = Column(length_mm=150, id_mm=4.6, particle_um=3.5)
    mp = MobilePhase(flow_ml_min=1.0, temperature_C=25)
    method = Method(mode="isocratic", isocratic_phi=0.4, run_time_min=15)
    compounds = [
        Compound(name="A", k_w=8, S=4),
        Compound(name="B", k_w=15, S=4.2),
    ]
    sim = HPLCSimulator(col, mp, method, Detector(), compounds)
    result = sim.run(n_points=1000, add_noise=False)
    assert len(result["peaks"]) == 2
    # B should elute after A (higher k_w)
    tR_A = [p["tR"] for p in result["peaks"] if p["name"] == "A"][0]
    tR_B = [p["tR"] for p in result["peaks"] if p["name"] == "B"][0]
    assert tR_B > tR_A


def test_live_edit_flow_changes_retention_time():
    col = Column()
    mp = MobilePhase(flow_ml_min=1.0)
    method = Method(mode="isocratic", isocratic_phi=0.4, run_time_min=20)
    sim = HPLCSimulator(col, mp, method, Detector(), [Compound(name="X", k_w=10, S=4)])
    tR_before = sim.compute_retention(sim.compounds[0])
    sim.patch_flow(2.0)  # double flow -> ~half retention time
    tR_after = sim.compute_retention(sim.compounds[0])
    assert tR_after < tR_before


def test_gradient_mode_runs():
    col = Column()
    mp = MobilePhase(flow_ml_min=1.0)
    method = Method(mode="gradient", run_time_min=20, dwell_time_min=0.3)
    sim = HPLCSimulator(col, mp, method, Detector(), [Compound(name="G1", k_w=20, S=5)])
    sim.set_gradient([{"time_min": 0, "phi_B": 0.05}, {"time_min": 15, "phi_B": 0.95}])
    result = sim.run(n_points=800, add_noise=False)
    assert result["peaks"][0]["tR"] > 0



import math
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hplc_sim import Column, MobilePhase, Method, Compound, Detector, HPLCSimulator
from hplc_sim.metrics import retention_factor, selectivity, plate_number, resolution


def test_isocratic_monotonic_with_phi():
    """Higher %organic (phi) must reduce retention (weaker retention on RP)."""
    from hplc_sim.models import k_isocratic
    k_low_phi = k_isocratic(k_w=10, S=4, phi=0.2)
    k_high_phi = k_isocratic(k_w=10, S=4, phi=0.6)
    assert k_high_phi < k_low_phi


def test_retention_factor_formula():
    assert math.isclose(retention_factor(tR=5.0, t0=1.0), 4.0)


def test_selectivity_ge_one():
    alpha = selectivity(k1=2.0, k2=3.0)
    assert alpha >= 1.0
    assert math.isclose(alpha, 1.5)


def test_plate_number_positive():
    N = plate_number(tR=5.0, width_half=0.2)
    assert N > 0


def test_resolution_zero_for_identical_peaks():
    rs = resolution(tR1=5.0, tR2=5.0, w_half1=0.2, w_half2=0.2)
    assert rs == 0


def test_simulator_end_to_end_isocratic():
    col = Column(length_mm=150, id_mm=4.6, particle_um=3.5)
    mp = MobilePhase(flow_ml_min=1.0, temperature_C=25)
    method = Method(mode="isocratic", isocratic_phi=0.4, run_time_min=15)
    compounds = [
        Compound(name="A", k_w=8, S=4),
        Compound(name="B", k_w=15, S=4.2),
    ]
    sim = HPLCSimulator(col, mp, method, Detector(), compounds)
    result = sim.run(n_points=1000, add_noise=False)
    assert len(result["peaks"]) == 2
    # B should elute after A (higher k_w)
    tR_A = [p["tR"] for p in result["peaks"] if p["name"] == "A"][0]
    tR_B = [p["tR"] for p in result["peaks"] if p["name"] == "B"][0]
    assert tR_B > tR_A


def test_live_edit_flow_changes_retention_time():
    col = Column()
    mp = MobilePhase(flow_ml_min=1.0)
    method = Method(mode="isocratic", isocratic_phi=0.4, run_time_min=20)
    sim = HPLCSimulator(col, mp, method, Detector(), [Compound(name="X", k_w=10, S=4)])
    tR_before = sim.compute_retention(sim.compounds[0])
    sim.patch_flow(2.0)  # double flow -> ~half retention time
    tR_after = sim.compute_retention(sim.compounds[0])
    assert tR_after < tR_before


def test_gradient_mode_runs():
    col = Column()
    mp = MobilePhase(flow_ml_min=1.0)
    method = Method(mode="gradient", run_time_min=20, dwell_time_min=0.3)
    sim = HPLCSimulator(col, mp, method, Detector(), [Compound(name="G1", k_w=20, S=5)])
    sim.set_gradient([{"time_min": 0, "phi_B": 0.05}, {"time_min": 15, "phi_B": 0.95}])
    result = sim.run(n_points=800, add_noise=False)
    assert result["peaks"][0]["tR"] > 0




import math
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hplc_sim import Column, MobilePhase, Method, Compound, Detector, HPLCSimulator
from hplc_sim.metrics import retention_factor, selectivity, plate_number, resolution


def test_isocratic_monotonic_with_phi():
    """Higher %organic (phi) must reduce retention (weaker retention on RP)."""
    from hplc_sim.models import k_isocratic
    k_low_phi = k_isocratic(k_w=10, S=4, phi=0.2)
    k_high_phi = k_isocratic(k_w=10, S=4, phi=0.6)
    assert k_high_phi < k_low_phi


def test_retention_factor_formula():
    assert math.isclose(retention_factor(tR=5.0, t0=1.0), 4.0)


def test_selectivity_ge_one():
    alpha = selectivity(k1=2.0, k2=3.0)
    assert alpha >= 1.0
    assert math.isclose(alpha, 1.5)


def test_plate_number_positive():
    N = plate_number(tR=5.0, width_half=0.2)
    assert N > 0


def test_resolution_zero_for_identical_peaks():
    rs = resolution(tR1=5.0, tR2=5.0, w_half1=0.2, w_half2=0.2)
    assert rs == 0


def test_simulator_end_to_end_isocratic():
    col = Column(length_mm=150, id_mm=4.6, particle_um=3.5)
    mp = MobilePhase(flow_ml_min=1.0, temperature_C=25)
    method = Method(mode="isocratic", isocratic_phi=0.4, run_time_min=15)
    compounds = [
        Compound(name="A", k_w=8, S=4),
        Compound(name="B", k_w=15, S=4.2),
    ]
    sim = HPLCSimulator(col, mp, method, Detector(), compounds)
    result = sim.run(n_points=1000, add_noise=False)
    assert len(result["peaks"]) == 2
    # B should elute after A (higher k_w)
    tR_A = [p["tR"] for p in result["peaks"] if p["name"] == "A"][0]
    tR_B = [p["tR"] for p in result["peaks"] if p["name"] == "B"][0]
    assert tR_B > tR_A


def test_live_edit_flow_changes_retention_time():
    col = Column()
    mp = MobilePhase(flow_ml_min=1.0)
    method = Method(mode="isocratic", isocratic_phi=0.4, run_time_min=20)
    sim = HPLCSimulator(col, mp, method, Detector(), [Compound(name="X", k_w=10, S=4)])
    tR_before = sim.compute_retention(sim.compounds[0])
    sim.patch_flow(2.0)  # double flow -> ~half retention time
    tR_after = sim.compute_retention(sim.compounds[0])
    assert tR_after < tR_before


def test_gradient_mode_runs():
    col = Column()
    mp = MobilePhase(flow_ml_min=1.0)
    method = Method(mode="gradient", run_time_min=20, dwell_time_min=0.3)
    sim = HPLCSimulator(col, mp, method, Detector(), [Compound(name="G1", k_w=20, S=5)])
    sim.set_gradient([{"time_min": 0, "phi_B": 0.05}, {"time_min": 15, "phi_B": 0.95}])
    result = sim.run(n_points=800, add_noise=False)
    assert result["peaks"][0]["tR"] > 0



bash

python3 -c "import numpy, scipy; print('numpy/scipy ok')"; pip list 2>/dev/null | grep -iE "numpy|scipy|pytest|fastapi"
Output

numpy/scipy ok
numpy                      2.4.4
scipy                      1.17.1




