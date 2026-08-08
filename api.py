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
