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
