"""
Demonstration of the validation workflow required by the project brief:
"compare simulation results against real experimental data (several dozen
test mixtures)".

IMPORTANT: tests/reference_mixture_data.csv ships with illustrative
placeholder values so this script runs end-to-end offline. For a real
validation campaign, replace it with:
  - your own lab's experimental results, or
  - values pulled from the NIST Retention Index Database / ACD Labs
    public DB (see hplc_sim/integrations/nist_webbook.py), or
  - a GitHub-hosted reference set (see hplc_sim/integrations/github_datasets.py)
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hplc_sim import Column, MobilePhase, Method, Compound, Detector, HPLCSimulator
from hplc_sim.integrations.github_datasets import validate_against_reference

CSV_PATH = os.path.join(os.path.dirname(__file__), "reference_mixture_data.csv")

with open(CSV_PATH, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

compounds = [
    Compound(name=r["compound_name"], k_w=float(r["k_w"]), S=float(r["S"]))
    for r in rows if r["compound_name"] != "Uracil"
]

col = Column(length_mm=150, id_mm=4.6, particle_um=3.5)
mp = MobilePhase(flow_ml_min=1.0, temperature_C=25)
method = Method(mode="isocratic", isocratic_phi=0.30, run_time_min=15)
sim = HPLCSimulator(col, mp, method, Detector(), compounds)
result = sim.run(n_points=3000, add_noise=False)

report = validate_against_reference(result["peaks"], rows,
                                     name_key="compound_name", tR_key="tR_observed")

print(f"{'Compound':<18}{'tR_sim':>10}{'tR_ref':>10}{'AbsErr':>10}{'RelErr%':>10}")
for row in report["per_compound"]:
    print(f"{row['name']:<18}{row['tR_sim']:>10.3f}{row['tR_ref']:>10.3f}"
          f"{row['abs_error_min']:>10.3f}{row['rel_error_pct']:>10.2f}")
print(f"\nRMSE = {report['rmse']} min | MAE = {report['mae']} min | n = {report['n']}")
