"""
Bridge to chrom-rs (Rust chromatography PDE simulator, Langmuir isotherms,
Euler/RK4 solvers). Repo: https://github.com/biface/chromatography

This project's Python model (models.py) uses fast closed-form / ODE-lite
LSS equations for interactive, real-time editing. When higher physical
fidelity is needed (nonlinear/competitive adsorption, overloaded columns,
preparative-scale effects), delegate to chrom-rs as an external process:

    1. Write a YAML config (see `build_chromrs_config`)
    2. Run the compiled chrom-rs binary (`cargo install --path .` from the repo,
       or a prebuilt binary you ship alongside this module)
    3. Parse its CSV/JSON output back into this module's peak schema.

This keeps chrom-rs as an OPTIONAL high-accuracy backend rather than a
hard dependency, matching the "choose model by required precision" design.
"""
import json
import subprocess
import tempfile
import os
from typing import Dict, List


def build_chromrs_config(column: Dict, compounds: List[Dict], flow_ml_min: float,
                          injection_profile: str = "gaussian") -> dict:
    return {
        "column": {
            "length_m": column["length_mm"] / 1000.0,
            "diameter_m": column["id_mm"] / 1000.0,
            "porosity": column.get("porosity", 0.65),
        },
        "flow_m3_s": flow_ml_min * 1e-6 / 60.0,
        "injection": {"profile": injection_profile},
        "solutes": [
            {
                "name": c["name"],
                "isotherm": "langmuir",
                "params": {"k_eq": c.get("k_w", 1.0), "q_max": c.get("q_max", 10.0)},
            }
            for c in compounds
        ],
        "solver": "rk4",
    }


def run_chromrs(config: dict, binary_path: str = "chrom-rs") -> dict:
    """
    Executes the chrom-rs binary against a temp YAML config and returns
    parsed JSON results. Requires the binary to be built/installed separately
    (see repo README: https://github.com/biface/chromatography).
    """
    import yaml  # PyYAML
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.safe_dump(config, f)
        cfg_path = f.name
    out_path = cfg_path.replace(".yaml", ".json")
    try:
        subprocess.run(
            [binary_path, "--config", cfg_path, "--output", out_path, "--format", "json"],
            check=True, capture_output=True, timeout=120,
        )
        with open(out_path) as f:
            return json.load(f)
    finally:
        for p in (cfg_path, out_path):
            if os.path.exists(p):
                os.remove(p)
