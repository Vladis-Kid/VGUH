"""
Loaders for free, open GitHub-hosted chromatography datasets/tools used to
validate the simulator against real experimental data (per project brief:
"compare simulation vs several dozen test mixtures").

All fetches use raw.githubusercontent.com — no auth needed for public repos.
Swap REPO_* constants for a fork/mirror if you vendor the data into your
own GitHub project (recommended for reproducibility / rate-limit safety).
"""
import requests
import io
import csv

# Known open repositories relevant to this project (see integration report):
KNOWN_REPOS = {
    "hplc_simulator_excel": "https://github.com/mlibby/hplc_simulator",
    "chrom_rs": "https://github.com/biface/chromatography",
    "retip": "https://github.com/oloBion/Retip",
    "rt_transformer": "https://github.com/Qiong-Yang/RT-Transformer",
    "mzrtsim": "https://github.com/yufree/mzrtsim",
    "openchrom": "https://github.com/OpenChrom/openchrom",
    "chromatography_modeling": "https://github.com/icredd-cheminfo/chromatography-modeling",
    "chromatographic_data": "https://github.com/mpho-mafata/Chromatographic-data",
    "analytical_chem_tools": "https://github.com/pb-cdunn/analytical-chemistry-tools-calculators-and-simulators",
    "gc2asm": "https://github.com/IFPen/GC2ASM",
}


def fetch_raw_csv(repo_raw_url: str) -> list[dict]:
    """
    repo_raw_url must be a raw.githubusercontent.com URL to a CSV file, e.g.:
    https://raw.githubusercontent.com/<org>/<repo>/<branch>/<path>/data.csv
    """
    r = requests.get(repo_raw_url, timeout=15)
    r.raise_for_status()
    reader = csv.DictReader(io.StringIO(r.text))
    return list(reader)


def validate_against_reference(simulated_peaks: list[dict], reference_csv_rows: list[dict],
                                name_key: str = "compound_name", tR_key: str = "tR_observed") -> dict:
    """
    Compares simulator tR predictions against an experimental reference table
    (loaded via fetch_raw_csv or nist_webbook.load_local_ri_database).
    Returns per-compound error and aggregate RMSE/MAE — the core of the
    "simulation vs real data" validation step required by the project brief.
    """
    ref_by_name = {row[name_key]: float(row[tR_key]) for row in reference_csv_rows if row.get(tR_key)}
    errors = []
    per_compound = []
    for p in simulated_peaks:
        name = p["name"]
        if name in ref_by_name:
            err = p["tR"] - ref_by_name[name]
            errors.append(err)
            per_compound.append({
                "name": name, "tR_sim": p["tR"], "tR_ref": ref_by_name[name],
                "abs_error_min": round(abs(err), 4),
                "rel_error_pct": round(abs(err) / ref_by_name[name] * 100, 2) if ref_by_name[name] else None,
            })
    if not errors:
        return {"per_compound": [], "rmse": None, "mae": None}
    rmse = (sum(e ** 2 for e in errors) / len(errors)) ** 0.5
    mae = sum(abs(e) for e in errors) / len(errors)
    return {"per_compound": per_compound, "rmse": round(rmse, 4), "mae": round(mae, 4), "n": len(errors)}
