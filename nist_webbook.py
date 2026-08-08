"""
NIST Chemistry WebBook (free, no auth) - reference physicochemical & GC
retention-index data for validation/benchmarking.
https://webbook.nist.gov/chemistry/

NIST does not expose a formal JSON REST API; this client does a best-effort
HTML fetch + light parsing of the public search page. For production use,
prefer downloading the NIST Retention Index Database bulk files and loading
them locally (see `load_local_ri_database`).
"""
import re
import requests

SEARCH_URL = "https://webbook.nist.gov/cgi/cbook.cgi"


def search_by_name(name: str) -> dict | None:
    params = {"Name": name, "Units": "SI"}
    r = requests.get(SEARCH_URL, params=params, timeout=10)
    if r.status_code != 200:
        return None
    html = r.text
    cas_match = re.search(r"CAS Registry Number:</strong>\s*([\d-]+)", html)
    formula_match = re.search(r"Formula:</strong>\s*([A-Za-z0-9]+)", html)
    return {
        "name": name,
        "cas": cas_match.group(1) if cas_match else None,
        "formula": formula_match.group(1) if formula_match else None,
        "source_url": r.url,
    }


def load_local_ri_database(csv_path: str):
    """
    Load a locally downloaded NIST Retention Index Database export (CSV)
    for offline validation of predicted retention against reference GC/LC RI values.
    Expected columns: compound_name, cas, retention_index, column_type, conditions
    """
    import csv
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows
