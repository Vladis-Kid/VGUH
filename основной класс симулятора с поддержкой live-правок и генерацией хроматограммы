"""
Generic LIMS integration stub.
Most LIMS (LabWare, STARLIMS, LabVantage, open-source Bika/SENAITE) accept
results either via REST/JSON, or via standard file interchange formats
(AnIML, JCAMP-DX, ASTM E1394/SiLA2). This module implements the two
lowest-friction, dependency-free paths:

  1. Generic REST push (works with SENAITE, LabWare REST connectors, or a
     thin custom adapter you register with your LIMS vendor).
  2. AnIML-lite XML export (Allotrope/AnIML is emerging as the open,
     vendor-neutral standard many chromatography data systems now support).
"""
import json
import requests
from xml.etree.ElementTree import Element, SubElement, tostring
from typing import Dict


def push_result_rest(lims_endpoint: str, api_key: str, run_result: Dict) -> requests.Response:
    payload = {
        "method": run_result["method"],
        "column": run_result["column"],
        "peaks": run_result["peaks"],
        "system_suitability": run_result.get("system_suitability"),
    }
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    return requests.post(lims_endpoint, headers=headers, data=json.dumps(payload), timeout=15)


def export_animl_lite(run_result: Dict) -> str:
    """Minimal AnIML-flavored XML for column/peak results (not full-schema-validated;
    intended as a starting point to map into your LIMS' actual AnIML importer)."""
    root = Element("AnIML")
    exp = SubElement(root, "Experiment")
    SubElement(exp, "Column", {
        "name": run_result["column"]["name"],
        "length_mm": str(run_result["column"]["length_mm"]),
        "id_mm": str(run_result["column"]["id_mm"]),
    })
    result_set = SubElement(exp, "PeakResults")
    for p in run_result["peaks"]:
        SubElement(result_set, "Peak", {
            "name": p["name"], "tR_min": str(p["tR"]),
            "area": str(p["area"]), "N": str(p["N"]),
        })
    return tostring(root, encoding="unicode")
