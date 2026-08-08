"""
PubChem PUG-REST client (free, no API key).
Docs: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest

Used to auto-populate Compound.logP and molecular weight, which feed
a simple QSAR-style estimate of k_w / S when no experimental value exists.
"""
import requests

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def fetch_cid_by_name(name: str) -> int | None:
    url = f"{BASE}/compound/name/{name}/cids/JSON"
    r = requests.get(url, timeout=10)
    if r.status_code != 200:
        return None
    return r.json().get("IdentifierList", {}).get("CID", [None])[0]


def fetch_properties(name: str) -> dict | None:
    """Returns XLogP, molecular weight, TPSA, H-bond donors/acceptors for a compound name."""
    cid = fetch_cid_by_name(name)
    if cid is None:
        return None
    props = "XLogP,MolecularWeight,TPSA,HBondDonorCount,HBondAcceptorCount,CanonicalSMILES"
    url = f"{BASE}/compound/cid/{cid}/property/{props}/JSON"
    r = requests.get(url, timeout=10)
    if r.status_code != 200:
        return None
    result = r.json()["PropertyTable"]["Properties"][0]
    result["CID"] = cid
    return result


def estimate_rp_params_from_logp(xlogp: float) -> dict:
    """
    Very rough QSAR heuristic (linear free-energy relationship) mapping
    XLogP -> starting k_w / S for reverse-phase C18, calibrated to
    typical literature ranges. Meant as a *seed* value for the simulator,
    to be refined by Retip/RT-Transformer or real calibration data.
    """
    k_w = max(0.5, 10 ** (0.4 * xlogp + 0.3))
    S = max(1.5, min(7.0, 2.5 + 0.35 * xlogp))
    return {"k_w": round(k_w, 3), "S": round(S, 3)}
