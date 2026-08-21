#!/usr/bin/env python3
"""
Crystallography MCP server — CIF inspection, d-spacings, powder patterns, and
polymorph fingerprinting, backed by pymatgen.

Built for the Laurelate bay-nut-fat work but domain-general. The fat-specific
part is the alpha / beta-prime / beta short-spacing diagnostic, which is
standard lipid crystallography rather than anything proprietary.

Why this exists rather than a VESTA MCP: VESTA is a GUI application with no
documented command-line, batch, or scripting interface, so an MCP over it could
only open files and screenshot windows. pymatgen exposes the actual science.

Run:  ~/.venvs/crystal/bin/python ~/mcp-servers/crystallography/server.py
"""

import json
import math
import warnings
from pathlib import Path
from typing import Any

warnings.filterwarnings("ignore")

import numpy as np
from mcp.server.mcpserver import MCPServer

from pymatgen.io.cif import CifParser
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

CU_KALPHA = 1.5406

# Diagnostic subcell short spacings, Angstrom. Literature values for the three
# common triacylglycerol polymorphs. These are the fingerprints powder XRD
# resolves; the long spacing tells you chain packing, the short spacings tell
# you which form.
POLYMORPH_SHORT_A = {
    "alpha": [4.15],                       # single broad line, hexagonal subcell
    "beta_prime": [4.34, 4.11, 3.85, 3.80],  # ~4.2 + ~3.8 doublet, orthorhombic
    "beta": [4.60, 3.85, 3.70],            # strong 4.6 singlet, triclinic
}

MAX_SITES = 20000   # refuse absurd files rather than hanging

server = MCPServer(
    name="crystallography",
    version="0.1.0",
    instructions=(
        "Crystallography tools over pymatgen: read CIF files, compute "
        "d-spacings and powder diffraction patterns, analyse symmetry, compare "
        "structures, and fingerprint triacylglycerol polymorphs (alpha / "
        "beta-prime / beta) from either a structure or an experimental powder "
        "pattern.\n\n"
        "Note on multi-datablock CIFs: some deposited files hold more than one "
        "structure (e.g. van Langevelde 2000 contains both CLC and MPM). Always "
        "call list_datablocks first on an unfamiliar file, then pass the "
        "datablock index you want. Silently taking block 0 is a common way to "
        "analyse the wrong structure."
    ),
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _resolve(cif_path: str) -> Path:
    p = Path(cif_path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"no such file: {p}")
    if not p.is_file():
        raise ValueError(f"not a file: {p}")
    return p


def _load(cif_path: str, datablock: int = 0):
    p = _resolve(cif_path)
    structures = CifParser(str(p)).parse_structures(primitive=False)
    if not structures:
        raise ValueError(f"{p.name}: no structures parsed")
    if not 0 <= datablock < len(structures):
        raise IndexError(
            f"{p.name} has {len(structures)} datablock(s); index {datablock} "
            f"is out of range. Call list_datablocks first."
        )
    s = structures[datablock]
    if len(s) > MAX_SITES:
        raise ValueError(f"{p.name} block {datablock} has {len(s)} sites, "
                         f"above the {MAX_SITES} limit")
    return s, p


def _d_to_two_theta(d: float, wavelength: float) -> float:
    x = wavelength / (2.0 * d)
    if abs(x) >= 1.0:
        return float("nan")
    return math.degrees(2.0 * math.asin(x))


def _two_theta_to_d(tt: float, wavelength: float) -> float:
    s = math.sin(math.radians(tt) / 2.0)
    return float("inf") if s == 0 else wavelength / (2.0 * s)


def _ok(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, default=str)


# ---------------------------------------------------------------------------
# tools
# ---------------------------------------------------------------------------

@server.tool(
    description=(
        "List the datablocks (separate structures) inside a CIF file. Call this "
        "FIRST on any unfamiliar file — deposited CIFs often contain several "
        "structures, and analysing the wrong one is a silent error."
    )
)
def list_datablocks(cif_path: str) -> str:
    p = _resolve(cif_path)
    parser = CifParser(str(p))
    keys = [k for k in parser.as_dict().keys() if k != "global"]
    structures = parser.parse_structures(primitive=False)
    blocks = []
    for i, s in enumerate(structures):
        lat = s.lattice
        blocks.append({
            "index": i,
            "name": keys[i] if i < len(keys) else None,
            "reduced_formula": s.composition.reduced_formula,
            "n_sites": len(s),
            "a": round(lat.a, 4), "b": round(lat.b, 4), "c": round(lat.c, 4),
            "alpha": round(lat.alpha, 3), "beta": round(lat.beta, 3),
            "gamma": round(lat.gamma, 3),
            "volume_A3": round(lat.volume, 2),
        })
    return _ok({"file": p.name, "n_datablocks": len(structures), "datablocks": blocks})


@server.tool(
    description=(
        "Read one structure from a CIF: unit cell, volume, composition, site "
        "count, density, and the space group as re-derived from the coordinates "
        "by spglib (which can legitimately differ from what the CIF declares)."
    )
)
def read_structure(cif_path: str, datablock: int = 0, symprec: float = 0.1) -> str:
    s, p = _load(cif_path, datablock)
    lat = s.lattice
    out: dict[str, Any] = {
        "file": p.name,
        "datablock": datablock,
        "reduced_formula": s.composition.reduced_formula,
        "full_formula": s.composition.formula,
        "n_sites": len(s),
        "element_counts": {str(el): int(n) for el, n in s.composition.items()},
        "cell": {
            "a": round(lat.a, 5), "b": round(lat.b, 5), "c": round(lat.c, 5),
            "alpha": round(lat.alpha, 4), "beta": round(lat.beta, 4),
            "gamma": round(lat.gamma, 4),
            "volume_A3": round(lat.volume, 3),
        },
        "density_g_cm3": round(float(s.density), 4),
    }
    try:
        sga = SpacegroupAnalyzer(s, symprec=symprec)
        out["symmetry"] = {
            "spglib_symbol": sga.get_space_group_symbol(),
            "spglib_number": sga.get_space_group_number(),
            "crystal_system": sga.get_crystal_system(),
            "point_group": sga.get_point_group_symbol(),
            "symprec": symprec,
            "note": ("spglib re-derives symmetry from coordinates. Published "
                     "long-chain structures with disordered or split chain-end "
                     "sites often need a looser symprec to recover the "
                     "deposited space group."),
        }
    except Exception as exc:                                   # noqa: BLE001
        out["symmetry"] = {"error": f"{type(exc).__name__}: {exc}", "symprec": symprec}
    return _ok(out)


@server.tool(
    description=(
        "Compute interplanar d-spacings for specific reflections. Pass hkl as a "
        "list of [h,k,l] triples, e.g. [[2,0,0],[0,0,1]]. Use this for a known "
        "lamellar long spacing — for a triacylglycerol beta-prime cell where a "
        "is the stacking axis, d(200) is the long spacing."
    )
)
def d_spacings(cif_path: str, hkl: list[list[int]], datablock: int = 0,
               wavelength: float = CU_KALPHA) -> str:
    s, p = _load(cif_path, datablock)
    rows = []
    for triple in hkl:
        if len(triple) != 3:
            raise ValueError(f"hkl entries must be 3 integers, got {triple}")
        h, k, l = (int(v) for v in triple)
        if (h, k, l) == (0, 0, 0):
            raise ValueError("(000) is not a reflection")
        d = float(s.lattice.d_hkl((h, k, l)))
        rows.append({
            "hkl": [h, k, l],
            "d_A": round(d, 5),
            "two_theta_deg": round(_d_to_two_theta(d, wavelength), 4),
        })
    return _ok({"file": p.name, "datablock": datablock,
                "wavelength_A": wavelength, "reflections": rows})


@server.tool(
    description=(
        "Calculate the powder X-ray diffraction pattern from the atomic "
        "coordinates, using tabulated atomic form factors. Returns the "
        "strongest peaks with 2theta, d-spacing, relative intensity and hkl "
        "indices. Default wavelength is Cu K-alpha (1.5406 A)."
    )
)
def powder_pattern(cif_path: str, datablock: int = 0,
                   wavelength: float = CU_KALPHA,
                   two_theta_min: float = 1.5, two_theta_max: float = 40.0,
                   top_n: int = 25) -> str:
    s, p = _load(cif_path, datablock)
    calc = XRDCalculator(wavelength=wavelength)
    pat = calc.get_pattern(s, two_theta_range=(two_theta_min, two_theta_max))
    if len(pat.x) == 0:
        return _ok({"file": p.name, "peaks": [],
                    "note": "no reflections in the requested 2theta range"})
    imax = float(max(pat.y))
    peaks = []
    for tt, inten, hkls in zip(pat.x, pat.y, pat.hkls):
        idx = hkls[0]["hkl"] if hkls else None
        peaks.append({
            "two_theta_deg": round(float(tt), 4),
            "d_A": round(_two_theta_to_d(float(tt), wavelength), 5),
            "rel_intensity": round(float(inten) / imax, 5),
            "hkl": [int(v) for v in idx] if idx is not None else None,
            "n_overlapping": len(hkls),
        })
    peaks.sort(key=lambda r: -r["rel_intensity"])
    return _ok({
        "file": p.name, "datablock": datablock, "wavelength_A": wavelength,
        "two_theta_range": [two_theta_min, two_theta_max],
        "n_peaks_total": len(pat.x),
        "peaks": peaks[:top_n],
        "note": ("Intensities use tabulated form factors with a Debye-Waller "
                 "term omitted; positions are exact, relative intensities are "
                 "good but not Rietveld-grade."),
    })


@server.tool(
    description=(
        "Fingerprint the polymorph of a triacylglycerol STRUCTURE by examining "
        "its short-spacing (subcell) reflections in the 3.4-5.0 A window, where "
        "alpha / beta-prime / beta differ diagnostically. Reports which form the "
        "computed pattern matches and why."
    )
)
def fingerprint_structure(cif_path: str, datablock: int = 0,
                          wavelength: float = CU_KALPHA) -> str:
    s, p = _load(cif_path, datablock)
    calc = XRDCalculator(wavelength=wavelength)
    lo_tt = _d_to_two_theta(5.0, wavelength)
    hi_tt = _d_to_two_theta(3.4, wavelength)
    pat = calc.get_pattern(s, two_theta_range=(lo_tt, hi_tt))
    if len(pat.x) == 0:
        return _ok({"file": p.name, "verdict": "INDETERMINATE",
                    "reason": "no reflections in the 3.4-5.0 A window"})
    imax = float(max(pat.y))
    lines = sorted(
        ({"d_A": round(_two_theta_to_d(float(tt), wavelength), 4),
          "two_theta_deg": round(float(tt), 3),
          "rel_intensity": round(float(i) / imax, 4)}
         for tt, i in zip(pat.x, pat.y)),
        key=lambda r: -r["rel_intensity"],
    )
    return _ok({
        "file": p.name, "datablock": datablock,
        "short_spacing_lines": lines[:10],
        "assessment": _assess(lines),
        "reference_fingerprints": POLYMORPH_SHORT_A,
    })


def _assess(lines: list[dict]) -> dict[str, Any]:
    """Shared polymorph logic for both the structure and pattern entry points."""
    if not lines:
        return {"verdict": "INDETERMINATE", "reason": "no short-spacing lines"}
    strongest = lines[0]["rel_intensity"]
    strong = [r for r in lines if r["rel_intensity"] >= 0.35 * strongest]

    def near(target: float, tol: float = 0.09) -> list[dict]:
        return [r for r in strong if abs(r["d_A"] - target) <= tol]

    has_46 = bool(near(4.60))
    bprime_doublet = bool(near(4.20, 0.12)) and bool(near(3.80, 0.10))
    alpha_only = len(strong) == 1 and bool(near(4.15, 0.10))

    if bprime_doublet and not has_46:
        return {"verdict": "BETA_PRIME",
                "reason": ("resolvable ~4.2 + ~3.8 A doublet with no strong "
                           "4.6 A line — orthorhombic perpendicular subcell")}
    if has_46 and not bprime_doublet:
        return {"verdict": "BETA",
                "reason": ("strong ~4.6 A singlet without the doublet — "
                           "triclinic parallel subcell")}
    if alpha_only:
        return {"verdict": "ALPHA",
                "reason": ("a single broad line near 4.15 A — hexagonal "
                           "subcell, the transient form. Check thermal history "
                           "before drawing equilibrium conclusions.")}
    if has_46 and bprime_doublet:
        return {"verdict": "MIXED",
                "reason": ("both a 4.6 A line and a 4.2/3.8 A doublet — a "
                           "beta + beta-prime mixture")}
    return {"verdict": "INDETERMINATE",
            "reason": "short-spacing lines match no single fingerprint cleanly"}


@server.tool(
    description=(
        "Fingerprint the polymorph from an EXPERIMENTAL powder pattern. Pass "
        "two_theta and intensity as equal-length lists, or a path to a "
        "two-column text/CSV file via pattern_path. Peak-picks, then matches "
        "the short spacings against the alpha / beta-prime / beta diagnostics. "
        "This is the tool for adjudicating a measured scan."
    )
)
def fingerprint_pattern(two_theta: list[float] | None = None,
                        intensity: list[float] | None = None,
                        pattern_path: str | None = None,
                        wavelength: float = CU_KALPHA,
                        min_prominence: float = 0.05) -> str:
    if pattern_path:
        tt_list, i_list = [], []
        with open(_resolve(pattern_path), encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith(("#", "'", '"', ";", "*")):
                    continue
                parts = [x for x in line.replace(",", " ").split() if x]
                if len(parts) < 2:
                    continue
                try:
                    tt_list.append(float(parts[0]))
                    i_list.append(float(parts[1]))
                except ValueError:
                    continue
        tt_arr, i_arr = np.array(tt_list), np.array(i_list)
    elif two_theta is not None and intensity is not None:
        if len(two_theta) != len(intensity):
            raise ValueError(f"two_theta has {len(two_theta)} points but "
                             f"intensity has {len(intensity)}")
        tt_arr, i_arr = np.array(two_theta, float), np.array(intensity, float)
    else:
        raise ValueError("provide either pattern_path, or both two_theta and intensity")

    if len(tt_arr) < 10:
        raise ValueError(f"need at least 10 data points, got {len(tt_arr)}")

    order = np.argsort(tt_arr)
    tt_arr, i_arr = tt_arr[order], i_arr[order]
    span = float(np.ptp(i_arr)) or 1.0
    norm = (i_arr - float(np.min(i_arr))) / span

    # Local-maximum peak pick with a flat-baseline prominence filter. Simple on
    # purpose: this triages which polymorph family a pattern belongs to, it is
    # not a Rietveld refinement.
    window = max(3, len(tt_arr) // 200)
    picked = []
    for i in range(window, len(norm) - window):
        seg = norm[i - window:i + window + 1]
        if norm[i] == seg.max() and norm[i] >= min_prominence:
            if norm[i] - min(seg[0], seg[-1]) >= min_prominence:
                picked.append({
                    "two_theta_deg": round(float(tt_arr[i]), 3),
                    "d_A": round(_two_theta_to_d(float(tt_arr[i]), wavelength), 4),
                    "rel_intensity": round(float(norm[i]), 4),
                })
    picked.sort(key=lambda r: -r["rel_intensity"])

    # collapse peaks within 0.05 A, keeping the stronger
    kept: list[dict] = []
    for pk in picked:
        if not any(abs(pk["d_A"] - q["d_A"]) < 0.05 for q in kept):
            kept.append(pk)

    short = [r for r in kept if 3.4 <= r["d_A"] <= 5.0]
    long_lines = [r for r in kept if r["d_A"] > 20.0]
    return _ok({
        "n_points": int(len(tt_arr)),
        "two_theta_range": [round(float(tt_arr.min()), 3), round(float(tt_arr.max()), 3)],
        "n_peaks_picked": len(kept),
        "short_spacing_lines": short[:10],
        "long_spacing_lines": long_lines[:5],
        "assessment": _assess(short),
        "reference_fingerprints": POLYMORPH_SHORT_A,
        "caveat": ("Simple local-maximum peak picking. For a weak or noisy scan, "
                   "lower min_prominence and sanity-check the picked list before "
                   "trusting the verdict."),
    })


@server.tool(
    description=(
        "Compare two structures for equivalence using pymatgen's "
        "StructureMatcher, reporting whether they match and the RMS "
        "displacement. Useful for checking whether a candidate model reproduces "
        "a solved template."
    )
)
def compare_structures(cif_path_a: str, cif_path_b: str,
                       datablock_a: int = 0, datablock_b: int = 0,
                       ltol: float = 0.2, stol: float = 0.3,
                       angle_tol: float = 5.0) -> str:
    sa, pa = _load(cif_path_a, datablock_a)
    sb, pb = _load(cif_path_b, datablock_b)
    matcher = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol)
    matched = bool(matcher.fit(sa, sb))
    rms = matcher.get_rms_dist(sa, sb)
    return _ok({
        "a": {"file": pa.name, "datablock": datablock_a,
              "formula": sa.composition.reduced_formula, "n_sites": len(sa)},
        "b": {"file": pb.name, "datablock": datablock_b,
              "formula": sb.composition.reduced_formula, "n_sites": len(sb)},
        "match": matched,
        "rms_dist": None if rms is None else [round(float(v), 5) for v in rms],
        "tolerances": {"ltol": ltol, "stol": stol, "angle_tol": angle_tol},
        "note": ("A False match between different chain lengths is expected — "
                 "StructureMatcher compares composition and geometry, so a C16 "
                 "template will never match a C10/C12 fat."),
    })


if __name__ == "__main__":
    server.run("stdio")
