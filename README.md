# Crystallography MCP server

CIF inspection, d-spacings, powder patterns, and triacylglycerol polymorph
fingerprinting, backed by [pymatgen](https://pymatgen.org/).

Built 21 August 2026 for the Laurelate bay-nut-fat work, but domain-general.

## Why this and not a VESTA MCP

VESTA has **no scripting interface** — no documented command-line arguments, no
batch mode, no macros, no headless operation. An MCP over it could only run
`open -a VESTA file.cif` and screenshot the window. pymatgen exposes the actual
science, so that is what this wraps.

Keep VESTA as a human viewer. It is genuinely good at that.

## Install / registration

Registered at **user** scope, so it loads in every Claude Code session:

```bash
claude mcp add -s user crystallography -- /Users/laurelate/.venvs/crystal/bin/python /Users/laurelate/mcp-servers/crystallography/server.py
```

Runs on its own venv (`~/.venvs/crystal`, Python 3.13) because current pymatgen
requires Python ≥3.10 and the system Python is 3.9.6. That venv is deliberately
separate from anything the R&D repo uses.

Check it:

```bash
claude mcp list
```

Remove it:

```bash
claude mcp remove -s user crystallography
```

## Tools

| Tool | What it does |
|---|---|
| `list_datablocks` | Enumerate the structures inside a CIF. **Call this first on an unfamiliar file.** |
| `read_structure` | Cell, volume, composition, density, spglib space group |
| `d_spacings` | d and 2θ for specific `[h,k,l]` reflections |
| `powder_pattern` | Full calculated pattern with hkl indices and relative intensities |
| `fingerprint_structure` | α / β′ / β verdict from a *structure's* short spacings |
| `fingerprint_pattern` | α / β′ / β verdict from an *experimental* scan (list or two-column file) |
| `compare_structures` | pymatgen `StructureMatcher` fit + RMS displacement |

### The multi-datablock trap

Some deposited CIFs hold more than one structure. `vanlangevelde2000_CLC_MPM_bprime.cif`
contains **both** CLC (`a = 57.368 Å`, the β′ template that matters) and MPM
(`a = 76.21 Å`). Taking block 0 silently is how you end up analysing the wrong
structure. `list_datablocks` exists to make that visible.

## Polymorph fingerprints

Diagnostic subcell short spacings, Å:

| Form | Lines | Signature |
|---|---|---|
| α | 4.15 | single broad line, hexagonal subcell |
| β′ | 4.34 / 4.11 + 3.85 / 3.80 | **doublet**, orthorhombic perpendicular |
| β | 4.60, 3.85, 3.70 | strong **4.6 singlet**, triclinic parallel |

Validated against both solved templates, which have known answers:

- `vanlangevelde2000_CLC_MPM_bprime.cif` block 0 → **BETA_PRIME** ✓
- `vanlangevelde1999_PPP_beta.cif` → **BETA** ✓

And against synthetic patterns built from each reference line set → correct
verdict in all three cases.

## Known limits

- **Intensities are not Rietveld-grade.** Tabulated form factors, Debye–Waller
  omitted. Peak *positions* are exact; relative intensities are good enough to
  rank lines, not to refine occupancies.
- **Peak picking is deliberately simple** — local maxima with a flat-baseline
  prominence filter. For a weak or noisy scan, lower `min_prominence` and eyeball
  the picked list before trusting the verdict.
- **spglib may disagree with the CIF's declared space group.** It re-derives
  symmetry from coordinates, and published long-chain structures often carry
  disordered or split chain-end sites. Loosen `symprec` before concluding
  anything.
- **`compare_structures` returning False across different chain lengths is
  expected**, not a bug — a C16 template will never match a C10/C12 fat.

## Provenance note

This server found a real bug in `R&D/modeling/baynutfat/structure.py`: the
I-centring translation was being applied twice for β′-CLC (the CIF already lists
all 8 centred operators), leaving 4 duplicated carbons — 284 C where there
should be 280. Peak positions were unaffected, structure factors were slightly
off. Guarded now by acceptance test T10 and `modeling/tests/crosscheck_pymatgen.py`.

That is the argument for having a second implementation available.
