# Higginson Lab fork — toppar changes (POM / VHP)

Fork-specific CNS topology/parameter work for peptide–PROTAC docking with **pomalidomide (POM)** and **VH032 (VHP)** C-terminal warheads. Upstream HADDOCK3 docs still apply; this file only records local toppar deltas.

## Files touched

| Path | Role |
|------|------|
| `src/haddock/cns/toppar/protein-allhdg5-4.top` | `residue POM`, `presidue PPOM`; `presidue PVHL` fix |
| `src/haddock/cns/toppar/protein-allhdg5-4.param` | POM bonded + AHX→NAT junction; VHP CTER `CAX`–`OC` |
| `src/haddock/cns/toppar/protein-allhdg5-4*.link` | `link PPOM … + POM`; CTER policy for POM/VHP |
| `src/haddock/cns/toppar/pom.top` / `pom.param` | Standalone POM copies (also embedded in protein toppar) |

## POM residue + PPOM presidue

- **Source**: ATB molid 927298.
- **Linker attachment**: atom **NAT** (ATB **N10**), with amide H **HN** and patch-deleted H **H04**.
- **Chemical types for ring carbons**: atom names `C5`/`C9` use types **`C51`/`C91`** so dihedrals in param do not collide with atom-name ambiguity.
- **`presidue PPOM`**: peptide-like link previous `C` → `+NAT`; deletes `+H04`; angles/dihedrals through `+C6` / `+C5`/`+C9` / `+HN` (mirrors PVHL pattern for VHP).

Standalone `pom.top` / `pom.param` mirror the embedded POM block for reference or ligand-only use.

## Param additions (`protein-allhdg5-4.param`)

- Full POM ATB bonded set (bonds/angles/dihedrals/impropers/nonbonded) using **`C51`/`C91`**.
- Junction (AHX/peptide `C` → POM `NAT`): `BOND C NAT`, omega-like dihedrals, amide impropers.
- **VHP charged CTER**: `BOND CAX OC` and angles `CBE–CAX–OC`, `OC–CAX–OC`, `N1–CAX–OC` (VHP carbonyl atom `C` has `TYPE=CAX`).

## Link files

- **`link PPOM head - <AA> tail + POM`** added in all four protein link variants (POM only as C-terminal/tail partner).
- **No `last CTER` for POM** in any link file (warhead is not a standard carboxylate C-term).
- **No `last CTER` for VHP** in `protein-allhdg5-4-noter.link` and `protein-allhdg5-4-noCter.link` (used when topoaa sets `charged_cter = false` on the PROTAC). VHP CTER remains in the full `.link` and `-noNter.link` for charged-C-term cases.

## PVHL fix (VHP junction)

`presidue PVHL` updated for CNS parsing and peptide atom names already present on **VHP**:

- Trailing **`end`** on `delete ATOM +H25` (without it CNS stays in `ATOM>` and mis-parses).
- Bond/angles/dihedrals use atoms **`N` / `CA` / `HN` / `C` / `CB`** (`N` has `TYPE=NAT`).

The **VHP residue block itself was not rewritten** for the POM work; only PVHL and the CAX–OC param support were adjusted.

## Operational note

PROTAC binary docks that use POM/VHP typically set mol2 `charged_nter = false` and `charged_cter = false` so the noter/noCter link sets apply (no CTER patch on the warhead).
