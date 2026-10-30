"""GAFF extra-FF XMLs from a HADDOCK PDB + a residue spec.

The spec names every residue to build. SMILES is the calculation molecule
(ACE–X–NME, or one-sided for a cap) and supplies bonds. `analog` is the
ff14SB residue used for N/C-term backbone charges. `cap` is false, N, or C.

  ACE  [CH3][C](=O)[N]  >> [NH3+]     N-term
  NME  [C](=O)[NH][CH3] >> [C](=O)[O-] C-term
  those atoms are dropped for ExternalBonds.
"""
from __future__ import annotations

import argparse
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("AMBERHOME", "/home/a/anaconda3/envs/rad51_protac")
os.environ["PATH"] = os.environ["AMBERHOME"] + "/bin:" + os.environ.get("PATH", "")

from lxml import etree
from rdkit import Chem
from rdkit.Chem import AllChem

BBQ = {"ns": -0.4157, "hn": 0.2719, "c": 0.5973, "o": -0.5679}
NTERM = {
    "VAL": {"N": 0.0577, "N1": 0.0577, "H1": 0.2272, "H2": 0.2272, "H3": 0.2272, "c": 0.6163, "o": -0.5722},
    "GLY": {"N": 0.2943, "N1": 0.2943, "H1": 0.1642, "H2": 0.1642, "H3": 0.1642, "c": 0.6163, "o": -0.5722},
    "SER": {"N": 0.1849, "N1": 0.1849, "H1": 0.1898, "H2": 0.1898, "H3": 0.1898, "c": 0.6163, "o": -0.5722},
}
CTERM = {
    "VAL": {"N": -0.3821, "N1": -0.3821, "H": 0.2681, "H1": 0.2719,
            "C": 0.835, "C3": 0.835, "C4": 0.835, "C5": 0.835,
            "O": -0.8173, "O1": -0.8173, "O3": -0.8173, "OXT": -0.8173},
    "GLY": {"N": -0.3821, "N1": -0.3821, "H": 0.2681, "H1": 0.2719,
            "C": 0.7231, "C3": 0.7231, "C4": 0.7231,
            "O": -0.7855, "O1": -0.7855, "O3": -0.7855, "OXT": -0.7855},
    "SER": {"N": -0.3821, "N1": -0.3821, "H": 0.2681, "H1": 0.2719,
            "C": 0.8113, "C3": 0.8113, "C4": 0.8113, "C5": 0.8113,
            "O": -0.8132, "O1": -0.8132, "O3": -0.8132, "OXT": -0.8132},
}
ELEM = {"c": "C", "c3": "C", "c6": "C", "o": "O", "os": "O", "ns": "N", "n4": "N", "n": "N",
        "hn": "H", "hx": "H", "hc": "H", "h1": "H"}
ACE_TO_NH3 = AllChem.ReactionFromSmarts("[CH3][C](=O)[N:1]>>[NH3+:1]")
NME_TO_COO = AllChem.ReactionFromSmarts("[C:1](=[O:2])[NH][CH3]>>[C:1](=[O:2])[O-]")
ACE_Q = Chem.MolFromSmarts("[CH3][C](=O)[N]")
NME_Q = Chem.MolFromSmarts("[C](=O)[NH][CH3]")


def _rxn(mol, rxn):
    hits = rxn.RunReactants((mol,))
    if not hits:
        return mol
    out = hits[0][0]
    try:
        Chem.SanitizeMol(out)
    except Exception:
        pass
    return out


def canonical(mol):
    return Chem.MolToSmiles(mol, canonical=True)


def nterm_smiles(capped):
    mol = Chem.MolFromSmiles(capped)
    Chem.SanitizeMol(mol)
    if not mol.HasSubstructMatch(ACE_Q):
        return None
    return canonical(_rxn(mol, ACE_TO_NH3))


def cterm_smiles(capped):
    mol = Chem.MolFromSmiles(capped)
    Chem.SanitizeMol(mol)
    if not mol.HasSubstructMatch(NME_Q):
        return None
    return canonical(_rxn(mol, NME_TO_COO))


def parse_cap(value):
    if value in (False, None, "", 0):
        return False
    s = str(value).strip().upper()
    if s in ("FALSE", "NO", "NONE", "0"):
        return False
    if s in ("N", "N-CAP", "NCAP"):
        return "N"
    if s in ("C", "C-CAP", "CCAP"):
        return "C"
    raise ValueError(f"cap must be false, N, or C; got {value!r}")


def load_spec(path):
    import yaml
    raw = yaml.safe_load(Path(path).read_text()) or {}
    if not raw:
        raise ValueError("empty auto_force_field spec")
    out = {}
    for name, rec in raw.items():
        if name in ("chain", "ligand") or not isinstance(rec, dict):
            continue
        if not {"smiles", "analog", "cap"} <= rec.keys():
            raise ValueError(f"{name} needs smiles, analog, cap")
        cap = parse_cap(rec["cap"])
        analog = str(rec["analog"]).upper()
        if cap is False and analog not in NTERM:
            raise ValueError(f"{name} analog {analog} not in {sorted(NTERM)}")
        smiles = rec["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"{name} bad SMILES {smiles}")
        Chem.SanitizeMol(mol)
        has_ace = bool(mol.HasSubstructMatch(ACE_Q))
        has_nme = bool(mol.HasSubstructMatch(NME_Q))
        if cap is False and not (has_ace and has_nme):
            raise ValueError(f"{name} cap=false SMILES must be ACE-X-NME")
        if cap == "N" and has_ace:
            raise ValueError(f"{name} N-cap SMILES is X-NHMe, not ACE-X-NME")
        if cap == "C" and has_nme:
            raise ValueError(f"{name} C-cap SMILES is ACE-X, not ACE-X-NME")
        out[name] = {
            "smiles": smiles,
            "analog": analog,
            "cap": cap,
            "ace": cap != "N",
            "nme": cap != "C",
            }
    return out


def _hs(rd, idx):
    return [n.GetIdx() for n in rd.GetAtomWithIdx(idx).GetNeighbors() if n.GetAtomicNum() == 1]


def drop_cap_indices(rd, ace=True, nme=True):
    drop = set()
    if ace:
        hit = rd.GetSubstructMatch(ACE_Q)
        if hit:
            me, c, o, _n = hit
            drop.update((me, c, o), _hs(rd, me))
    if nme:
        hit = rd.GetSubstructMatch(NME_Q)
        if hit:
            _c, _o, n, me = hit
            drop.update((n, me), _hs(rd, n), _hs(rd, me))
    return drop


def gaff_xml_path():
    from importlib.resources import files
    return str(files("openmmforcefields") / "ffxml/amber/gaff/ffxml/gaff-2.11.xml")


def _cns_top():
    try:
        from haddock import toppar_path
        return Path(toppar_path) / "protein-allhdg5-4.top"
    except Exception:
        return Path("/home/a/haddock3/src/haddock/cns/toppar/protein-allhdg5-4.top")


def cns_bonds(resname):
    grab, bonds = False, []
    path = _cns_top()
    if not path.is_file():
        return bonds
    for line in path.read_text().splitlines():
        if line.startswith(f"residue {resname} ") or line.strip() == f"residue {resname}":
            grab = True
            continue
        if grab and (line.strip() == "end" or line.startswith("residue ") or line.startswith("presidue ")):
            break
        if grab and "BOND" in line:
            toks = line.replace("BOND", " ").split()
            for i in range(0, len(toks) - 1, 2):
                a, b = toks[i], toks[i + 1]
                bonds.append(("H" if a == "HN" else a, "H" if b == "HN" else b))
    return bonds


def parse_pdb(path):
    from haddock.modules.refinement.openmm.openmm import LINKER_RES, remap_cns_atom_name
    chains = defaultdict(list)
    cur = None
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            if line.startswith("TER") and cur:
                cur = None
            continue
        name = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip() or " "
        resid = line[22:26].strip()
        nterm = cur is None or (cur[0] != chain)
        name = remap_cns_atom_name(
            name, nterm=nterm, linker=resn in LINKER_RES,
            )
        elem = (name[0] if name[0].isalpha() else name[1]).upper()
        xyz = tuple(float(line[30 + 8 * i:38 + 8 * i]) for i in range(3))
        key = (chain, resid, resn)
        if cur != key:
            chains[chain].append({"name": resn, "atoms": []})
            cur = key
        chains[chain][-1]["atoms"].append((name, elem, xyz))
    return chains


def residue_graph(res):
    atoms = res["atoms"]
    g, elem = defaultdict(set), {}
    for n, e, _ in atoms:
        elem[n] = e
    names = set(elem)
    cbonds = cns_bonds(res["name"])
    if cbonds:
        for a, b in cbonds:
            a, b = ("H" if a == "HN" else a), ("H" if b == "HN" else b)
            if a in names and b in names:
                g[a].add(b)
                g[b].add(a)
        return g, elem
    for i, (n1, e1, x) in enumerate(atoms):
        for n2, e2, y in atoms[i + 1:]:
            d = sum((x0 - y0) ** 2 for x0, y0 in zip(x, y)) ** 0.5
            if d < (1.2 if "H" in (e1, e2) else 1.8):
                g[n1].add(n2)
                g[n2].add(n1)
    return g, elem


def strip_leaving(atoms, g, elem):
    drop = set()
    for n, e in elem.items():
        if e != "C":
            continue
        oxy = [x for x in g[n] if elem[x] == "O"]
        if len(oxy) == 2:
            for o in oxy:
                hs = [x for x in g[o] if elem[x] == "H"]
                if hs:
                    drop.add(o)
                    drop.update(hs)
    for n, e in elem.items():
        if e != "N":
            continue
        hs = [x for x in g[n] if elem[x] == "H"]
        if len(hs) >= 2:
            drop.add(sorted(hs)[-1])
    atoms = [a for a in atoms if a[0] not in drop]
    g2 = defaultdict(set)
    for a, nbrs in g.items():
        if a in drop:
            continue
        g2[a] = {b for b in nbrs if b not in drop}
    elem = {k: v for k, v in elem.items() if k not in drop}
    return atoms, g2, elem


def pdb_residues(pdb, names):
    found = {}
    for residues in parse_pdb(pdb).values():
        for res in residues:
            name = res["name"]
            if name in names and name not in found:
                g, elem = residue_graph(res)
                _atoms, g, elem = strip_leaving(res["atoms"], g, elem)
                found[name] = (g, elem)
    missing = [n for n in names if n not in found]
    if missing:
        print("skip not in PDB:", missing)
    return found


def xml_graph(res):
    g = defaultdict(set)
    for b in res.findall("Bond"):
        g[b.get("atomName1")].add(b.get("atomName2"))
        g[b.get("atomName2")].add(b.get("atomName1"))
    elem = {a.get("name"): ELEM.get(a.get("type"), a.get("type")[0].upper())
            for a in res.findall("Atom")}
    return g, elem


def match(g_xml, e_xml, g_pdb, e_pdb, kind):
    h_xml = sorted(n for n in e_xml if e_xml[n] != "H")
    h_pdb = [n for n in e_pdb if e_pdb[n] != "H"]
    m, used = {}, set()

    def rec(i, xmls):
        if i == len(xmls):
            return True
        a = xmls[i]
        for b in h_pdb:
            if b in used or e_xml[a] != e_pdb[b]:
                continue
            if any((m[x] in g_pdb[b]) != (x in g_xml[a]) for x in m if e_xml[x] != "H"):
                continue
            m[a] = b
            used.add(b)
            if rec(i + 1, xmls):
                return True
            del m[a]
            used.remove(b)
        return False

    if not rec(0, h_xml):
        skip = []
        if sum(e_xml[n] == "N" for n in h_xml) > sum(e_pdb[n] == "N" for n in h_pdb):
            skip += [n for n in h_xml if e_xml[n] == "N"
                     and sum(e_xml[x] != "H" for x in g_xml[n]) == 1]
        if sum(e_xml[n] == "O" for n in h_xml) > sum(e_pdb[n] == "O" for n in h_pdb):
            extra_o = [n for n in h_xml if e_xml[n] == "O"]
            skip += extra_o[-1:]
        m, used = {}, set()
        if not skip or not rec(0, [n for n in h_xml if n not in skip]):
            raise RuntimeError("heavy-atom match failed")
    for a, b in list(m.items()):
        for x, y in zip(sorted(n for n in g_xml[a] if e_xml[n] == "H"),
                        sorted(n for n in g_pdb[b] if e_pdb[n] == "H")):
            m[x] = y
    unmapped = [n for n in e_xml if n not in m]
    extra_o = [n for n in unmapped if e_xml[n] == "O"]
    extra_n = [n for n in unmapped if e_xml[n] == "N"]
    extra_h = [n for n in unmapped if e_xml[n] == "H"]
    if extra_o:
        m[extra_o[0]] = "OXT"
        extra_h = [n for n in extra_h if n not in g_xml[extra_o[0]]]
    if extra_n:
        m[extra_n[0]] = "NT"
        for x, y in zip(sorted(n for n in g_xml[extra_n[0]] if e_xml[n] == "H"), ("HNT1", "HNT2")):
            m[x] = y
        extra_h = [n for n in extra_h if n not in m]
    if kind.startswith("N"):
        n_xml = next(a for a, b in m.items() if b in ("N", "N1") and e_xml[a] == "N")
        for x, y in zip(sorted(n for n in g_xml[n_xml] if e_xml[n] == "H"), ("H1", "H2", "H3")):
            m[x] = y
    leftover = [n for n in e_xml if n not in m]
    if leftover:
        if all(e_xml[n] in "NH" for n in leftover):
            for n in leftover:
                m[n] = n
        else:
            raise RuntimeError(f"unmapped {leftover}")
    return m


def set_charges(res, kind, target, analog_name):
    table = None
    if kind.startswith("N") and len(kind) > 1 and analog_name in NTERM:
        table = NTERM[analog_name]
    elif kind.startswith("C") and len(kind) > 1 and analog_name in CTERM:
        table = CTERM[analog_name]
    for atom in res.findall("Atom"):
        n, t = atom.get("name"), atom.get("type")
        if table and n in table:
            atom.set("charge", f"{table[n]:.6f}")
        elif table and t in table:
            atom.set("charge", f"{table[t]:.6f}")
        elif t in BBQ:
            atom.set("charge", f"{BBQ[t]:.6f}")
    keyed = set(table or ())
    side = [a for a in res.findall("Atom")
            if a.get("name") not in keyed and a.get("type") not in (table or BBQ)]
    leftover = target - sum(float(a.get("charge")) for a in res.findall("Atom"))
    if side:
        for a in side:
            a.set("charge", f"{float(a.get('charge')) + leftover / len(side):.6f}")


def _write_residue(name, smiles, ace, nme, charge, analog_name, pdb_g, pdb_e, out_res_parent):
    from openff.toolkit import Molecule
    from openmmforcefields.generators import GAFFTemplateGenerator

    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    rd = mol.to_rdkit()
    drop = drop_cap_indices(rd, ace=ace, nme=nme)
    kept = [a for i, a in enumerate(mol.atoms) if i not in drop]
    gen = GAFFTemplateGenerator(forcefield="gaff-2.11")
    root = etree.fromstring(gen.generate_residue_template(mol, residue_atoms=kept).encode())
    res = root.find(".//Residues/Residue")
    res.set("name", name)
    keep = {a.name for a in kept}
    for el in list(res):
        if el.tag == "Atom" and el.get("name") not in keep:
            res.remove(el)
        elif el.tag == "Bond" and (el.get("atomName1") not in keep or el.get("atomName2") not in keep):
            res.remove(el)
    for atom in res.findall("Atom"):
        if atom.get("type") == "c6":
            atom.set("type", "c3")
    new = match(*xml_graph(res), pdb_g, pdb_e, name)
    for atom in res.findall("Atom"):
        atom.set("name", new[atom.get("name")])
    for el in res.findall("Bond") + res.findall("ExternalBond"):
        for k in ("atomName", "atomName1", "atomName2"):
            if el.get(k) in new:
                el.set(k, new[el.get(k)])
    set_charges(res, name, charge, analog_name)
    out_res_parent.append(res)
    return res


def build_from_pdb(pdb, outdir, spec):
    spec = load_spec(spec) if not isinstance(spec, dict) else spec
    pdb = Path(pdb)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    graphs = pdb_residues(pdb, spec)
    ff = etree.Element("ForceField")
    residues_el = etree.SubElement(ff, "Residues")
    term_res = []
    for name, rec in spec.items():
        if name not in graphs:
            continue
        g, elem = graphs[name]
        print(f"{name} smiles {rec['smiles']} analog {rec['analog']} cap {rec['cap']}")
        _write_residue(name, rec["smiles"], rec["ace"], rec["nme"], 0.0,
                       rec["analog"], g, elem, residues_el)
        if rec["cap"] is False:
            term_res.append(name)
            nt, ct = nterm_smiles(rec["smiles"]), cterm_smiles(rec["smiles"])
            print(f"  N{name} {nt}")
            print(f"  C{name} {ct}")
            _write_residue("N" + name, nt, False, True, 1.0,
                           rec["analog"], g, elem, residues_el)
            _write_residue("C" + name, ct, True, False, -1.0,
                           rec["analog"], g, elem, residues_el)
    path = outdir / "auto_force_field.xml"
    path.write_text(etree.tostring(ff, pretty_print=True).decode())
    return path, term_res


def main(argv=None):
    p = argparse.ArgumentParser(description="GAFF XMLs from HADDOCK PDB + residue spec")
    p.add_argument("pdb")
    p.add_argument("--spec", required=True, help="YAML: resname: {smiles, analog, cap}")
    p.add_argument("-o", "--outdir", default="openmm_ff")
    args = p.parse_args(argv)
    path, term = build_from_pdb(args.pdb, args.outdir, args.spec)
    print("wrote", path, "termini", term)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
