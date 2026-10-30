"""SMIRKS + spec parsing for auto_force_field (no antechamber)."""
import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem

from haddock.modules.refinement.openmm.auto_force_field import (
    cterm_smiles,
    load_spec,
    nterm_smiles,
    parse_cap,
    )


def _s(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)


def test_iva_smirks():
    capped = "CC(=O)N[C@@](C)(CC)C(=O)NC"
    assert _s(nterm_smiles(capped)) == _s("[NH3+][C@@](C)(CC)C(=O)NC")
    assert _s(cterm_smiles(capped)) == _s("CC(=O)N[C@@](C)(CC)C(=O)[O-]")


def test_ncap_has_no_nterm():
    assert nterm_smiles("CNC(=O)CC12CC3CC(CC(C3)C1)C2") is None


def test_ccap_has_no_cterm():
    assert cterm_smiles("CC(=O)NCC12CC3CC(CC(C3)C1)C2") is None


def test_parse_cap():
    assert parse_cap(False) is False
    assert parse_cap("N") == "N"
    assert parse_cap("c") == "C"
    with pytest.raises(ValueError):
        parse_cap(True)


def test_load_spec(tmp_path):
    p = tmp_path / "s.yaml"
    p.write_text(
        "IVA:\n  smiles: CC(=O)NCC(=O)NC\n  analog: VAL\n  cap: false\n"
        "ADA:\n  smiles: CNC(=O)CC12CC3CC(CC(C3)C1)C2\n  analog: ACE\n  cap: N\n"
        )
    spec = load_spec(p)
    assert spec["IVA"]["cap"] is False and spec["IVA"]["ace"] and spec["IVA"]["nme"]
    assert spec["ADA"]["cap"] == "N" and spec["ADA"]["ace"] is False and spec["ADA"]["nme"]
    p.write_text("chain: B\nligand: ADA_RB1_CTN\n" + p.read_text())
    spec2 = load_spec(p)
    assert "chain" not in spec2 and spec2["IVA"]["analog"] == "VAL"


def test_load_spec_rejects_ncap_with_ace(tmp_path):
    p = tmp_path / "s.yaml"
    p.write_text("ADA:\n  smiles: CC(=O)NCC(=O)NC\n  analog: ACE\n  cap: N\n")
    with pytest.raises(ValueError, match="N-cap"):
        load_spec(p)


def test_remap_cns_atom_name():
    pytest.importorskip("openmm")
    from haddock.modules.refinement.openmm.openmm import remap_cns_atom_name
    assert remap_cns_atom_name("HN") == "H"
    assert remap_cns_atom_name("OT1") == "O"
    assert remap_cns_atom_name("OT2") == "OXT"
    assert remap_cns_atom_name("C5") == "C5"
    assert remap_cns_atom_name("C5", linker=True) == "C"
    assert remap_cns_atom_name("O7", linker=True) == "OXT"
    assert remap_cns_atom_name("HT1") == "HT1"
    assert remap_cns_atom_name("HT1", nterm=True) == "H1"
    assert remap_cns_atom_name("HT3", nterm=True) == "H3"
    assert remap_cns_atom_name("CA") == "CA"


def test_resolve_platform_cpu():
    pytest.importorskip("openmm")
    from haddock.core.exceptions import ModuleError
    from haddock.modules.refinement.openmm.openmm import OPENMM

    job = OPENMM.__new__(OPENMM)
    job.params = {"platform": "CPU", "ncores": 2, "constraints": "HBonds"}
    job.log = type("L", (), {"info": lambda *a, **k: None})()
    platform, props = job.resolve_platform()
    assert platform.getName() == "CPU"
    assert props["Threads"] == "2"
    job.params["platform"] = "not-a-platform"
    with pytest.raises(ModuleError, match="not available"):
        job.resolve_platform()


def test_gpu_precision():
    openmm = pytest.importorskip("openmm")
    from openmm import Platform
    from haddock.modules.refinement.openmm.openmm import OPENMM

    names = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    gpu = next((n for n in ("CUDA", "OpenCL", "HIP") if n in names), None)
    if gpu is None:
        pytest.skip("no GPU OpenMM platform")
    job = OPENMM.__new__(OPENMM)
    job.params = {"platform": gpu, "precision": "double", "ncores": 1}
    job.log = type("L", (), {"info": lambda *a, **k: None})()
    _platform, props = job.resolve_platform()
    assert props["Precision"] == "double"
