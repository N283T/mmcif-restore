"""Test fixtures for mmcif-sync tests."""

import os
import tempfile
from collections.abc import Generator
from pathlib import Path

import gemmi
import pytest


@pytest.fixture
def pdb_mmcif_path() -> Path | None:
    """Get PDB mmCIF path from environment."""
    path = os.environ.get("PDB_MMCIF_PATH")
    if path:
        return Path(path)
    return None


@pytest.fixture
def sample_cif_content() -> str:
    """Minimal valid mmCIF content for testing."""
    return """\
data_TEST
#
_entry.id TEST
#
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'Test protein'
2 non-polymer 'Ligand'
3 water 'water'
#
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_strand_id
1 'polypeptide(L)' A
#
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
C 3
#
loop_
_pdbx_nonpoly_scheme.asym_id
_pdbx_nonpoly_scheme.entity_id
_pdbx_nonpoly_scheme.mon_id
_pdbx_nonpoly_scheme.pdb_seq_num
_pdbx_nonpoly_scheme.auth_seq_num
_pdbx_nonpoly_scheme.pdb_mon_id
_pdbx_nonpoly_scheme.auth_mon_id
_pdbx_nonpoly_scheme.pdb_strand_id
_pdbx_nonpoly_scheme.auth_asym_id
B 2 LIG 1 1 LIG LIG B B
C 3 HOH 1 1 HOH HOH C C
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 1 1 ? 0.0 0.0 0.0 1.0 10.0 1 ALA A N 1
ATOM 2 CA CA . ALA A 1 1 ? 1.5 0.0 0.0 1.0 10.0 1 ALA A CA 1
HETATM 3 C1 C1 . LIG B 2 . ? 5.0 5.0 5.0 1.0 15.0 1 LIG B C1 1
HETATM 4 O O . HOH C 3 . ? 10.0 10.0 10.0 1.0 20.0 1 HOH C O 1
#
"""


@pytest.fixture
def sample_cif_file(sample_cif_content: str) -> Generator[Path, None, None]:
    """Create a temporary CIF file for testing."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
        f.write(sample_cif_content)
        f.flush()
        yield Path(f.name)
    os.unlink(f.name)


@pytest.fixture
def sample_structure(sample_cif_file: Path) -> gemmi.Structure:
    """Load sample structure for testing."""
    return gemmi.read_structure(str(sample_cif_file))


@pytest.fixture
def sample_cif_document(sample_cif_file: Path) -> gemmi.cif.Document:
    """Load sample CIF document for testing."""
    return gemmi.cif.read(str(sample_cif_file))


def get_real_cif_path(pdb_id: str, pdb_mmcif_path: Path | None) -> Path | None:
    """Get path to real PDB mmCIF file if available."""
    if pdb_mmcif_path is None:
        return None

    # PDB files are stored in subdirectories by middle two characters
    # e.g., 1abc -> ab/1abc.cif.gz
    middle = pdb_id[1:3].lower()
    cif_path = pdb_mmcif_path / middle / f"{pdb_id.lower()}.cif.gz"

    if cif_path.exists():
        return cif_path
    return None


@pytest.fixture
def real_cif_2y2a(pdb_mmcif_path: Path | None) -> Path | None:
    """Get path to 2y2a.cif.gz for integration testing."""
    return get_real_cif_path("2y2a", pdb_mmcif_path)
