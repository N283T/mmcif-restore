"""Test fixtures for mmcif-restore tests."""

from pathlib import Path

import gemmi
import pytest

# Path to test fixtures
FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def sample_cif_file() -> Path:
    """Path to sample mmCIF file (5i55.cif from Gemmi tests).

    Structure contains:
    - Entity 1: polymer (chain A)
    - Entity 2: non-polymer (chain B)
    - Entity 3: non-polymer (chain C)
    - Entity 4: water (chain D)
    - 1 struct_conn record
    """
    return FIXTURES_DIR / "5i55.cif"


@pytest.fixture
def sample_structure(sample_cif_file: Path) -> gemmi.Structure:
    """Load sample structure for testing."""
    return gemmi.read_structure(str(sample_cif_file))


@pytest.fixture
def sample_cif_document(sample_cif_file: Path) -> gemmi.cif.Document:
    """Load sample CIF document for testing."""
    return gemmi.cif.read(str(sample_cif_file))
