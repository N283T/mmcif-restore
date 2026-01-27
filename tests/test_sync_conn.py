"""Tests for connection synchronization module."""

import gemmi
import pytest

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync.conn import sync_conn_categories


@pytest.fixture
def cif_with_struct_conn() -> gemmi.cif.Document:
    """Create a CIF document with _struct_conn category."""
    cif_content = """\
data_TEST
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
disulf1 disulf A CYS 10 A CYS 20
metalc1 metalc A HIS 50 B ZN .
hydrog1 hydrog A ALA 5 C HOH .
covale1 covale A SER 30 A SER 35
#
"""
    return gemmi.cif.read_string(cif_content)


class TestSyncConnCategories:
    """Tests for sync_conn_categories function."""

    def test_removes_connection_when_partner_removed(
        self, cif_with_struct_conn: gemmi.cif.Document
    ) -> None:
        """Test that connection is removed when one partner's chain is removed."""
        block = cif_with_struct_conn[0]

        # Keep only chains A and C (remove chain B)
        info = StructureInfo(
            entity_ids=frozenset(["1", "3"]),
            chain_ids=frozenset(["A", "C"]),
        )

        sync_conn_categories(block, info)

        # Check _struct_conn table
        conn_table = block.find("_struct_conn.", ["id", "ptnr1_label_asym_id"])
        conn_ids = [row[0] for row in conn_table]

        # metalc1 should be removed (involves chain B)
        assert "metalc1" not in conn_ids
        # disulf1 and covale1 should remain (both partners in chain A)
        assert "disulf1" in conn_ids
        assert "covale1" in conn_ids

    def test_removes_connection_to_water(
        self, cif_with_struct_conn: gemmi.cif.Document
    ) -> None:
        """Test that connection to water is removed when water removed."""
        block = cif_with_struct_conn[0]

        # Keep only chains A and B (remove water chain C)
        info = StructureInfo(
            entity_ids=frozenset(["1", "2"]),
            chain_ids=frozenset(["A", "B"]),
        )

        sync_conn_categories(block, info)

        # Check _struct_conn table
        conn_table = block.find("_struct_conn.", ["id"])
        conn_ids = [row[0] for row in conn_table]

        # hydrog1 should be removed (involves chain C/water)
        assert "hydrog1" not in conn_ids
        # Others should remain
        assert "disulf1" in conn_ids
        assert "metalc1" in conn_ids
        assert "covale1" in conn_ids

    def test_removes_multiple_connections(
        self, cif_with_struct_conn: gemmi.cif.Document
    ) -> None:
        """Test removing multiple connections at once."""
        block = cif_with_struct_conn[0]

        # Keep only chain A (remove both B and C)
        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        sync_conn_categories(block, info)

        # Check _struct_conn table
        conn_table = block.find("_struct_conn.", ["id"])
        conn_ids = [row[0] for row in conn_table]

        # metalc1 and hydrog1 should be removed
        assert "metalc1" not in conn_ids
        assert "hydrog1" not in conn_ids
        # intra-chain A connections should remain
        assert "disulf1" in conn_ids
        assert "covale1" in conn_ids

    def test_no_changes_when_all_kept(
        self, cif_with_struct_conn: gemmi.cif.Document
    ) -> None:
        """Test that nothing changes when all chains are kept."""
        block = cif_with_struct_conn[0]

        # Get original count
        original_count = len(list(block.find("_struct_conn.", ["id"])))

        # Keep all chains
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3"]),
            chain_ids=frozenset(["A", "B", "C"]),
        )

        sync_conn_categories(block, info)

        # Count should be unchanged
        new_count = len(list(block.find("_struct_conn.", ["id"])))
        assert new_count == original_count
