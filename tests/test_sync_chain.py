"""Tests for chain synchronization module."""

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync.chain import sync_chain_categories


class TestSyncChainCategories:
    """Tests for sync_chain_categories function."""

    def test_removes_water_chain(self, sample_cif_document: gemmi.cif.Document) -> None:
        """Test that water chain is removed from _struct_asym."""
        block = sample_cif_document[0]

        # Keep only chains A and B (remove water chain C)
        info = StructureInfo(
            entity_ids=frozenset(["1", "2"]),
            chain_ids=frozenset(["A", "B"]),
        )

        sync_chain_categories(block, info)

        # Check _struct_asym table
        struct_asym = block.find("_struct_asym.", ["id", "entity_id"])
        chain_ids = [row[0] for row in struct_asym]

        assert "C" not in chain_ids
        assert "A" in chain_ids  # polymer chain still present
        assert "B" in chain_ids  # non-polymer chain still present

    def test_removes_multiple_chains(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test removing multiple chains at once."""
        block = sample_cif_document[0]

        # Keep only chain A (remove both non-polymer and water chains)
        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        sync_chain_categories(block, info)

        # Only polymer chain should remain
        struct_asym = block.find("_struct_asym.", ["id", "entity_id"])
        chain_ids = [row[0] for row in struct_asym]

        assert chain_ids == ["A"]

    def test_no_changes_when_all_kept(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that nothing changes when all chains are kept."""
        block = sample_cif_document[0]

        # Get original count
        original_count = len(list(block.find("_struct_asym.", ["id"])))

        # Keep all chains
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3"]),
            chain_ids=frozenset(["A", "B", "C"]),
        )

        sync_chain_categories(block, info)

        # Count should be unchanged
        new_count = len(list(block.find("_struct_asym.", ["id"])))
        assert new_count == original_count
