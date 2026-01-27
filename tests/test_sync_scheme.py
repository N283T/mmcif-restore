"""Tests for scheme synchronization module."""

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync.scheme import sync_scheme_categories


class TestSyncSchemeCategories:
    """Tests for sync_scheme_categories function."""

    def test_removes_nonpoly_scheme_rows(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that _pdbx_nonpoly_scheme rows are removed for removed chains."""
        block = sample_cif_document[0]

        # Keep only chains A and B (remove water chain C)
        info = StructureInfo(
            entity_ids=frozenset(["1", "2"]),
            chain_ids=frozenset(["A", "B"]),
        )

        sync_scheme_categories(block, info)

        # Check _pdbx_nonpoly_scheme table
        nonpoly_scheme = block.find("_pdbx_nonpoly_scheme.", ["asym_id", "mon_id"])
        asym_ids = [row[0] for row in nonpoly_scheme]

        assert "C" not in asym_ids
        assert "B" in asym_ids  # non-polymer still present

    def test_removes_multiple_scheme_rows(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test removing multiple scheme rows."""
        block = sample_cif_document[0]

        # Keep only chain A (remove both non-polymer and water)
        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        sync_scheme_categories(block, info)

        # _pdbx_nonpoly_scheme should be empty
        nonpoly_scheme = block.find("_pdbx_nonpoly_scheme.", ["asym_id"])
        asym_ids = [row[0] for row in nonpoly_scheme]

        assert len(asym_ids) == 0

    def test_no_changes_when_all_kept(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that nothing changes when all chains are kept."""
        block = sample_cif_document[0]

        # Get original count
        original_count = len(list(block.find("_pdbx_nonpoly_scheme.", ["asym_id"])))

        # Keep all chains
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3"]),
            chain_ids=frozenset(["A", "B", "C"]),
        )

        sync_scheme_categories(block, info)

        # Count should be unchanged
        new_count = len(list(block.find("_pdbx_nonpoly_scheme.", ["asym_id"])))
        assert new_count == original_count
