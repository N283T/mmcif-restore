"""Tests for entity synchronization module."""

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync.entity import sync_entity_categories


class TestSyncEntityCategories:
    """Tests for sync_entity_categories function."""

    def test_removes_water_entity(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that water entity is removed from _entity."""
        block = sample_cif_document[0]

        # Keep only entities 1 and 2 (remove water entity 3)
        info = StructureInfo(
            entity_ids=frozenset(["1", "2"]),
            chain_ids=frozenset(["A", "B"]),
        )

        sync_entity_categories(block, info)

        # Check _entity table
        entity_table = block.find("_entity.", ["id", "type"])
        entity_ids = [row[0] for row in entity_table]

        assert "3" not in entity_ids
        assert "1" in entity_ids  # polymer still present
        assert "2" in entity_ids  # non-polymer still present

    def test_removes_entity_poly_row(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that entity_poly row is removed for removed polymer."""
        block = sample_cif_document[0]

        # Keep only entities 2 and 3 (remove polymer entity 1)
        info = StructureInfo(
            entity_ids=frozenset(["2", "3"]),
            chain_ids=frozenset(["B", "C"]),
        )

        sync_entity_categories(block, info)

        # Check _entity_poly table - should be empty after removing the only polymer
        entity_poly = block.find("_entity_poly.", ["entity_id"])
        entity_ids = [row[0] for row in entity_poly]

        assert "1" not in entity_ids

    def test_no_changes_when_all_kept(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that nothing changes when all entities are kept."""
        block = sample_cif_document[0]

        # Get original counts
        original_entity_count = len(list(block.find("_entity.", ["id"])))

        # Keep all entities
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3"]),
            chain_ids=frozenset(["A", "B", "C"]),
        )

        sync_entity_categories(block, info)

        # Counts should be unchanged
        new_entity_count = len(list(block.find("_entity.", ["id"])))
        assert new_entity_count == original_entity_count

    def test_removes_multiple_entities(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test removing multiple entities at once."""
        block = sample_cif_document[0]

        # Keep only entity 1 (remove both non-polymer and water)
        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        sync_entity_categories(block, info)

        # Only polymer should remain
        entity_table = block.find("_entity.", ["id", "type"])
        entity_ids = [row[0] for row in entity_table]

        assert entity_ids == ["1"]
