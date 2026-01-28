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

        # 5i55.cif: entities 1 (polymer), 2 (non-polymer), 3 (non-polymer), 4 (water)
        # Keep entities 1, 2, 3 (remove water entity 4)
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3"]),
            chain_ids=frozenset(["A", "B", "C"]),
        )

        sync_entity_categories(block, info)

        # Check _entity table
        entity_table = block.find("_entity.", ["id", "type"])
        entity_ids = [row[0] for row in entity_table]

        assert "4" not in entity_ids  # water removed
        assert "1" in entity_ids  # polymer still present
        assert "2" in entity_ids  # non-polymer still present
        assert "3" in entity_ids  # non-polymer still present

    def test_removes_entity_poly_row(self) -> None:
        """Test that entity_poly row is removed for removed polymer."""
        # Use synthetic CIF with _entity_poly in loop format
        # (real files may use pair format for single-row categories)
        cif_content = """\
data_TEST
#
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_strand_id
1 'polypeptide(L)' A
2 'polypeptide(L)' B
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Keep only entity 2 (remove polymer entity 1)
        info = StructureInfo(
            entity_ids=frozenset(["2"]),
            chain_ids=frozenset(["B"]),
        )

        sync_entity_categories(block, info)

        # Check _entity_poly table - entity 1 should be removed
        entity_poly = block.find("_entity_poly.", ["entity_id"])
        entity_ids = [row[0] for row in entity_poly]

        assert "1" not in entity_ids
        assert "2" in entity_ids

    def test_no_changes_when_all_kept(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that nothing changes when all entities are kept."""
        block = sample_cif_document[0]

        # Get original counts
        original_entity_count = len(list(block.find("_entity.", ["id"])))

        # 5i55.cif has 4 entities and 4 chains
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3", "4"]),
            chain_ids=frozenset(["A", "B", "C", "D"]),
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


class TestSyncEntityPolySeq:
    """Tests for _entity_poly_seq synchronization."""

    def test_removes_entity_poly_seq_rows(self) -> None:
        """Test that _entity_poly_seq rows are removed for removed entities."""
        cif_content = """\
data_TEST
#
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 ALA n
1 2 GLY n
1 3 SER n
2 1 VAL n
2 2 LEU n
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Keep only entity 1
        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        sync_entity_categories(block, info)

        # Check _entity_poly_seq - entity 2 rows should be removed
        seq_table = block.find("_entity_poly_seq.", ["entity_id", "num"])
        entity_ids = [row[0] for row in seq_table]

        assert "2" not in entity_ids
        assert entity_ids.count("1") == 3

    def test_keeps_all_entity_poly_seq_when_all_kept(self) -> None:
        """Test that all _entity_poly_seq rows remain when all entities kept."""
        cif_content = """\
data_TEST
#
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 ALA n
1 2 GLY n
2 1 VAL n
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Keep both entities
        info = StructureInfo(
            entity_ids=frozenset(["1", "2"]),
            chain_ids=frozenset(["A", "B"]),
        )

        original_count = len(list(block.find("_entity_poly_seq.", ["entity_id"])))
        sync_entity_categories(block, info)
        new_count = len(list(block.find("_entity_poly_seq.", ["entity_id"])))

        assert new_count == original_count
        assert new_count == 3

    def test_handles_partial_entity_sync(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that sync works with partial entity set."""
        block = sample_cif_document[0]

        # Keep only entity 1
        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        # Should not raise
        sync_entity_categories(block, info)

        # Verify only entity 1 remains
        entity_ids = [row[0] for row in block.find("_entity.", ["id"])]
        assert entity_ids == ["1"]
