"""Tests for restore_categories function."""

from pathlib import Path

import gemmi
import pytest

from mmcif_restore import restore_categories


class TestRestoreCategories:
    """Tests for restore_categories function."""

    def test_restores_entity(self, sample_cif_file: Path, tmp_path: Path) -> None:
        """Test restoring _entity category."""
        # Create minimal CIF (simulate external editor output)
        structure = gemmi.read_structure(str(sample_cif_file))
        structure.remove_waters()

        groups = gemmi.MmcifOutputGroups(False)
        groups.atoms = True
        groups.entry = True
        minimal_doc = structure.make_mmcif_document(groups)
        minimal_path = tmp_path / "minimal.cif"
        minimal_doc.write_file(str(minimal_path))

        # Restore _entity from reference
        doc = restore_categories(
            minimal_path,
            sample_cif_file,
            categories=["_entity."],
        )

        # Check _entity was restored and synced (no water)
        block = doc[0]
        entity_types = [row[1] for row in block.find("_entity.", ["id", "type"])]

        assert "polymer" in entity_types
        assert "non-polymer" in entity_types
        assert "water" not in entity_types

    def test_restores_struct_asym(self, sample_cif_file: Path, tmp_path: Path) -> None:
        """Test restoring _struct_asym category."""
        # Create minimal CIF
        structure = gemmi.read_structure(str(sample_cif_file))
        structure.remove_waters()

        groups = gemmi.MmcifOutputGroups(False)
        groups.atoms = True
        groups.entry = True
        minimal_doc = structure.make_mmcif_document(groups)
        minimal_path = tmp_path / "minimal.cif"
        minimal_doc.write_file(str(minimal_path))

        # Restore _struct_asym from reference
        doc = restore_categories(
            minimal_path,
            sample_cif_file,
            categories=["_struct_asym."],
        )

        # Check _struct_asym was restored and synced (no water chain)
        block = doc[0]
        chain_ids = [row[0] for row in block.find("_struct_asym.", ["id"])]

        assert "A" in chain_ids
        assert "B" in chain_ids
        assert "C" not in chain_ids  # Water chain removed

    def test_restores_multiple_categories(
        self, sample_cif_file: Path, tmp_path: Path
    ) -> None:
        """Test restoring multiple categories at once."""
        # Create minimal CIF
        structure = gemmi.read_structure(str(sample_cif_file))
        structure.remove_waters()

        groups = gemmi.MmcifOutputGroups(False)
        groups.atoms = True
        groups.entry = True
        minimal_doc = structure.make_mmcif_document(groups)
        minimal_path = tmp_path / "minimal.cif"
        minimal_doc.write_file(str(minimal_path))

        # Restore multiple categories
        doc = restore_categories(
            minimal_path,
            sample_cif_file,
            categories=["_entity.", "_struct_asym."],
        )

        # Check both categories exist
        block = doc[0]
        entity_count = len(list(block.find("_entity.", ["id"])))
        asym_count = len(list(block.find("_struct_asym.", ["id"])))

        assert entity_count == 2  # polymer and non-polymer (no water)
        assert asym_count == 2  # A and B (no water chain C)

    def test_normalizes_category_prefix(
        self, sample_cif_file: Path, tmp_path: Path
    ) -> None:
        """Test that category prefixes are normalized (with or without dot)."""
        # Create minimal CIF
        structure = gemmi.read_structure(str(sample_cif_file))
        structure.remove_waters()

        groups = gemmi.MmcifOutputGroups(False)
        groups.atoms = True
        groups.entry = True
        minimal_doc = structure.make_mmcif_document(groups)
        minimal_path = tmp_path / "minimal.cif"
        minimal_doc.write_file(str(minimal_path))

        # Restore with prefix without dot
        doc = restore_categories(
            minimal_path,
            sample_cif_file,
            categories=["_entity"],  # No trailing dot
        )

        # Should still work
        block = doc[0]
        entity_count = len(list(block.find("_entity.", ["id"])))
        assert entity_count == 2
