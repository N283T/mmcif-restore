"""Tests for StructureInfo class."""

from pathlib import Path

import gemmi

from mmcif_restore.structure_info import StructureInfo


class TestStructureInfoFromStructure:
    """Tests for from_structure class method."""

    def test_handles_empty_structure(self) -> None:
        """Test behavior with empty structure."""
        structure = gemmi.Structure()
        info = StructureInfo.from_structure(structure)
        assert info.entity_ids == frozenset()
        assert info.chain_ids == frozenset()

    def test_extracts_entity_ids(self, sample_structure: gemmi.Structure) -> None:
        """Test that entity IDs are extracted from structure."""
        info = StructureInfo.from_structure(sample_structure)

        # 5i55.cif has entities 1, 2, 3, 4 (polymer, 2x non-polymer, water)
        assert "1" in info.entity_ids
        assert "2" in info.entity_ids
        assert "3" in info.entity_ids
        assert "4" in info.entity_ids

    def test_extracts_chain_ids(self, sample_structure: gemmi.Structure) -> None:
        """Test that chain IDs are extracted from structure."""
        info = StructureInfo.from_structure(sample_structure)

        # 5i55.cif has chains A, B, C, D
        assert "A" in info.chain_ids
        assert "B" in info.chain_ids
        assert "C" in info.chain_ids
        assert "D" in info.chain_ids


class TestStructureInfoFromStructureWithReference:
    """Tests for from_structure_with_reference class method."""

    def test_extracts_entity_ids_from_reference(
        self,
        sample_structure: gemmi.Structure,
        sample_cif_document: gemmi.cif.Document,
    ) -> None:
        """Test entity ID extraction using reference CIF."""
        # Remove water from structure
        sample_structure.remove_waters()

        # Get info with reference
        info = StructureInfo.from_structure_with_reference(
            sample_structure,
            sample_cif_document[0],
        )

        # 5i55.cif: after removing water, should have entities 1, 2, 3 (not 4/water)
        assert "1" in info.entity_ids
        assert "2" in info.entity_ids
        assert "3" in info.entity_ids
        assert "4" not in info.entity_ids

    def test_extracts_chain_ids_without_water(
        self,
        sample_structure: gemmi.Structure,
        sample_cif_document: gemmi.cif.Document,
    ) -> None:
        """Test chain ID extraction after water removal."""
        # Remove water from structure
        sample_structure.remove_waters()

        # Get info with reference
        info = StructureInfo.from_structure_with_reference(
            sample_structure,
            sample_cif_document[0],
        )

        # 5i55.cif: after removing water, should have chains A, B, C (not D/water)
        assert "A" in info.chain_ids
        assert "B" in info.chain_ids
        assert "C" in info.chain_ids
        assert "D" not in info.chain_ids

    def test_handles_minimal_cif_without_entities(self, tmp_path: Path) -> None:
        """Test with a minimal CIF that has no entity info."""
        # Create minimal CIF content (only _atom_site)
        minimal_content = """\
data_TEST
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
#
"""
        # Reference CIF with full info
        ref_content = """\
data_TEST
#
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
#
"""
        # Write minimal CIF to temp file and read structure
        minimal_path = tmp_path / "minimal.cif"
        minimal_path.write_text(minimal_content)

        ref_doc = gemmi.cif.read_string(ref_content)

        # Read structure from minimal file
        structure = gemmi.read_structure(str(minimal_path))

        # Get info with reference
        info = StructureInfo.from_structure_with_reference(structure, ref_doc[0])

        # Should have chain A and entity 1
        assert "A" in info.chain_ids
        assert "1" in info.entity_ids
