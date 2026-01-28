"""Tests for modified residue synchronization module."""

import gemmi
import pytest

from mmcif_restore.sync.modres import sync_modres_categories


def _make_atom(name: str) -> gemmi.Atom:
    """Create an atom with the given name."""
    atom = gemmi.Atom()
    atom.name = name
    return atom


@pytest.fixture
def cif_with_modres() -> gemmi.cif.Document:
    """Create a CIF document with _pdbx_struct_mod_residue category."""
    cif_content = """\
data_TEST
#
loop_
_pdbx_struct_mod_residue.id
_pdbx_struct_mod_residue.auth_asym_id
_pdbx_struct_mod_residue.auth_seq_id
_pdbx_struct_mod_residue.PDB_ins_code
_pdbx_struct_mod_residue.auth_comp_id
_pdbx_struct_mod_residue.label_comp_id
_pdbx_struct_mod_residue.parent_comp_id
_pdbx_struct_mod_residue.details
1 A 10 ? MSE MSE MET 'SELENOMETHIONINE'
2 A 25 ? MSE MSE MET 'SELENOMETHIONINE'
3 B 5 ? PTR PTR TYR 'PHOSPHOTYROSINE'
#
"""
    return gemmi.cif.read_string(cif_content)


@pytest.fixture
def structure_with_all_residues() -> gemmi.Structure:
    """Create a structure with all residues referenced in _pdbx_struct_mod_residue."""
    st = gemmi.Structure()
    model = gemmi.Model(1)

    # Chain A with MSE 10 and MSE 25
    chain_a = gemmi.Chain("A")

    res_mse10 = gemmi.Residue()
    res_mse10.name = "MSE"
    res_mse10.seqid = gemmi.SeqId("10")
    res_mse10.add_atom(_make_atom("CA"))
    chain_a.add_residue(res_mse10)

    res_mse25 = gemmi.Residue()
    res_mse25.name = "MSE"
    res_mse25.seqid = gemmi.SeqId("25")
    res_mse25.add_atom(_make_atom("CA"))
    chain_a.add_residue(res_mse25)

    model.add_chain(chain_a)

    # Chain B with PTR 5
    chain_b = gemmi.Chain("B")

    res_ptr5 = gemmi.Residue()
    res_ptr5.name = "PTR"
    res_ptr5.seqid = gemmi.SeqId("5")
    res_ptr5.add_atom(_make_atom("CA"))
    chain_b.add_residue(res_ptr5)

    model.add_chain(chain_b)

    st.add_model(model)
    return st


class TestSyncModresCategories:
    """Tests for sync_modres_categories function."""

    def test_keeps_all_modres_when_all_residues_present(
        self,
        cif_with_modres: gemmi.cif.Document,
        structure_with_all_residues: gemmi.Structure,
    ) -> None:
        """Test that all MODRES records are kept when all residues exist."""
        block = cif_with_modres[0]
        original_count = len(list(block.find("_pdbx_struct_mod_residue.", ["id"])))

        sync_modres_categories(block, structure_with_all_residues)

        new_count = len(list(block.find("_pdbx_struct_mod_residue.", ["id"])))
        assert new_count == original_count
        assert new_count == 3

    def test_removes_modres_when_chain_removed(
        self,
        cif_with_modres: gemmi.cif.Document,
        structure_with_all_residues: gemmi.Structure,
    ) -> None:
        """Test that MODRES is removed when its chain is removed."""
        block = cif_with_modres[0]

        # Remove chain B
        model = structure_with_all_residues[0]
        model.remove_chain("B")

        sync_modres_categories(block, structure_with_all_residues)

        modres_ids = [row[0] for row in block.find("_pdbx_struct_mod_residue.", ["id"])]

        # PTR (id=3) should be removed
        assert "3" not in modres_ids
        # MSE entries should remain
        assert "1" in modres_ids
        assert "2" in modres_ids

    def test_removes_modres_when_specific_residue_removed(
        self,
        cif_with_modres: gemmi.cif.Document,
        structure_with_all_residues: gemmi.Structure,
    ) -> None:
        """Test that MODRES is removed when the specific residue is removed."""
        block = cif_with_modres[0]

        # Remove MSE 10 from chain A (keep MSE 25)
        chain_a = structure_with_all_residues[0]["A"]
        for i, res in enumerate(chain_a):
            if res.name == "MSE" and res.seqid.num == 10:
                del chain_a[i]
                break

        sync_modres_categories(block, structure_with_all_residues)

        modres_ids = [row[0] for row in block.find("_pdbx_struct_mod_residue.", ["id"])]

        # MSE 10 (id=1) should be removed
        assert "1" not in modres_ids
        # Others should remain
        assert "2" in modres_ids
        assert "3" in modres_ids

    def test_removes_multiple_modres(
        self,
        cif_with_modres: gemmi.cif.Document,
        structure_with_all_residues: gemmi.Structure,
    ) -> None:
        """Test removing multiple MODRES records at once."""
        block = cif_with_modres[0]

        # Remove chain A entirely
        model = structure_with_all_residues[0]
        model.remove_chain("A")

        sync_modres_categories(block, structure_with_all_residues)

        modres_ids = [row[0] for row in block.find("_pdbx_struct_mod_residue.", ["id"])]

        # Both MSE entries should be removed
        assert "1" not in modres_ids
        assert "2" not in modres_ids
        # PTR should remain
        assert "3" in modres_ids


class TestSyncModresWithInsertionCode:
    """Tests for handling insertion codes in MODRES."""

    def test_handles_insertion_codes(self) -> None:
        """Test that insertion codes are properly matched."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_struct_mod_residue.id
_pdbx_struct_mod_residue.auth_asym_id
_pdbx_struct_mod_residue.auth_seq_id
_pdbx_struct_mod_residue.PDB_ins_code
_pdbx_struct_mod_residue.auth_comp_id
_pdbx_struct_mod_residue.parent_comp_id
1 A 10 ? MSE MET
2 A 10 A MSE MET
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Create structure with residue 10 (no insertion) and 10A
        st = gemmi.Structure()
        model = gemmi.Model(1)
        chain = gemmi.Chain("A")

        # Residue 10 without insertion code
        res1 = gemmi.Residue()
        res1.name = "MSE"
        res1.seqid = gemmi.SeqId("10")
        res1.add_atom(_make_atom("CA"))
        chain.add_residue(res1)

        # Residue 10A with insertion code
        res2 = gemmi.Residue()
        res2.name = "MSE"
        res2.seqid = gemmi.SeqId(10, "A")
        res2.add_atom(_make_atom("CA"))
        chain.add_residue(res2)

        model.add_chain(chain)
        st.add_model(model)

        sync_modres_categories(block, st)

        modres_ids = [row[0] for row in block.find("_pdbx_struct_mod_residue.", ["id"])]

        # Both should remain
        assert "1" in modres_ids
        assert "2" in modres_ids

    def test_removes_modres_with_missing_insertion_code(self) -> None:
        """Test that MODRES is removed if insertion code doesn't match."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_struct_mod_residue.id
_pdbx_struct_mod_residue.auth_asym_id
_pdbx_struct_mod_residue.auth_seq_id
_pdbx_struct_mod_residue.PDB_ins_code
_pdbx_struct_mod_residue.auth_comp_id
_pdbx_struct_mod_residue.parent_comp_id
1 A 10 A MSE MET
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Create structure with only residue 10 (no insertion code)
        st = gemmi.Structure()
        model = gemmi.Model(1)
        chain = gemmi.Chain("A")

        res1 = gemmi.Residue()
        res1.name = "MSE"
        res1.seqid = gemmi.SeqId("10")  # No insertion code
        res1.add_atom(_make_atom("CA"))
        chain.add_residue(res1)

        model.add_chain(chain)
        st.add_model(model)

        sync_modres_categories(block, st)

        modres_ids = [row[0] for row in block.find("_pdbx_struct_mod_residue.", ["id"])]

        # id=1 should be removed (10A doesn't exist, only 10)
        assert "1" not in modres_ids


class TestSyncModresEdgeCases:
    """Tests for edge cases in sync_modres_categories."""

    def test_handles_empty_structure(self) -> None:
        """Test that sync handles empty structure gracefully."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_struct_mod_residue.id
_pdbx_struct_mod_residue.auth_asym_id
_pdbx_struct_mod_residue.auth_seq_id
_pdbx_struct_mod_residue.auth_comp_id
1 A 10 MSE
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        empty_structure = gemmi.Structure()
        sync_modres_categories(block, empty_structure)

        # Should return early, loop unchanged
        modres_ids = [row[0] for row in block.find("_pdbx_struct_mod_residue.", ["id"])]
        assert "1" in modres_ids

    def test_handles_no_modres_loop(self) -> None:
        """Test that sync handles CIF without _pdbx_struct_mod_residue loop."""
        cif_content = """\
data_TEST
#
loop_
_entity.id
_entity.type
1 polymer
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        st = gemmi.Structure()
        model = gemmi.Model(1)
        chain = gemmi.Chain("A")
        res = gemmi.Residue()
        res.name = "MSE"
        res.seqid = gemmi.SeqId("1")
        res.add_atom(_make_atom("CA"))
        chain.add_residue(res)
        model.add_chain(chain)
        st.add_model(model)

        # Should not raise
        sync_modres_categories(block, st)

    def test_handles_modres_without_insertion_code_column(self) -> None:
        """Test sync works when PDB_ins_code column is absent."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_struct_mod_residue.id
_pdbx_struct_mod_residue.auth_asym_id
_pdbx_struct_mod_residue.auth_seq_id
_pdbx_struct_mod_residue.auth_comp_id
_pdbx_struct_mod_residue.parent_comp_id
1 A 10 MSE MET
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        st = gemmi.Structure()
        model = gemmi.Model(1)
        chain = gemmi.Chain("A")

        res = gemmi.Residue()
        res.name = "MSE"
        res.seqid = gemmi.SeqId("10")
        res.add_atom(_make_atom("CA"))
        chain.add_residue(res)

        model.add_chain(chain)
        st.add_model(model)

        sync_modres_categories(block, st)

        # Should remain (residue exists)
        modres_ids = [row[0] for row in block.find("_pdbx_struct_mod_residue.", ["id"])]
        assert "1" in modres_ids
