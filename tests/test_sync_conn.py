"""Tests for connection synchronization module."""

import gemmi
import pytest

from mmcif_restore.sync.conn import sync_conn_categories


def _make_atom(name: str) -> gemmi.Atom:
    """Create an atom with the given name."""
    atom = gemmi.Atom()
    atom.name = name
    return atom


@pytest.fixture
def cif_with_struct_conn() -> gemmi.cif.Document:
    """Create a CIF document with _struct_conn category using auth columns."""
    cif_content = """\
data_TEST
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.pdbx_ptnr2_PDB_ins_code
disulf1 disulf A 10 CYS SG A 20 CYS SG ? ?
metalc1 metalc A 50 HIS NE2 B 1 ZN ZN ? ?
hydrog1 hydrog A 5 ALA O C 1 HOH O ? ?
covale1 covale A 30 SER OG A 35 SER CB ? ?
#
"""
    return gemmi.cif.read_string(cif_content)


@pytest.fixture
def structure_with_all_atoms() -> gemmi.Structure:
    """Create a structure with all atoms referenced in _struct_conn."""
    st = gemmi.Structure()
    model = gemmi.Model("1")

    # Chain A with CYS 10, CYS 20, HIS 50, ALA 5, SER 30, SER 35
    chain_a = gemmi.Chain("A")

    res_ala5 = gemmi.Residue()
    res_ala5.name = "ALA"
    res_ala5.seqid = gemmi.SeqId("5")
    res_ala5.add_atom(_make_atom("O"))
    chain_a.add_residue(res_ala5)

    res_cys10 = gemmi.Residue()
    res_cys10.name = "CYS"
    res_cys10.seqid = gemmi.SeqId("10")
    res_cys10.add_atom(_make_atom("SG"))
    chain_a.add_residue(res_cys10)

    res_cys20 = gemmi.Residue()
    res_cys20.name = "CYS"
    res_cys20.seqid = gemmi.SeqId("20")
    res_cys20.add_atom(_make_atom("SG"))
    chain_a.add_residue(res_cys20)

    res_ser30 = gemmi.Residue()
    res_ser30.name = "SER"
    res_ser30.seqid = gemmi.SeqId("30")
    res_ser30.add_atom(_make_atom("OG"))
    chain_a.add_residue(res_ser30)

    res_ser35 = gemmi.Residue()
    res_ser35.name = "SER"
    res_ser35.seqid = gemmi.SeqId("35")
    res_ser35.add_atom(_make_atom("CB"))
    chain_a.add_residue(res_ser35)

    res_his50 = gemmi.Residue()
    res_his50.name = "HIS"
    res_his50.seqid = gemmi.SeqId("50")
    res_his50.add_atom(_make_atom("NE2"))
    chain_a.add_residue(res_his50)

    model.add_chain(chain_a)

    # Chain B with ZN 1
    chain_b = gemmi.Chain("B")
    res_zn = gemmi.Residue()
    res_zn.name = "ZN"
    res_zn.seqid = gemmi.SeqId("1")
    res_zn.add_atom(_make_atom("ZN"))
    chain_b.add_residue(res_zn)
    model.add_chain(chain_b)

    # Chain C with HOH 1
    chain_c = gemmi.Chain("C")
    res_hoh = gemmi.Residue()
    res_hoh.name = "HOH"
    res_hoh.seqid = gemmi.SeqId("1")
    res_hoh.add_atom(_make_atom("O"))
    chain_c.add_residue(res_hoh)
    model.add_chain(chain_c)

    st.add_model(model)
    return st


class TestSyncConnCategories:
    """Tests for sync_conn_categories function."""

    def test_keeps_all_connections_when_all_atoms_present(
        self,
        cif_with_struct_conn: gemmi.cif.Document,
        structure_with_all_atoms: gemmi.Structure,
    ) -> None:
        """Test that all connections are kept when all atoms exist."""
        block = cif_with_struct_conn[0]
        original_count = len(list(block.find("_struct_conn.", ["id"])))

        sync_conn_categories(block, structure_with_all_atoms)

        new_count = len(list(block.find("_struct_conn.", ["id"])))
        assert new_count == original_count

    def test_removes_connection_when_partner_chain_removed(
        self,
        cif_with_struct_conn: gemmi.cif.Document,
        structure_with_all_atoms: gemmi.Structure,
    ) -> None:
        """Test that connection is removed when one partner's chain is removed."""
        block = cif_with_struct_conn[0]

        # Remove chain B (ZN)
        model = structure_with_all_atoms[0]
        model.remove_chain("B")

        sync_conn_categories(block, structure_with_all_atoms)

        conn_ids = [row[0] for row in block.find("_struct_conn.", ["id"])]

        # metalc1 should be removed (involves chain B)
        assert "metalc1" not in conn_ids
        # Others should remain
        assert "disulf1" in conn_ids
        assert "hydrog1" in conn_ids
        assert "covale1" in conn_ids

    def test_removes_connection_when_specific_residue_removed(
        self,
        cif_with_struct_conn: gemmi.cif.Document,
        structure_with_all_atoms: gemmi.Structure,
    ) -> None:
        """Test that connection is removed when specific residue is removed."""
        block = cif_with_struct_conn[0]

        # Remove CYS 10 from chain A (keep CYS 20)
        chain_a = structure_with_all_atoms[0]["A"]
        for i, res in enumerate(chain_a):
            if res.name == "CYS" and res.seqid.num == 10:
                del chain_a[i]
                break

        sync_conn_categories(block, structure_with_all_atoms)

        conn_ids = [row[0] for row in block.find("_struct_conn.", ["id"])]

        # disulf1 should be removed (CYS 10 removed)
        assert "disulf1" not in conn_ids
        # Others should remain
        assert "metalc1" in conn_ids
        assert "hydrog1" in conn_ids
        assert "covale1" in conn_ids

    def test_removes_connection_to_water(
        self,
        cif_with_struct_conn: gemmi.cif.Document,
        structure_with_all_atoms: gemmi.Structure,
    ) -> None:
        """Test that connection to water is removed when water removed."""
        block = cif_with_struct_conn[0]

        # Remove water chain C
        model = structure_with_all_atoms[0]
        model.remove_chain("C")

        sync_conn_categories(block, structure_with_all_atoms)

        conn_ids = [row[0] for row in block.find("_struct_conn.", ["id"])]

        # hydrog1 should be removed (involves chain C/water)
        assert "hydrog1" not in conn_ids
        # Others should remain
        assert "disulf1" in conn_ids
        assert "metalc1" in conn_ids
        assert "covale1" in conn_ids

    def test_removes_multiple_connections(
        self,
        cif_with_struct_conn: gemmi.cif.Document,
        structure_with_all_atoms: gemmi.Structure,
    ) -> None:
        """Test removing multiple connections at once."""
        block = cif_with_struct_conn[0]

        # Remove chains B and C
        model = structure_with_all_atoms[0]
        model.remove_chain("B")
        model.remove_chain("C")

        sync_conn_categories(block, structure_with_all_atoms)

        conn_ids = [row[0] for row in block.find("_struct_conn.", ["id"])]

        # metalc1 and hydrog1 should be removed
        assert "metalc1" not in conn_ids
        assert "hydrog1" not in conn_ids
        # intra-chain A connections should remain
        assert "disulf1" in conn_ids
        assert "covale1" in conn_ids


class TestSyncConnWithInsertionCode:
    """Tests for handling insertion codes."""

    def test_handles_insertion_codes(self) -> None:
        """Test that insertion codes are properly matched."""
        cif_content = """\
data_TEST
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.pdbx_ptnr2_PDB_ins_code
conn1 disulf A 10 CYS SG A 10 CYS SG ? A
conn2 disulf A 10 CYS SG A 10 CYS SG A ?
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Create structure with residue 10 (no insertion) and 10A
        st = gemmi.Structure()
        model = gemmi.Model("1")
        chain = gemmi.Chain("A")

        # Residue 10 without insertion code
        res1 = gemmi.Residue()
        res1.name = "CYS"
        res1.seqid = gemmi.SeqId("10")
        res1.add_atom(_make_atom("SG"))
        chain.add_residue(res1)

        # Residue 10A with insertion code
        res2 = gemmi.Residue()
        res2.name = "CYS"
        res2.seqid = gemmi.SeqId("10", "A")
        res2.add_atom(_make_atom("SG"))
        chain.add_residue(res2)

        model.add_chain(chain)
        st.add_model(model)

        sync_conn_categories(block, st)

        conn_ids = [row[0] for row in block.find("_struct_conn.", ["id"])]

        # Both connections should remain (10 and 10A both exist)
        assert "conn1" in conn_ids
        assert "conn2" in conn_ids

    def test_removes_connection_with_missing_insertion_code(self) -> None:
        """Test that connection is removed if insertion code doesn't match."""
        cif_content = """\
data_TEST
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.pdbx_ptnr2_PDB_ins_code
conn1 disulf A 10 CYS SG A 20 CYS SG A ?
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Create structure with only residue 10 (no insertion) and 20
        st = gemmi.Structure()
        model = gemmi.Model("1")
        chain = gemmi.Chain("A")

        res1 = gemmi.Residue()
        res1.name = "CYS"
        res1.seqid = gemmi.SeqId("10")  # No insertion code
        res1.add_atom(_make_atom("SG"))
        chain.add_residue(res1)

        res2 = gemmi.Residue()
        res2.name = "CYS"
        res2.seqid = gemmi.SeqId("20")
        res2.add_atom(_make_atom("SG"))
        chain.add_residue(res2)

        model.add_chain(chain)
        st.add_model(model)

        sync_conn_categories(block, st)

        conn_ids = [row[0] for row in block.find("_struct_conn.", ["id"])]

        # conn1 should be removed (10A doesn't exist, only 10)
        assert "conn1" not in conn_ids
