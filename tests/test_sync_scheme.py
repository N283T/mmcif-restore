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

        # 5i55.cif: chains A (polymer), B, C (non-polymer), D (water)
        # Keep chains A, B, C (remove water chain D)
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3"]),
            chain_ids=frozenset(["A", "B", "C"]),
        )

        sync_scheme_categories(block, info)

        # Check _pdbx_nonpoly_scheme table
        nonpoly_scheme = block.find("_pdbx_nonpoly_scheme.", ["asym_id", "mon_id"])
        asym_ids = [row[0] for row in nonpoly_scheme]

        assert "D" not in asym_ids  # water chain removed
        assert "B" in asym_ids  # non-polymer still present
        assert "C" in asym_ids  # non-polymer still present

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

        # 5i55.cif has 4 chains
        info = StructureInfo(
            entity_ids=frozenset(["1", "2", "3", "4"]),
            chain_ids=frozenset(["A", "B", "C", "D"]),
        )

        sync_scheme_categories(block, info)

        # Count should be unchanged
        new_count = len(list(block.find("_pdbx_nonpoly_scheme.", ["asym_id"])))
        assert new_count == original_count


class TestSyncPolySeqScheme:
    """Tests for _pdbx_poly_seq_scheme synchronization."""

    def test_removes_poly_seq_scheme_rows(self) -> None:
        """Test that _pdbx_poly_seq_scheme rows are removed for removed chains."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.auth_seq_num
_pdbx_poly_seq_scheme.pdb_mon_id
_pdbx_poly_seq_scheme.auth_mon_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_ins_code
_pdbx_poly_seq_scheme.hetero
A 1 1 ALA 1 1 ALA ALA A . n
A 1 2 GLY 2 2 GLY GLY A . n
B 2 1 VAL 1 1 VAL VAL B . n
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Keep only chain A
        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        sync_scheme_categories(block, info)

        # Chain B rows should be removed
        scheme = block.find("_pdbx_poly_seq_scheme.", ["asym_id"])
        asym_ids = [row[0] for row in scheme]

        assert "B" not in asym_ids
        assert asym_ids.count("A") == 2

    def test_keeps_all_poly_seq_scheme_when_all_kept(self) -> None:
        """Test that all rows remain when all chains kept."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
A 1 1 ALA
A 1 2 GLY
B 2 1 VAL
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        info = StructureInfo(
            entity_ids=frozenset(["1", "2"]),
            chain_ids=frozenset(["A", "B"]),
        )

        original_count = len(list(block.find("_pdbx_poly_seq_scheme.", ["asym_id"])))
        sync_scheme_categories(block, info)
        new_count = len(list(block.find("_pdbx_poly_seq_scheme.", ["asym_id"])))

        assert new_count == original_count

    def test_handles_missing_poly_seq_scheme(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that sync works when _pdbx_poly_seq_scheme doesn't exist."""
        block = sample_cif_document[0]

        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        # Should not raise
        sync_scheme_categories(block, info)


class TestSyncBranchScheme:
    """Tests for _pdbx_branch_scheme synchronization."""

    def test_removes_branch_scheme_rows(self) -> None:
        """Test that _pdbx_branch_scheme rows are removed for removed chains."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_branch_scheme.asym_id
_pdbx_branch_scheme.entity_id
_pdbx_branch_scheme.num
_pdbx_branch_scheme.mon_id
_pdbx_branch_scheme.pdb_asym_id
_pdbx_branch_scheme.pdb_seq_num
_pdbx_branch_scheme.auth_asym_id
_pdbx_branch_scheme.auth_seq_num
_pdbx_branch_scheme.auth_mon_id
_pdbx_branch_scheme.hetero
C 3 1 NAG C 1 C 1 NAG n
C 3 2 GAL C 2 C 2 GAL n
D 4 1 MAN D 1 D 1 MAN n
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Keep only chain C
        info = StructureInfo(
            entity_ids=frozenset(["3"]),
            chain_ids=frozenset(["C"]),
        )

        sync_scheme_categories(block, info)

        # Chain D rows should be removed
        scheme = block.find("_pdbx_branch_scheme.", ["asym_id"])
        asym_ids = [row[0] for row in scheme]

        assert "D" not in asym_ids
        assert asym_ids.count("C") == 2

    def test_keeps_all_branch_scheme_when_all_kept(self) -> None:
        """Test that all rows remain when all chains kept."""
        cif_content = """\
data_TEST
#
loop_
_pdbx_branch_scheme.asym_id
_pdbx_branch_scheme.entity_id
_pdbx_branch_scheme.num
_pdbx_branch_scheme.mon_id
C 3 1 NAG
C 3 2 GAL
D 4 1 MAN
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        info = StructureInfo(
            entity_ids=frozenset(["3", "4"]),
            chain_ids=frozenset(["C", "D"]),
        )

        original_count = len(list(block.find("_pdbx_branch_scheme.", ["asym_id"])))
        sync_scheme_categories(block, info)
        new_count = len(list(block.find("_pdbx_branch_scheme.", ["asym_id"])))

        assert new_count == original_count

    def test_handles_missing_branch_scheme(
        self, sample_cif_document: gemmi.cif.Document
    ) -> None:
        """Test that sync works when _pdbx_branch_scheme doesn't exist."""
        block = sample_cif_document[0]

        info = StructureInfo(
            entity_ids=frozenset(["1"]),
            chain_ids=frozenset(["A"]),
        )

        # Should not raise
        sync_scheme_categories(block, info)
