"""Scheme category synchronization for mmCIF files."""

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync._utils import keep_rows_by_column


def sync_scheme_categories(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Synchronize scheme-related categories to match current structure.

    Keeps only rows for chains that exist in the structure.

    Categories synced:
    - _pdbx_poly_seq_scheme
    - _pdbx_nonpoly_scheme
    - _pdbx_branch_scheme

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    keep_chain_ids = info.chain_ids

    # Sync _pdbx_poly_seq_scheme
    keep_rows_by_column(block, "_pdbx_poly_seq_scheme.", "asym_id", keep_chain_ids)

    # Sync _pdbx_nonpoly_scheme
    keep_rows_by_column(block, "_pdbx_nonpoly_scheme.", "asym_id", keep_chain_ids)

    # Sync _pdbx_branch_scheme (for branched entities like carbohydrates)
    keep_rows_by_column(block, "_pdbx_branch_scheme.", "asym_id", keep_chain_ids)
