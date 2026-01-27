"""Chain category synchronization for mmCIF files."""

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync._utils import keep_rows_by_column


def sync_chain_categories(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Synchronize chain-related categories to match current structure.

    Keeps only rows for chains that exist in the structure.

    Categories synced:
    - _struct_asym

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    # Sync _struct_asym (chain ID is in "id" column)
    keep_rows_by_column(block, "_struct_asym.", "id", info.chain_ids)
