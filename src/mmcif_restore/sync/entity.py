"""Entity category synchronization for mmCIF files."""

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync._utils import keep_rows_by_column


def sync_entity_categories(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Synchronize entity-related categories to match current structure.

    Keeps only rows for entities that exist in the structure.

    Categories synced:
    - _entity
    - _entity_poly
    - _pdbx_entity_nonpoly

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    # Sync _entity (entity ID is in "id" column)
    keep_rows_by_column(block, "_entity.", "id", info.entity_ids)

    # Sync _entity_poly (entity_id column)
    keep_rows_by_column(block, "_entity_poly.", "entity_id", info.entity_ids)

    # Sync _pdbx_entity_nonpoly (entity_id column)
    keep_rows_by_column(block, "_pdbx_entity_nonpoly.", "entity_id", info.entity_ids)
