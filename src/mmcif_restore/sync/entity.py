"""Entity category synchronization for mmCIF files."""

import logging
from typing import Any

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync._utils import keep_rows_by_column

logger = logging.getLogger(__name__)


def sync_entity_categories(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Synchronize entity-related categories to match current structure.

    Keeps only rows for entities that exist in the structure.

    Categories synced:
    - _entity
    - _entity_poly
    - _entity_poly_seq
    - _pdbx_entity_nonpoly

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    # Sync _entity (entity ID is in "id" column)
    keep_rows_by_column(block, "_entity.", "id", info.entity_ids)

    # Sync _entity_poly (entity_id column)
    keep_rows_by_column(block, "_entity_poly.", "entity_id", info.entity_ids)

    # Sync _entity_poly.pdbx_strand_id (comma-separated auth chain IDs)
    _sync_strand_ids(block, info)

    # Sync _entity_poly_seq (entity_id column)
    keep_rows_by_column(block, "_entity_poly_seq.", "entity_id", info.entity_ids)

    # Sync _pdbx_entity_nonpoly (entity_id column)
    keep_rows_by_column(block, "_pdbx_entity_nonpoly.", "entity_id", info.entity_ids)


def _sync_strand_ids(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Sync pdbx_strand_id in _entity_poly to match current auth chain IDs.

    Each pdbx_strand_id value is a comma-separated list of auth_asym_id values.
    Removes IDs no longer present in the structure.

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    try:
        data: dict[str, list[Any]] = block.get_mmcif_category("_entity_poly.", raw=True)
    except RuntimeError:
        logger.debug("Category _entity_poly. not found")
        return

    if not data or "pdbx_strand_id" not in data:
        return

    strand_ids = data["pdbx_strand_id"]
    updated = False
    new_strand_ids: list[str] = []

    for value in strand_ids:
        chains = [c.strip() for c in value.split(",")]
        filtered = [c for c in chains if c in info.auth_chain_ids]
        new_value = ",".join(filtered) if filtered else "?"
        new_strand_ids.append(new_value)
        if new_value != value:
            updated = True

    if not updated:
        return

    filtered_data = {**data, "pdbx_strand_id": new_strand_ids}
    block.set_mmcif_category("_entity_poly.", filtered_data, raw=True)
