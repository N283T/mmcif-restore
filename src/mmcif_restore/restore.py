"""Main restore function for mmcif-restore."""

import logging
from collections.abc import Callable
from pathlib import Path

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync.chain import sync_chain_categories
from mmcif_restore.sync.conn import sync_conn_categories
from mmcif_restore.sync.entity import sync_entity_categories
from mmcif_restore.sync.modres import sync_modres_categories
from mmcif_restore.sync.scheme import sync_scheme_categories

logger = logging.getLogger(__name__)

# Registry of category sync handlers (takes StructureInfo)
CATEGORY_SYNC_HANDLERS: dict[str, Callable[[gemmi.cif.Block, StructureInfo], None]] = {
    "_entity.": sync_entity_categories,
    "_entity_poly.": sync_entity_categories,
    "_entity_poly_seq.": sync_entity_categories,
    "_pdbx_entity_nonpoly.": sync_entity_categories,
    "_struct_asym.": sync_chain_categories,
    "_pdbx_nonpoly_scheme.": sync_scheme_categories,
    "_pdbx_poly_seq_scheme.": sync_scheme_categories,
    "_pdbx_branch_scheme.": sync_scheme_categories,
}

# Categories that need structure instead of info (residue/atom-level sync)
STRUCTURE_SYNC_HANDLERS: dict[
    str, Callable[[gemmi.cif.Block, gemmi.Structure], None]
] = {
    "_struct_conn.": sync_conn_categories,
    "_pdbx_struct_mod_residue.": sync_modres_categories,
}


class RestoreError(Exception):
    """Error during CIF restoration."""


def restore_categories(
    edited_path: str | Path,
    reference_path: str | Path,
    categories: list[str],
) -> gemmi.cif.Document:
    """Restore mmCIF categories from reference CIF to an edited structure.

    This function takes a CIF file that has been edited (e.g., by ChimeraX,
    PyMOL, or Gemmi) and restores specified categories from the original
    reference CIF, synchronized to match the current structure.

    Args:
        edited_path: Path to the edited CIF file (can have minimal categories)
        reference_path: Path to the original reference CIF with all categories
        categories: List of category prefixes to restore
                   (e.g., ["_entity.", "_struct_asym.", "_struct_conn."])

    Returns:
        A new CIF Document with the restored and synchronized categories

    Raises:
        RestoreError: If files cannot be read or are invalid

    Example:
        # After editing a structure in ChimeraX and saving as minimal CIF:
        doc = restore_categories(
            "edited.cif",
            "original.cif.gz",
            categories=["_entity.", "_struct_asym.", "_struct_conn."]
        )
        doc.write_file("restored.cif")
    """
    edited_path = Path(edited_path)
    reference_path = Path(reference_path)

    # Load reference CIF (needed for entity mapping)
    try:
        ref_doc = gemmi.cif.read(str(reference_path))
    except Exception as e:
        raise RestoreError(
            f"Failed to read reference CIF '{reference_path}': {e}"
        ) from e

    if len(ref_doc) == 0:
        raise RestoreError(f"Reference CIF '{reference_path}' contains no data blocks")

    ref_block = ref_doc[0]

    # Read structure from edited CIF and extract info using reference
    try:
        structure = gemmi.read_structure(str(edited_path))
    except Exception as e:
        raise RestoreError(f"Failed to read edited CIF '{edited_path}': {e}") from e

    # Validate structure is not empty
    if len(structure) == 0 or all(len(model) == 0 for model in structure):
        raise RestoreError(f"Edited CIF '{edited_path}' contains no atoms")

    info = StructureInfo.from_structure_with_reference(structure, ref_block)

    # Load edited CIF as base document
    try:
        doc = gemmi.cif.read(str(edited_path))
    except Exception as e:
        raise RestoreError(f"Failed to read edited CIF '{edited_path}': {e}") from e

    if len(doc) == 0:
        raise RestoreError(f"Edited CIF '{edited_path}' contains no data blocks")

    block = doc[0]

    # Normalize category prefixes (ensure they end with ".")
    categories = [c if c.endswith(".") else c + "." for c in categories]

    # Restore each category from reference with sync
    for category_prefix in categories:
        _restore_and_sync_category(ref_block, block, category_prefix, info, structure)

    return doc


def _restore_and_sync_category(
    source_block: gemmi.cif.Block,
    target_block: gemmi.cif.Block,
    category_prefix: str,
    info: StructureInfo,
    structure: gemmi.Structure,
) -> None:
    """Restore a category from source to target, synced to structure info.

    Args:
        source_block: Reference CIF block to copy from
        target_block: Edited CIF block to copy to
        category_prefix: Category prefix (e.g., "_struct_conn.")
        info: Current structure information for filtering
        structure: Current structure for atom-level filtering
    """
    # Find the loop in source
    source_loop = None
    for item in source_block:
        if item.loop is not None:
            tags = list(item.loop.tags)
            if tags and tags[0].startswith(category_prefix):
                source_loop = item.loop
                break

    if source_loop is None:
        logger.debug("Category %s not found in reference CIF", category_prefix)
        return

    # Copy to target
    tags = list(source_loop.tags)
    col_names = [tag.replace(category_prefix, "") for tag in tags]

    new_loop = target_block.init_mmcif_loop(category_prefix, col_names)

    values = list(source_loop.values)
    width = len(tags)
    num_rows = len(values) // width
    columns = [
        [values[row * width + col] for row in range(num_rows)] for col in range(width)
    ]
    new_loop.set_all_values(columns)

    # Apply sync based on category type
    if category_prefix in STRUCTURE_SYNC_HANDLERS:
        # Residue/atom-level sync (needs structure)
        handler = STRUCTURE_SYNC_HANDLERS[category_prefix]
        handler(target_block, structure)
    else:
        # Chain/entity-level sync (uses info)
        handler = CATEGORY_SYNC_HANDLERS.get(category_prefix)
        if handler:
            handler(target_block, info)
        else:
            logger.debug("No sync handler for category %s", category_prefix)
