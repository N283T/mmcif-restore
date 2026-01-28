"""Extract current entity/chain information from a Gemmi Structure."""

from dataclasses import dataclass

import gemmi


@dataclass(frozen=True)
class StructureInfo:
    """Current entity and chain information from a structure.

    Attributes:
        entity_ids: Set of entity IDs currently in the structure
        chain_ids: Set of chain IDs (label_asym_id/subchain) currently in the structure
        auth_chain_ids: Set of auth chain IDs (auth_asym_id/PDB chain name)
    """

    entity_ids: frozenset[str]
    chain_ids: frozenset[str]
    auth_chain_ids: frozenset[str]

    @classmethod
    def from_structure(cls, structure: gemmi.Structure) -> "StructureInfo":
        """Extract entity and chain info by scanning actual residues.

        Args:
            structure: Gemmi Structure to analyze

        Returns:
            StructureInfo with current entity/chain IDs
        """
        entity_ids: set[str] = set()
        chain_ids: set[str] = set()
        auth_chain_ids: set[str] = set()

        # Build subchain to entity mapping from structure.entities
        subchain_to_entity: dict[str, str] = {}
        for entity in structure.entities:
            for subchain in entity.subchains:
                subchain_to_entity[subchain] = entity.name

        for model in structure:
            for chain in model:
                has_residues = False
                for residue in chain:
                    if residue.subchain:
                        chain_ids.add(residue.subchain)
                        entity_id = subchain_to_entity.get(residue.subchain)
                        if entity_id:
                            entity_ids.add(entity_id)
                        has_residues = True
                if has_residues:
                    auth_chain_ids.add(chain.name)

        return cls(
            entity_ids=frozenset(entity_ids),
            chain_ids=frozenset(chain_ids),
            auth_chain_ids=frozenset(auth_chain_ids),
        )

    @classmethod
    def from_structure_with_reference(
        cls,
        structure: gemmi.Structure,
        reference_block: gemmi.cif.Block,
    ) -> "StructureInfo":
        """Extract info from structure using reference CIF for entity mapping.

        This is useful when the structure was read from a minimal CIF that
        doesn't have entity information. The reference CIF provides the
        chain-to-entity mapping.

        Args:
            structure: Gemmi Structure (may lack entity info)
            reference_block: Reference CIF block with _struct_asym

        Returns:
            StructureInfo with current entity/chain IDs
        """
        chain_ids: set[str] = set()
        auth_chain_ids: set[str] = set()

        # Get chain IDs from structure's actual residues
        for model in structure:
            for chain in model:
                has_residues = False
                for residue in chain:
                    if residue.subchain:
                        chain_ids.add(residue.subchain)
                        has_residues = True
                if has_residues:
                    auth_chain_ids.add(chain.name)

        # Get entity mapping from reference _struct_asym
        chain_to_entity: dict[str, str] = {}
        struct_asym = reference_block.find("_struct_asym.", ["id", "entity_id"])
        for row in struct_asym:
            chain_to_entity[row[0]] = row[1]

        # Map chain IDs to entity IDs
        entity_ids: set[str] = set()
        for chain_id in chain_ids:
            entity_id = chain_to_entity.get(chain_id)
            if entity_id:
                entity_ids.add(entity_id)

        return cls(
            entity_ids=frozenset(entity_ids),
            chain_ids=frozenset(chain_ids),
            auth_chain_ids=frozenset(auth_chain_ids),
        )
