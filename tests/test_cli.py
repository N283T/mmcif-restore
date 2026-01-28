"""Tests for CLI interface."""

from pathlib import Path

import gemmi
from typer.testing import CliRunner

from mmcif_restore.cli import app

runner = CliRunner()


class TestCli:
    """Tests for mmcif-restore CLI."""

    def test_restore_entity(self, sample_cif_file: Path, tmp_path: Path) -> None:
        """Test restoring _entity from reference."""
        # Create minimal CIF
        structure = gemmi.read_structure(str(sample_cif_file))
        structure.remove_waters()

        groups = gemmi.MmcifOutputGroups(False)
        groups.atoms = True
        groups.entry = True
        minimal_doc = structure.make_mmcif_document(groups)
        minimal_path = tmp_path / "minimal.cif"
        minimal_doc.write_file(str(minimal_path))

        output = tmp_path / "restored.cif"

        result = runner.invoke(
            app,
            [
                str(minimal_path),
                str(sample_cif_file),
                "-o",
                str(output),
                "-c",
                "_entity.",
            ],
        )

        assert result.exit_code == 0
        assert output.exists()

        # Verify _entity restored without water
        doc = gemmi.cif.read(str(output))
        block = doc[0]
        entity_types = [row[1] for row in block.find("_entity.", ["id", "type"])]
        assert "polymer" in entity_types
        assert "water" not in entity_types

    def test_restore_multiple_categories(
        self, sample_cif_file: Path, tmp_path: Path
    ) -> None:
        """Test restoring multiple categories."""
        # Create minimal CIF
        structure = gemmi.read_structure(str(sample_cif_file))
        structure.remove_waters()

        groups = gemmi.MmcifOutputGroups(False)
        groups.atoms = True
        groups.entry = True
        minimal_doc = structure.make_mmcif_document(groups)
        minimal_path = tmp_path / "minimal.cif"
        minimal_doc.write_file(str(minimal_path))

        output = tmp_path / "restored.cif"

        result = runner.invoke(
            app,
            [
                str(minimal_path),
                str(sample_cif_file),
                "-o",
                str(output),
                "-c",
                "_entity.,_struct_asym.",
            ],
        )

        assert result.exit_code == 0
        assert output.exists()

        # Verify both categories restored
        # 5i55.cif has 4 entities/chains, after removing water: 3 each
        doc = gemmi.cif.read(str(output))
        block = doc[0]
        entity_count = len(list(block.find("_entity.", ["id"])))
        asym_count = len(list(block.find("_struct_asym.", ["id"])))
        assert entity_count == 3
        assert asym_count == 3

    def test_error_no_categories(self, sample_cif_file: Path, tmp_path: Path) -> None:
        """Test error when no categories specified."""
        output = tmp_path / "output.cif"

        result = runner.invoke(
            app,
            [
                str(sample_cif_file),
                str(sample_cif_file),
                "-o",
                str(output),
                "-c",
                "",
            ],
        )

        assert result.exit_code == 1
        assert "At least one category" in result.output

    def test_error_nonexistent_file(self, tmp_path: Path) -> None:
        """Test error when input file doesn't exist."""
        output = tmp_path / "output.cif"

        result = runner.invoke(
            app,
            [
                "/nonexistent/file.cif",
                "/nonexistent/ref.cif",
                "-o",
                str(output),
                "-c",
                "_entity.",
            ],
        )

        assert result.exit_code != 0
