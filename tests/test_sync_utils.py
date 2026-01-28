"""Tests for sync utility functions."""

import gemmi

from mmcif_restore.sync._utils import keep_rows_by_column


class TestKeepRowsByColumnEdgeCases:
    """Tests for edge cases in keep_rows_by_column function."""

    def test_empty_keep_values_no_changes(self) -> None:
        """Test that empty keep_values results in early return with no changes."""
        cif_content = """\
data_TEST
#
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Call with empty keep_values
        keep_rows_by_column(
            block,
            category_prefix="_entity.",
            filter_column="id",
            keep_values=frozenset(),  # Empty
        )

        # Data should be unchanged
        entity_ids = [row[0] for row in block.find("_entity.", ["id"])]
        assert entity_ids == ["1", "2"]

    def test_missing_category_no_error(self) -> None:
        """Test that missing category is handled gracefully without error."""
        cif_content = """\
data_TEST
#
_entry.id TEST
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Call with non-existent category - should not raise
        keep_rows_by_column(
            block,
            category_prefix="_nonexistent.",
            filter_column="id",
            keep_values=frozenset(["1"]),
        )

        # Block should be unchanged
        assert block.find_value("_entry.id") == "TEST"

    def test_missing_filter_column_no_error(self) -> None:
        """Test that missing filter column is handled gracefully."""
        cif_content = """\
data_TEST
#
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Call with non-existent filter column - should not raise
        keep_rows_by_column(
            block,
            category_prefix="_entity.",
            filter_column="nonexistent_column",
            keep_values=frozenset(["1"]),
        )

        # Data should be unchanged
        entity_ids = [row[0] for row in block.find("_entity.", ["id"])]
        assert entity_ids == ["1", "2"]

    def test_empty_category_data_no_error(self) -> None:
        """Test that empty category data is handled gracefully."""
        # Create a block with an empty loop
        cif_content = """\
data_TEST
#
loop_
_empty_cat.id
#
"""
        doc = gemmi.cif.read_string(cif_content)
        block = doc[0]

        # Call on empty category - should not raise
        keep_rows_by_column(
            block,
            category_prefix="_empty_cat.",
            filter_column="id",
            keep_values=frozenset(["1"]),
        )
