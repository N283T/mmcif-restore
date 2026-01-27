"""Command-line interface for mmcif-restore."""

from pathlib import Path
from typing import Annotated

import typer

from mmcif_restore import RestoreError, restore_categories

app = typer.Typer(
    name="mmcif-restore",
    help="Restore mmCIF categories from reference CIF to edited structures.",
    add_completion=False,
)


@app.command()
def main(
    edited_file: Annotated[
        Path,
        typer.Argument(
            help="Edited CIF file (e.g., output from ChimeraX/PyMOL)",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    reference_file: Annotated[
        Path,
        typer.Argument(
            help="Reference CIF file with all categories (can be gzipped)",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "-o",
            "--output",
            help="Output file path",
            writable=True,
        ),
    ],
    categories: Annotated[
        str,
        typer.Option(
            "-c",
            "--categories",
            help="Categories to restore (comma-separated)",
        ),
    ],
) -> None:
    """Restore mmCIF categories from reference to edited structure.

    This tool takes a CIF file that has been edited (e.g., by removing waters
    or ligands in ChimeraX/PyMOL) and restores specified categories from the
    original reference CIF, synchronized to match the current structure.

    Example usage:

        # Restore _entity and _struct_conn from original
        mmcif-restore edited.cif original.cif.gz -o output.cif \\
            -c "_entity.,_struct_conn."

        # Restore multiple categories
        mmcif-restore edited.cif original.cif.gz -o output.cif \\
            -c "_entity.,_struct_asym.,_struct_conn.,_pdbx_nonpoly_scheme."
    """
    # Parse categories
    cats = [c.strip() for c in categories.split(",") if c.strip()]
    if not cats:
        typer.echo("Error: At least one category must be specified", err=True)
        raise typer.Exit(1)

    # Validate output path
    if not output.parent.exists():
        typer.echo(
            f"Error: Output directory '{output.parent}' does not exist", err=True
        )
        raise typer.Exit(1)

    typer.echo(f"Edited CIF: {edited_file}")
    typer.echo(f"Reference CIF: {reference_file}")
    typer.echo(f"Categories to restore: {', '.join(cats)}")

    # Restore
    try:
        doc = restore_categories(
            edited_file,
            reference_file,
            categories=cats,
        )
    except RestoreError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1) from e

    # Write output
    doc.write_file(str(output))
    typer.echo(f"Written to {output}")


if __name__ == "__main__":
    app()
