import logging
import os
import pathlib
from pkg_resources import resource_filename
import shutil

import biom
import click
import pandas as pd

from qadabra import __version__
from qadabra.utils import _validate_input

SNKFILE_TEXT = """from pkg_resources import resource_filename

from snakemake.utils import min_version
min_version("6.0")

qadabra_snakefile = resource_filename("qadabra", "workflow/Snakefile")
configfile: "config/config.yaml"

module qadabra:
    snakefile:
        qadabra_snakefile
    config:
        config

use rule * from qadabra
"""


@click.group()
@click.version_option(__version__)
def qadabra():
    """Differential abundance workflow"""
    pass


@qadabra.command()
@click.option(
    "--workflow-dest",
    type=click.Path(exists=True),
    required=True,
    help="Location of workflow",
    default="."
)
@click.option(
    "--table",
    type=click.Path(exists=True),
    required=True,
    help="Feature table in BIOM format"
)
@click.option(
    "--metadata",
    type=click.Path(exists=True),
    required=True,
    help="Metadata in TSV format"
)
@click.option(
    "--tree",
    type=click.Path(exists=True),
    required=False,
    help="Phylogenetic tree in Newick format"
)
@click.option(
    "--name",
    type=str,
    required=True,
    help="Name of dataset"
)
@click.option(
    "--factor-name",
    type=str,
    required=True,
    help="Name of factor grouping in metadata"
)
@click.option(
    "--target-level",
    type=str,
    required=True,
    help="Grouping level on which to perform differential abundance"
)
@click.option(
    "--reference-level",
    type=str,
    required=True,
    help="Grouping level to use as reference"
)
@click.option(
    "--confounder",
    type=str,
    required=False,
    multiple=True,
    help="Confounder variable to consider (can provide multiple)"
)
@click.option(
    "--validate-input/--no-validate-input",
    show_default=True,
    default=True,
    help="Whether to validate input prior to adding dataset"
)
@click.option(
    "--verbose",
    is_flag=True,
    show_default=True,
    default=False,
    help="Whether to output progress to console"
)
def add_dataset(
    workflow_dest,
    table,
    metadata,
    tree,
    name,
    factor_name,
    target_level,
    reference_level,
    confounder,
    validate_input,
    verbose
):
    """Add dataset on which to run Qadabra"""
    wflow_path = pathlib.Path(workflow_dest)
    wflow_dir = wflow_path / "workflow"
    cfg_dir = wflow_path / "config"

    if not wflow_dir.exists:
        raise ValueError("Workflow has not been created!")

    if not cfg_dir.exists:
        raise ValueError("Config directory has not been created!")

    dataset_sheet = cfg_dir / "datasets.tsv"
    logger = logging.getLogger(__name__)
    log_level = logging.INFO if verbose else logging.WARNING
    logger.setLevel(log_level)
    sh = logging.StreamHandler()
    sh.setLevel(log_level)
    formatter = logging.Formatter(
        "[%(asctime)s - %(levelname)s] :: %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    if validate_input:
        logger.info("Validating input...")
        _validate_input(logger, table, metadata, factor_name, target_level,
                        reference_level, tree, confounder)

    dataset_sheet = pathlib.Path(dataset_sheet)
    new_ds = pd.Series({
        "table": pathlib.Path(table).resolve(),
        "metadata": pathlib.Path(metadata).resolve(),
        "factor_name": factor_name,
        "target_level": target_level,
        "reference_level": reference_level,
    }, name=name).to_frame().T

    if tree is not None:
        new_ds["tree"] = pathlib.Path(tree).resolve()
    else:
        new_ds["tree"] = None

    if confounder:
        new_ds["confounders"] = ";".join(confounder)
    else:
        new_ds["confounders"] = None

    if dataset_sheet.exists():
        logger.info("Loading datasheet...")
        ds_sheet = pd.read_table(dataset_sheet, sep="\t", index_col=0)
        if name in ds_sheet.index:
            raise ValueError(f"{name} already exists in dataset sheet!")
        ds_sheet = pd.concat([ds_sheet, new_ds], axis=0)
    else:
        logger.info("Dataset does not exist. Creating...")
        ds_sheet = new_ds

    ds_sheet.to_csv(dataset_sheet, sep="\t", index=True)
    logger.info(f"Saved dataset sheet to {dataset_sheet}")


@qadabra.command()
@click.option(
    "--workflow-dest",
    type=click.Path(exists=False),
    default="."
)
def create_workflow(workflow_dest):
    """Create new Qadabra workflow structure"""
    wflow_dest = pathlib.Path(workflow_dest)
    wflow_dir = wflow_dest / "workflow"
    cfg_dir = wflow_dest/ "config"
    os.makedirs(wflow_dir)
    os.makedirs(cfg_dir)

    cfg_file = resource_filename("qadabra", "config/config.yaml")
    shutil.copy(cfg_file, cfg_dir / "config.yaml")

    style_file = resource_filename("qadabra", "config/qadabra.mplstyle")
    shutil.copy(style_file, cfg_dir / "qadabra.mplstyle")

    snkfile_path = wflow_dir / "Snakefile"
    with open(snkfile_path, "w") as f:
        f.write(SNKFILE_TEXT)


if __name__ == "__main__":
    qadabra()
