import logging
import pathlib
from typing import List

import biom
import click
import pandas as pd

from qadabra import __version__


@click.group()
@click.version_option(__version__)
def qadabra():
    """Differential abundance workflow"""
    pass


@qadabra.command()
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
    "--validate-input",
    is_flag=True,
    show_default=True,
    default=True
)
@click.option(
    "--dataset-sheet",
    type=click.Path(),
    required=False,
    default="config/datasets.tsv",
    help="Path to dataset sheet for Qadabra (created if not found)"
)
@click.option(
    "--verbose",
    is_flag=True,
    show_default=True,
    default=False,
    help="Whether to output progress to console"
)
def add_dataset(
    table,
    metadata,
    tree,
    name,
    factor_name,
    target_level,
    reference_level,
    confounder,
    validate_input,
    dataset_sheet,
    verbose
):
    """Add dataset on which to run Qadabra"""
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
        "table": table,
        "metadata": metadata,
        "tree": tree,
        "factor_name": factor_name,
        "target_level": target_level,
        "reference_level": reference_level,
        "confounders": ";".join(confounder),
    }, name=name).to_frame().T

    if dataset_sheet.exists():
        logger.info("Loading datasheet...")
        ds_sheet = pd.read_table(dataset_sheet, sep="\t", index_col=0)
        if name in ds_sheet.index:
            raise ValueError(f"{name} already exists in dataset sheet!")
        ds_sheet = pd.concat([ds_sheet, new_ds], axis=0)
    else:
        logger.info("Dataset does not exist. Creating...")
        ds_sheet = new_ds

    print(ds_sheet)
    ds_sheet.to_csv(dataset_sheet, sep="\t", index=True)
    logger.info(f"Saved dataset sheet to {dataset_sheet}")


def _validate_input(
    logger: logging.Logger,
    table: pathlib.Path,
    metadata: pathlib.Path,
    factor_name: str,
    target_level: str,
    reference_level: str,
    tree: pathlib.Path = None,
    confounders: List[str] = None,
):
    logger.info("Loading metadata...")
    md = pd.read_table(metadata, sep="\t", index_col=0)

    if factor_name not in md.columns:
        raise ValueError(f"{factor_name} not found in metadata!")

    logger.info("Making sure factor & levels are all present in metadata...")
    factor_counts = md[factor_name].value_counts()
    logger.info(f"Factor counts: \n{factor_counts}")
    for level in [target_level, reference_level]:
        if level not in factor_counts:
            raise ValueError(f"{level} not found in {factor_name} values!")

    logger.info("Making sure confounders are all metadata columns...")
    for conf in confounders:
        if conf not in md.columns:
            raise ValueError(f"{conf} not found in metadata!")

    logger.info("Loading table...")
    tbl = biom.load_table(table)
    logger.info(f"Table shape: {tbl.shape}")
    tbl_idx = set(tbl.ids())
    md_idx = set(md.index)

    if not tbl_idx.issubset(md_idx):
        raise ValueError("Table IDs are not a subset of metadata IDs!")

    if tbl_idx != md_idx:
        logger.warn("Table IDs and metadata IDs are not exactly the same.")

    logger.info("Looking for completely discriminatory taxa...")
    tbl_df = tbl.to_dataframe(dense=True).T
    joint_df = tbl_df.join(md)
    gb = joint_df.groupby(factor_name).sum()
    feat_presence = gb.apply(lambda x: x.all())
    if not feat_presence.all():
        raise ValueError(
            "Some taxa in the table perfectly discriminate factor groups. "
            "Please filter out these taxa before running Qadabra."
        )

    if tree:
        from bp import parse_newick, to_skbio_treenode

        logger.info("Reading phylogenetic tree...")
        with open(tree) as f:
            tr = parse_newick(f.readline())
        tr = to_skbio_treenode(tree)
        tips = set(tr.tips())

        if not tips.issubset(tbl.ids("observation")):
            raise ValueError("Tree tips are not a subset of table features!")


if __name__ == "__main__":
    qadabra()
