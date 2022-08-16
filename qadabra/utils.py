import logging
import pathlib
from typing import List

import biom
import pandas as pd


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
