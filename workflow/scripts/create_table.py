import logging

import pandas as pd

from bokeh.plotting import output_file, save
from bokeh.models import DataTable, TableColumn, ColumnDataSource


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

output_file(filename=snakemake.output[0], title="Table")
diff_df = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
diff_df = diff_df.reset_index()

source = ColumnDataSource(diff_df)
columns = [
    TableColumn(field=x, title=x)
    for x in diff_df.columns
]
data_table = DataTable(source=source, columns=columns,
                       autosize_mode="fit_viewport",
                       frozen_columns=1,
                       sizing_mode="stretch_height")
save(data_table)
logger.info(f"Saved table to {snakemake.output[0]}")
