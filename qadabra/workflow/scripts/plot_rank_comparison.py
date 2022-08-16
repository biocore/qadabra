import logging

import numpy as np
import pandas as pd
import seaborn as sns

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Select, HoverTool
from bokeh.models.callbacks import CustomJS


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

logging.captureWarnings(True)
logging.getLogger("py.warnings").addHandler(fh)

output_file(filename=snakemake.output[0], title="Rank Comparison")

diff_df = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
rank_df = (
    diff_df.rank(ascending=False)
    .rename(columns={f"{x}": f"{x}_rank" for x in diff_df.columns})
    .join(diff_df)
    .reset_index()
)

tool_list = diff_df.columns.tolist()
tool_rank_list = [x + "_rank" for x in tool_list]

# Rank Comparison
logger.info("Creating rank comparison panel...")
chosen_tool_1 = Select(options=tool_rank_list, value=tool_rank_list[0], title="x-axis")
chosen_tool_2 = Select(options=tool_rank_list, value=tool_rank_list[1], title="y-axis")

rank_df["x"] = rank_df[chosen_tool_1.value]
rank_df["y"] = rank_df[chosen_tool_2.value]

source = ColumnDataSource(rank_df)

hover = HoverTool(mode="mouse", names=["points"], attachment="below")
hover.tooltips = (
    [("Feature ID", "@feature_id")] +
    [(f"{x}", f"@{x}") for x in rank_df.columns if x not in ["feature_id", "x", "y"]]
)

tools = [
    hover,
    "box_zoom",
    "reset",
    "pan",
    "save"
]

plot = figure(tools=tools)
plot.circle(
    source=source,
    x="x",
    y="y",
    color="lightgreen",
    line_width=0.5,
    line_color="black",
    size=10,
    name="points"
)
plot.line(
    x=np.arange(0, rank_df.shape[0]),
    y=np.arange(0, rank_df.shape[0]),
    line_width=3,
    line_dash="dashed",
    color="black"
)

tool_callback = CustomJS(args=dict(source=source, plot=plot), code="""
const data = source.data[cb_obj.value];
if (cb_obj.title == 'x-axis') {
    plot.below[0].axis_label = cb_obj.value;
    source.data.x = data;
} else {
    plot.left[0].axis_label = cb_obj.value;
    source.data.y = data;
}

source.change.emit();
""")
chosen_tool_1.js_on_change("value", tool_callback)
chosen_tool_2.js_on_change("value", tool_callback)

plot.title.text = "Differential Rank Comparison"
plot.title.text_font_size = "20pt"

plot.xaxis.axis_label = chosen_tool_1.value
plot.yaxis.axis_label = chosen_tool_2.value

for ax in [plot.xaxis, plot.yaxis]:
    ax.axis_label_text_font_size = "14pt"
    ax.axis_label_text_font_style = "normal"
    ax.major_label_text_font_size = "10pt"

layout = row(column(chosen_tool_1, chosen_tool_2), plot)
save(layout)
logger.info(f"Saved rank comparisons to {snakemake.output[0]}")
