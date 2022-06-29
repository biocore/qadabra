import logging

import numpy as np
import pandas as pd
import seaborn as sns

from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, Arrow, VeeHead, HoverTool


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

output_file(filename=snakemake.output[0], title="PCA")

feat_df = pd.read_table(snakemake.input["features"], sep="\t", index_col=0)
tool_df = pd.read_table(snakemake.input["tools"], sep="\t", index_col=0)
prop_exp = pd.read_table(snakemake.input["prop_exp"], sep="\t", index_col=0)
prop_exp = prop_exp.squeeze()

pc_cols = tool_df.columns.tolist()
tool_list = tool_df.index.tolist()

feat_df = feat_df.reset_index()
tool_df = tool_df.reset_index()

min_max_vals = pd.concat(
    [tool_df[pc_cols], feat_df[pc_cols]]
).agg(["min", "max"])

x_min, x_max = min_max_vals["PC1"]
y_min, y_max = min_max_vals["PC2"]
x_range = x_max - x_min
y_range = y_max - y_min

offset = 0.05
x_min = x_min - offset*x_range
x_max = x_max + offset*x_range

y_min = y_min - offset*y_range
y_max = y_max + offset*y_range

arrow_palette = dict(zip(
    tool_list,
    sns.color_palette("Dark2", len(tool_list)).as_hex()
))

source = ColumnDataSource(feat_df)

hover_points = HoverTool(mode="mouse", names=["points"], attachment="below")
hover_points.tooltips = [(f"{x}", f"@{x}") for x in feat_df.columns]

tools = [
    "box_zoom",
    "reset",
    "pan",
    "save"
]

plot = figure(tools=tools + [hover_points],
              x_range=(x_min, x_max), y_range=(y_min, y_max))
plot.circle(
    source=source,
    x="PC1",
    y="PC2",
    color="lightgreen",
    line_width=0.5,
    line_color="black",
    size=10,
    name="points"
)

for i, row in tool_df.iterrows():
    tool_name = row["tool"]
    color = arrow_palette[tool_name]
    plot.line(x=[0, row["PC1"]], y=[0, row["PC2"]], line_width=3, line_color=color,
              legend_label=tool_name)
    # plot arrow on top since we can't add arrows to legend
    arr = Arrow(
        end=VeeHead(line_width=3, line_color=color, fill_color=color),
        x_start=0, y_start=0, x_end=row["PC1"], y_end=row["PC2"],
        line_color=color, line_width=4,
        name="arrows"
    )
    plot.add_layout(arr)

prop_exp_labels = [
    f"{i} ({x*100:.2f}%)"
    for i, x in prop_exp.iteritems()
]
plot.xaxis.axis_label = prop_exp_labels[0]
plot.yaxis.axis_label = prop_exp_labels[1]
plot.xaxis.major_label_text_font_size = "10pt"
plot.yaxis.major_label_text_font_size = "10pt"

for ax in [plot.xaxis, plot.yaxis]:
    ax.axis_label_text_font_size = "14pt"
    ax.axis_label_text_font_style = "normal"
    ax.major_label_text_font_size = "10pt"

save(plot)
logger.info(f"Saved PCA to {snakemake.output[0]}")
