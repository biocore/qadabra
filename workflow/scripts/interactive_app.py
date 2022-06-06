import logging

import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, Select, HoverTool, Arrow,
                          VeeHead, DataTable, TableColumn)
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Tabs, Panel


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

output_file(filename=snakemake.output[0], title="qadabra")

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
panel_1 = Panel(child=layout, title="Rank Comparisons")

# PCA
logger.info("Creating PCA panel...")
scaled_diff_values = StandardScaler().fit_transform(diff_df.values)

pc_cols = [f"PC{x+1}" for x in range(len(tool_list))]
pca = PCA(n_components=len(tool_list)).fit(scaled_diff_values)
feat_df = pd.DataFrame(
    pca.transform(scaled_diff_values),
    index=diff_df.index,
    columns=pc_cols
).join(rank_df.set_index("feature_id").drop(columns=["x", "y"])).reset_index()
max_vals = feat_df[pc_cols].max()

tool_df = pd.DataFrame(
    np.dot(pca.components_, np.diag(max_vals)),
    index=diff_df.columns,
    columns=pc_cols
).reset_index()
tool_df = tool_df.rename(columns={"index": "tool"})

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

prop_exp = pca.explained_variance_ratio_
prop_exp_labels = [f"PC{i+1} ({x*100:.2f}%)" for i, x in enumerate(prop_exp)]
plot.xaxis.axis_label = prop_exp_labels[0]
plot.yaxis.axis_label = prop_exp_labels[1]
plot.xaxis.major_label_text_font_size = "10pt"
plot.yaxis.major_label_text_font_size = "10pt"

for ax in [plot.xaxis, plot.yaxis]:
    ax.axis_label_text_font_size = "14pt"
    ax.axis_label_text_font_style = "normal"
    ax.major_label_text_font_size = "10pt"

panel_2 = Panel(child=plot, title="PCA")

# Table
logger.info("Creating table panel...")
source = ColumnDataSource(diff_df.reset_index())

columns = [
    TableColumn(field=x, title=x)
    for x in diff_df.reset_index().columns
]
data_table = DataTable(source=source, columns=columns,
                       autosize_mode="fit_viewport")
panel_3 = Panel(child=data_table, title="Table")

tabs = Tabs(tabs=[panel_1, panel_2, panel_3], tabs_location="above")
save(tabs)

logger.info(f"Saved to {snakemake.output[0]}")
