import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bokeh.plotting import figure, save
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.io import output_file

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

plt.style.use(snakemake.config["stylesheet"])

logger.info("Loading differentials...")
diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0).squeeze()
values = diffs.sort_values().values

feature_names = diffs.sort_values().index

logger.info("Plotting differentials...")

# Create a ColumnDataSource to store the data
source = ColumnDataSource(data=dict(x=np.arange(len(feature_names)), feature_names=feature_names, values=values))

# Create the figure
# p = figure(x_range=(0, len(values)), height=400, width=600, toolbar_location=None)
p = figure(x_range=(0, len(values)), toolbar_location=None)

# Plot the bars
p.vbar(x='x', top='values', width=0.8, source=source)

# Add HoverTool to display information on hover
hover = HoverTool(tooltips=[("Feature Name", "@feature_names")])
p.add_tools(hover)

# Set plot properties
p.xaxis.axis_label = "Rank"
p.xgrid.grid_line_color = 'gray'
p.ygrid.grid_line_color = 'gray'

# Configure the title
tool_name = snakemake.wildcards["tool"]
p.title.text = f"{tool_name} Differentials"
p.title.align = 'center'
p.title.text_font_size = '18pt'

# Save plot
output_file(snakemake.output[0])
save(p)

logger.info(f"Saved interactive rank plot to {snakemake.output[0]}")

### the below code makes the non-interactive rank plots

# fig, ax = plt.subplots(1, 1)
# plt.bar(
#     np.arange(len(values)),
#     values,
#     width=1
# )
# tool_name = snakemake.wildcards["tool"]
# ax.set_title(f"{tool_name} Differentials")
# ax.set_xlabel("Rank")
# ax.grid("both")
# plt.savefig(snakemake.output[0])
# logger.info(f"Saved to {snakemake.output[0]}")


# logger.info("Plotting differentials...")
# fig, ax = plt.subplots(1, 1)
# plt.bar(
#     np.arange(len(values)),
#     values,
#     width=1
# )
# tool_name = snakemake.wildcards["tool"]
# ax.set_title(f"{tool_name} Differentials")
# ax.set_xlabel("Rank")
# ax.grid("both")
# plt.savefig(snakemake.output[0])
# logger.info(f"Saved to {snakemake.output[0]}")

