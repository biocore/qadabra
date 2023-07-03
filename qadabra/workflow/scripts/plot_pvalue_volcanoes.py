import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from bokeh.plotting import figure, save
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.io import output_file

# Add logging 
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

# Read in tsv containing p-values and LFCs
logger.info("Loading p-values and LFCs...")
data = pd.read_table(snakemake.input[0], sep="\t", index_col=0).squeeze()
pval_col = snakemake.params["pval_col"]
diff_ab_col = snakemake.params["diff_ab_col"]

# Take negative log of pval_col column and add as new column
data['negative_log_pval'] = -np.log10(data[pval_col])

# Create a ColumnDataSource object
source_blue = ColumnDataSource(data=dict(
    diff_ab_col=data[diff_ab_col],
    pval=data['negative_log_pval'],
    feature=data.index
))

source_red = ColumnDataSource(data=dict(
    diff_ab_col=data[diff_ab_col],
    pval=data['negative_log_pval'],
    feature=data.index
))

# Create a figure
tool_name = snakemake.wildcards["tool"]
p = figure(title=f"{tool_name} Volcano Plot", x_axis_label='coefficient', y_axis_label='-log10(pval_col)')

p.title.text_font_size = '20pt'

# Define the desired y-axis range
y_min = data['negative_log_pval'].min() - 1  # Minimum value
y_max = data['negative_log_pval'].max() + 1  # Maximum value

# Set the y-axis range
p.y_range.start = y_min
p.y_range.end = y_max

# Define the desired x-axis range
x_min = data[diff_ab_col].min() - 1  # Minimum value
x_max = data[diff_ab_col].max() + 1  # Maximum value

# Set the y-axis range
p.x_range.start = x_min
p.x_range.end = x_max

# # Add a vertical line at x = -1
p.line(x=[-1, -1], y=[y_min, y_max], line_color='gray', line_width=1, line_dash='dashed')
p.line(x=[1, 1], y=[y_min, y_max], line_color='gray', line_width=1, line_dash='dashed')
p.line(x=[x_min,x_max], y=[-np.log10(0.05), -np.log10(0.05)], line_color='gray', line_width=1, line_dash='dashed')

# Add hover tool
hover = HoverTool(
    tooltips=[
        ('Feature', '@feature'),
        ('coefficient', '@diff_ab_col'),
        ('-log10(pval_col)', '@pval'),
    ]
)
p.add_tools(hover)

# Select points where y >= 1 and x >= -log10(0.05)
pos_sig_diff = (data[diff_ab_col] >= 1) & (data['negative_log_pval'] >= -np.log10(0.05))
# Select points where y <= -1 and x >= -log10(0.05)
neg_sig_diff = (data[diff_ab_col] <= -1) & (data['negative_log_pval'] >= -np.log10(0.05))

# Apply selection to the data source
source_blue.selected.indices = np.where(neg_sig_diff)[0]
# Apply selection to the data source for negative significant differences (blue)
p.circle('diff_ab_col', 'pval', size=5, fill_alpha=0.6, nonselection_line_color='black', nonselection_fill_color='black', source=source_blue, selection_color='blue')

# Apply selection to the data source
source_red.selected.indices = np.where(pos_sig_diff)[0]
# Apply selection to the data source for positive significant differences (red)
p.circle('diff_ab_col', 'pval', size=5, fill_alpha=0.6, nonselection_line_color='black', nonselection_fill_color='black', source=source_red, selection_color='red')


output_file(snakemake.output[0])
save(p)

logger.info(f"Saved volocano plot to {snakemake.output[0]}")
