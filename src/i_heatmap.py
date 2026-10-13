#-------------------------------------
#Heatmap of IS distributions from a .rcp stats file (output of h_is_stats.py).
#
#IS are split into two panels by their maximum percentage over all conditions:
#  high panel: IS with > 15 % of the IS_DB insertions in at least one condition
#  low panel : IS with <= 15 % in every condition (own colour scale)
#Only non-empty panels are drawn (B42).
#Only IS_DB insertions are plotted: the Genome / Unannotated / Contamination counts
#are '#' metadata lines in the .rcp file and are ignored here.
#-------------------------------------
import logging
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

HIGH_THRESHOLD = 15


# Pastel teal colour map
def create_pastel_teal_cmap():
    cdict = {
        'red': [(0.0, 1.0, 1.0), (0.2, 0.8, 0.8), (0.8, 0.3, 0.3), (1.0, 0.3, 0.3)],
        'green': [(0.0, 1.0, 1.0), (0.2, 0.9, 0.9), (0.8, 0.7, 0.7), (1.0, 0.7, 0.7)],
        'blue': [(0.0, 1.0, 1.0), (0.2, 0.9, 0.9), (0.8, 0.8, 0.8), (1.0, 0.8, 0.8)],
    }
    return mcolors.LinearSegmentedColormap('pastel_teal_cmap', segmentdata=cdict, N=256)


# Pastel orange colour map
def create_pastel_orange_cmap(color):
    cdict = {
        'red': [(0.0, 1.0, 1.0), (0.2, 1.0, 1.0), (0.8, 1.0, 1.0), (1.0, color[0], color[0])],
        'green': [(0.0, 1.0, 1.0), (0.2, 0.8, 0.8), (0.8, 0.5, 0.5), (1.0, 0.4, 0.4)],
        'blue': [(0.0, 1.0, 1.0), (0.2, 0.5, 0.5), (0.8, 0.3, 0.3), (1.0, 0.1, 0.1)],
    }
    return mcolors.LinearSegmentedColormap('pastel_cmap', segmentdata=cdict, N=256)


# Wide (Species x Condition) -> long format
def prepare_data_for_facetgrid(df):
    df = df.copy()
    df.index.name = 'Species'
    return df.reset_index().melt(id_vars='Species', var_name='Condition', value_name='Percentage')


def prepare_combined_data(df):
    return prepare_data_for_facetgrid(df)


def save_combined_matrix(df, filename):
    df.pivot_table(index='Species', columns='Condition', values='Percentage').to_csv(filename, sep='\t')


def _draw_heatmap(ax, data_long, cmap, vmax, title, order=None):
    matrix = data_long.pivot(index='Species', columns='Condition', values='Percentage')
    if order:
        matrix = matrix.reindex(columns=[c for c in order if c in matrix.columns])
    sns.heatmap(matrix, cmap=cmap, vmin=0, vmax=vmax, annot=True, fmt='.1f', linewidths=0.5,
                linecolor='gray', ax=ax, cbar_kws={'label': 'Percentage of IS_DB insertions'})
    ax.set_title(title)


def plot_heatmaps(data_high, data_low, output, vmax_low, threshold=HIGH_THRESHOLD, order=None):
    """Draw the non-empty panels of the high/low split into `output` (B42)."""
    panels = []
    if not data_high.empty:
        panels.append((data_high, create_pastel_orange_cmap((1.0, 0.7, 0.3)), 100,
                       f"IS with > {threshold}% of IS_DB insertions in at least one condition"))
    if not data_low.empty:
        vmax = float(vmax_low) if vmax_low and vmax_low > 0 else 1.0
        panels.append((data_low, create_pastel_teal_cmap(), vmax,
                       f"IS with <= {threshold}% of IS_DB insertions in every condition "
                       f"(colour scale 0-{vmax:.1f}%)"))

    if not panels:
        logger.warning("No IS data to plot; writing an empty heatmap figure.")
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, "No IS insertions", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    else:
        heights = [max(2.5, 0.45 * d['Species'].nunique() + 1.5) for d, *_ in panels]
        fig, axes = plt.subplots(nrows=len(panels), figsize=(10, sum(heights)),
                                 gridspec_kw={'height_ratios': heights}, squeeze=False)
        for ax, (data, cmap, vmax, title) in zip(axes[:, 0], panels):
            _draw_heatmap(ax, data, cmap, vmax, title, order)
    try:
        fig.tight_layout()
        fig.savefig(output)
    finally:
        plt.close(fig)
    logger.info("Heatmap saved as %s", output)


def parse_rcp(file_path):
    """Parse a .rcp file into a long DataFrame (Species, Percentage, Condition).

    Returns (df, condition_order). '#' metadata lines, '?' count lines and
    malformed lines are skipped.
    """
    records = []
    order = []
    current_condition = None
    with open(file_path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('?'):
                continue
            if line.startswith(">"):
                current_condition = line[1:]
                if current_condition not in order:
                    order.append(current_condition)
            elif line == "---":
                current_condition = None
            elif current_condition:
                parts = [p.strip() for p in line.split(',')]
                try:
                    percentage = float(parts[2])
                except (IndexError, ValueError):
                    logger.warning("Skipping malformed line in %s: %r", file_path, line)
                    continue
                records.append(('_'.join((parts[0], parts[1])), percentage, current_condition))
    return pd.DataFrame(records, columns=["Species", "Percentage", "Condition"]), order


def main_heatmap(file_path, output, combined_matrix_file_path):
    if not os.path.exists(file_path):
        logger.warning("File '%s' does not exist.", file_path)
        return 0

    df, order = parse_rcp(file_path)

    if df.empty:
        heatmap_data = pd.DataFrame(index=pd.Index([], name="Species"), dtype=float)
    else:
        # pivot_table: a species listed twice in one condition is summed instead of crashing
        heatmap_data = df.pivot_table(index="Species", columns="Condition", values="Percentage",
                                      aggfunc="sum").fillna(0)
        heatmap_data = heatmap_data[[c for c in order if c in heatmap_data.columns]]

    row_max = heatmap_data.max(axis=1)
    heatmap_data_high = heatmap_data[row_max > HIGH_THRESHOLD]
    heatmap_data_low = heatmap_data[row_max <= HIGH_THRESHOLD]
    vmax_low = heatmap_data_low.max().max() if not heatmap_data_low.empty else 0

    plot_heatmaps(prepare_data_for_facetgrid(heatmap_data_high),
                  prepare_data_for_facetgrid(heatmap_data_low), output, vmax_low, order=order)

    save_combined_matrix(prepare_combined_data(heatmap_data), combined_matrix_file_path)
    return 1
