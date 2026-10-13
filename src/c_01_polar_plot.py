#-------------------------------------
#Circular (polar) peak plot of a per-base genome coverage array.
#Not wired into the CLI; can be fed with a coverage array from
#c_map_dist.compute_coverage(...)[contig].
#-------------------------------------
import logging

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)


def main_polar(genome_array, output="polar_plot.pdf", n_bins=99, top_n=5):
    """Save a polar plot of the max coverage in `n_bins` equal genome windows to `output`."""
    genome_array = np.asarray(genome_array)
    genome_length = len(genome_array)
    if genome_length == 0:
        raise ValueError("main_polar: empty genome array")

    n_bins = max(1, min(n_bins, genome_length))
    locations = np.linspace(0, genome_length, n_bins + 1)
    edges = locations.astype(int)
    counts = np.array([genome_array[a:b].max() if b > a else 0
                       for a, b in zip(edges[:-1], edges[1:])] + [0], dtype=float)
    max_count = counts.max() if counts.max() > 0 else 1.0

    theta = 2 * np.pi * locations / genome_length
    colors = plt.cm.viridis(counts / max_count)

    k = min(top_n, n_bins)
    top_idx = np.argsort(counts)[-k:][::-1]
    top_locations, top_counts, top_theta = locations[top_idx], counts[top_idx], theta[top_idx]

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(12, 6))
    try:
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        for i in range(len(locations)):
            ax.plot([theta[i], theta[i]], [max_count - counts[i], max_count], color=colors[i], linewidth=1.5)

        step = 1_000_000 if genome_length >= 2_000_000 else max(1, genome_length // 5)
        tick_positions = np.arange(0, genome_length, step)
        unit, div = ("Mb", 1e6) if step >= 1_000_000 else ("kb", 1e3)
        theta_ticks = 2 * np.pi * tick_positions / genome_length
        ax.set_xticks(theta_ticks)
        ax.set_xticklabels([f"{pos / div:g} {unit}" for pos in tick_positions])
        for t in theta_ticks:
            ax.plot([t, t], [0, max_count], color='black', lw=1)

        for i in range(k):
            ax.plot(top_theta[i], max_count, 'bo', markersize=10)
            ax.text(top_theta[i], max_count - max_count / 10, f"P{i + 1}", ha='center', va='center',
                    fontsize=10, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        ax2 = fig.add_axes([0.85, 0.1, 0.1, 0.8])
        ax2.axis('off')
        for i in range(k):
            ax2.text(0, 0.9 - i * 0.07, f"{i + 1}. Peak: {int(top_locations[i] / 1e3)} kbp, Count: {top_counts[i]:g}",
                     fontsize=12, ha='center')

        ax.set_title("Coverage peaks along the genome", va='bottom')
        ax.set_ylim(0, max_count)
        ax.yaxis.set_visible(False)

        sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, max_count))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.1, shrink=0.8)
        cbar.set_label("Max coverage per window")
        fig.savefig(output)
    finally:
        plt.close(fig)
    logger.info("Polar plot saved as %s", output)
    return output
