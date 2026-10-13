#!/usr/bin/env python3
"""
Paper figures from the benchmark grid (simulate/run_grid.py -> grid_results.tsv).

    python simulate/plot_grid.py --results publication/benchmark/grid/grid_results.tsv \\
        [--out_dir publication/benchmark/grid/figures] [--dev_default evalue=1e-10,blast_task=megablast]

Figures (PDF + PNG each; mean over seeds, band/error bar = ± 1 SD):
  fig1_detection_vs_accuracy     read-level F1 (insert finder, annotation) and
                                 annotation recall vs mean read accuracy
  fig2_position_vs_accuracy      % insertion reads placed exactly / within 5 bp
  fig3_classification_vs_evalue  IS-name accuracy, novel-insert accuracy, % true
                                 insertions labelled Contamination / Unannotated,
                                 per E-value (blue ramp) and BLAST task (rows);
                                 v0.1.0 (E == 0 rule, AUDIT B10) as grey reference;
                                 series are dodged by 0.12 % on x so equal ones
                                 do not hide each other
  fig4_tsd_vs_accuracy           TSD exact-length and consistent-repeat accuracy
  fig5_runtime                   per-step and total wall time (INDICATIVE: shared
                                 machine, see run_grid.py manifest for threads)

x axis: nominal mean read accuracy of the ONT model (accuracy sweep, r9.4.1-like
composition). Cells of the chemistry_r10 sweep (R10.4.1 preset) are drawn as
stars at the preset's nominal accuracy. Colours: v0.1.0 neutral grey, dev blue;
E-values a one-hue blue ramp (darker = stricter). Every series also has its own
marker/line style, so identity never depends on colour alone.
"""

import argparse
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
C_BASE, C_DEV = "#52514e", "#2a78d6"
EVALUE_RAMP = ["#86b6ef", "#2a78d6", "#1c5cab", "#0d366b"]   # lenient -> strict
EVALUE_MARKERS = ["o", "s", "D", "^", "v", "P"]
STEP_ORDER = ["filter", "map", "insert_finder_batch", "blast_annot_batch", "dr_finder"]

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "lines.linewidth": 2, "lines.markersize": 5, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 200, "figure.facecolor": "white", "axes.facecolor": "white",
})


def load(path):
    df = pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=False)
    for c in ("pipe.evalue",):
        if c in df:
            df[c] = df[c].astype(str)
    df["sweeps"] = df["sweeps"].fillna("").astype(str)
    return df


def in_sweep(df, name):
    return df[df["sweeps"].str.split(",").apply(lambda s: name in s)]


def parse_default(s):
    out = {}
    for item in s.split(","):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def same_value(a, b):
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-9)
    except (TypeError, ValueError):
        return str(a) == str(b)


def is_dev_default(row, default):
    if row["pipe.name"] == "v0.1.0":
        return False
    return all(same_value(row.get(f"pipe.{k}"), v) for k, v in default.items())


def series_frames(df, dev_default):
    """[(label, colour, marker, linestyle, frame)] for v0.1.0 and the dev default."""
    out = []
    base = df[df["pipe.name"] == "v0.1.0"]
    if len(base):
        out.append(("v0.1.0", C_BASE, "s", "--", base))
    dev = df[df.apply(lambda r: is_dev_default(r, dev_default), axis=1)] if len(df) else df
    if len(dev):
        lab = "dev (" + ", ".join(f"{k} {v}" for k, v in dev_default.items()) + ")"
        out.append((lab, C_DEV, "o", "-", dev))
    return out


def stats_by_x(frame, metric, x="sim.accuracy_nominal"):
    if metric not in frame or frame[metric].notna().sum() == 0:
        return None
    g = frame.groupby(x)[metric].agg(["mean", "std", "count"]).reset_index().sort_values(x)
    return g


def draw_series(ax, frame, metric, colour, marker, ls, label, scale=1.0, r10=None, dx=0.0):
    g = stats_by_x(frame, metric)
    if g is None:
        return False
    xs = g["sim.accuracy_nominal"] * 100 + dx  # dx: small dodge so identical series stay visible
    ys = g["mean"] * scale
    sd = g["std"].fillna(0) * scale
    ax.plot(xs, ys, color=colour, marker=marker, linestyle=ls, label=label, markeredgecolor="white",
            markeredgewidth=0.8, markersize=6, zorder=3)
    if (g["count"] > 1).any():
        ax.fill_between(xs, ys - sd, ys + sd, color=colour, alpha=0.15, linewidth=0, zorder=2)
    if r10 is not None and metric in r10 and r10[metric].notna().any():
        m = r10[metric].mean() * scale
        s = (r10[metric].std() if r10[metric].count() > 1 else 0) * scale
        ax.errorbar(r10["sim.accuracy_nominal"].iloc[0] * 100 + 0.15, m, yerr=s or None, color=colour,
                    marker="*", markersize=11, markeredgecolor="white", markeredgewidth=0.6,
                    linestyle="none", zorder=4, capsize=2)
    return True


def finish_axis(ax, title, ylabel, ylim=None):
    ax.set_title(title, loc="left", color=INK)
    ax.set_xlabel("Mean read accuracy (%)")
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)


def empty_panel(ax, title, why="metric not available"):
    ax.set_title(title, loc="left", color=INK)
    ax.text(0.5, 0.5, why, ha="center", va="center", transform=ax.transAxes, color=MUTED)
    ax.set_xticks([])
    ax.set_yticks([])


def add_r10_legend(fig, has_r10, handles=None):
    from matplotlib.lines import Line2D
    handles = list(handles or [])
    if has_r10:
        handles.append(Line2D([], [], color=INK2, marker="*", linestyle="none", markersize=10,
                              label="R10.4.1 preset (star)"))
    return handles


def save(fig, out_dir, name, note=None):
    if note:
        fig.text(0.01, 0.005, note, fontsize=7, color=MUTED, ha="left", va="bottom")
    paths = []
    for ext in ("pdf", "png"):
        p = os.path.join(out_dir, f"{name}.{ext}")
        fig.savefig(p, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


def accuracy_figure(df, r10_all, dev_default, panels, name, out_dir, note):
    acc = in_sweep(df, "accuracy")
    series = series_frames(acc, dev_default)
    r10_series = dict((lab, fr) for lab, *_x, fr in series_frames(r10_all, dev_default)) if len(r10_all) else {}
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels), 3.2), squeeze=False)
    any_r10 = False
    for ax, (metric, title, ylabel, scale, ylim) in zip(axes[0], panels):
        drawn = False
        for lab, col, mk, ls, fr in series:
            r10 = r10_series.get(lab)
            any_r10 |= r10 is not None and metric in r10 and r10[metric].notna().any()
            drawn |= draw_series(ax, fr, metric, col, mk, ls, lab, scale, r10)
        if drawn:
            finish_axis(ax, title, ylabel, ylim)
        else:
            empty_panel(ax, title)
    h, labels = axes[0][0].get_legend_handles_labels()
    if h:
        fig.legend(handles=add_r10_legend(fig, any_r10, h), loc="upper center", ncol=4, frameon=False,
                   bbox_to_anchor=(0.5, 1.08))
    fig.tight_layout()
    return save(fig, out_dir, name, note)


def classification_figure(df, out_dir, note):
    acc = in_sweep(df, "accuracy")
    dev = acc[acc["pipe.name"] != "v0.1.0"]
    base = acc[acc["pipe.name"] == "v0.1.0"]
    tasks = sorted(dev["pipe.blast_task"].dropna().unique(), key=lambda t: (t != "megablast", t)) if "pipe.blast_task" in dev else []
    tasks = tasks or ["default"]
    evals = sorted(dev["pipe.evalue"].dropna().unique(), key=lambda v: -float(v)) if len(dev) else []
    if len(evals) > len(EVALUE_RAMP):
        sys.stderr.write("plot_grid: more E-values than ramp steps; extra ones share the darkest step\n")
    panels = [("classification.is_insertion.accuracy", "IS-name accuracy (annotated IS reads)", "Accuracy", 1.0, (0, 1.02)),
              ("classification.novel_insertion.accuracy", "Novel inserts labelled Genome", "Accuracy", 1.0, (0, 1.02)),
              ("derived.true_ins.pct_Contamination", "True insertions → Contamination", "% of called reads", 1.0, None),
              ("derived.true_ins.pct_Unannotated", "True insertions → Unannotated", "% of called reads", 1.0, None)]
    fig, axes = plt.subplots(len(tasks), len(panels), figsize=(3.3 * len(panels), 2.9 * len(tasks)),
                             squeeze=False)
    for i, task in enumerate(tasks):
        sub = dev[dev["pipe.blast_task"] == task] if task != "default" else dev
        for j, (metric, title, ylabel, scale, ylim) in enumerate(panels):
            ax = axes[i][j]
            dodge = 0.12
            drawn = draw_series(ax, base, metric, C_BASE, "s", "--", "v0.1.0 (E == 0 rule)", scale,
                                dx=-dodge * (len(evals) + 1) / 2)
            for k, ev in enumerate(evals):
                fr = sub[sub["pipe.evalue"] == ev]
                drawn |= draw_series(ax, fr, metric, EVALUE_RAMP[min(k, len(EVALUE_RAMP) - 1)],
                                     EVALUE_MARKERS[k % len(EVALUE_MARKERS)], "-", f"dev, E ≤ {float(ev):g}", scale,
                                     dx=dodge * (k - (len(evals) - 1) / 2))
            t = f"{title}\n{task}" if len(tasks) > 1 or task != "default" else title
            if drawn:
                finish_axis(ax, t, ylabel, ylim)
                if metric.startswith("derived"):
                    ax.set_ylim(bottom=0)
            else:
                empty_panel(ax, t)
    h, labels = axes[0][0].get_legend_handles_labels()
    if h:
        fig.legend(h, labels, loc="upper center", ncol=len(h), frameon=False, bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout()
    return save(fig, out_dir, "fig3_classification_vs_evalue", note)


def runtime_figure(df, dev_default, out_dir):
    acc = in_sweep(df, "accuracy")
    series = series_frames(acc, dev_default)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.3), gridspec_kw={"width_ratios": [1.4, 1]})
    steps = [s for s in STEP_ORDER if f"runtime.{s}_s" in df]
    width = 0.8 / max(1, len(series))
    for k, (lab, col, mk, ls, fr) in enumerate(series):
        means = [fr[f"runtime.{s}_s"].mean() for s in steps]
        sds = [fr[f"runtime.{s}_s"].std() if fr[f"runtime.{s}_s"].count() > 1 else 0 for s in steps]
        xs = [i + (k - (len(series) - 1) / 2) * width for i in range(len(steps))]
        a1.bar(xs, means, width=width * 0.92, color=col, yerr=sds, capsize=2, label=lab,
               error_kw={"ecolor": INK2, "elinewidth": 0.8}, hatch="//" if ls == "--" else None,
               edgecolor="white", linewidth=0.5)
        draw_series(a2, fr, "runtime.total_s", col, mk, ls, lab)
    a1.set_xticks(range(len(steps)))
    a1.set_xticklabels([s.replace("_batch", "").replace("_", " ") for s in steps], rotation=20, ha="right")
    a1.set_ylabel("Wall time (s)")
    a1.set_title("Per step (mean ± SD over accuracy-sweep cells)", loc="left")
    a1.grid(axis="x", visible=False)
    finish_axis(a2, "Total vs read accuracy", "Wall time (s)")
    a2.set_ylim(bottom=0)
    h, labels = a1.get_legend_handles_labels()
    if h:
        fig.legend(h, labels, loc="upper center", ncol=len(h), frameon=False, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    threads = df["threads"].dropna().unique() if "threads" in df else []
    note = ("Runtimes are INDICATIVE: measured on a shared workstation with parallel cells "
            f"(threads per run: {', '.join(str(int(t)) for t in threads) or '?'}); time on a quiet machine for the paper.")
    return save(fig, out_dir, "fig5_runtime", note)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="grid_results.tsv from run_grid.py")
    ap.add_argument("--out_dir", help="Figure directory (default: figures/ next to --results)")
    ap.add_argument("--dev_default", default="evalue=1e-10,blast_task=megablast",
                    help="Settings of the dev series in figures 1, 2, 4, 5")
    args = ap.parse_args(argv)
    out_dir = args.out_dir or os.path.join(os.path.dirname(os.path.abspath(args.results)), "figures")
    os.makedirs(out_dir, exist_ok=True)
    df = load(args.results)
    dev_default = parse_default(args.dev_default)
    r10 = in_sweep(df, "chemistry_r10")
    note = "Mean ± SD over seeds; ONT-like simulation, r9.4.1-like error composition."
    made = []
    made += accuracy_figure(df, r10, dev_default, [
        ("detection_finder.f1", "Insert finder: read-level F1", "F1", 1.0, None),
        ("detection_annot.f1", "Annotated calls: read-level F1", "F1", 1.0, None),
        ("detection_annot.recall", "Annotated calls: recall", "Recall", 1.0, None)],
        "fig1_detection_vs_accuracy", out_dir, note)
    made += accuracy_figure(df, r10, dev_default, [
        ("position.error.pct_within_0bp", "Insertion position exact", "% of called insertion reads", 1.0, None),
        ("position.error.pct_within_5bp", "Insertion position ≤ 5 bp", "% of called insertion reads", 1.0, None)],
        "fig2_position_vs_accuracy", out_dir, note)
    made += classification_figure(df, out_dir, note + " Annotation settings apply to dev only.")
    made += accuracy_figure(df, r10, dev_default, [
        ("tsd.length_accuracy", "TSD exact length", "Fraction of insertion reads", 1.0, None),
        ("tsd.consistent_accuracy", "TSD consistent with observable repeat", "Fraction of insertion reads",
         1.0, None)],
        "fig4_tsd_vs_accuracy", out_dir, note)
    made += runtime_figure(df, dev_default, out_dir)
    for p in made:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
