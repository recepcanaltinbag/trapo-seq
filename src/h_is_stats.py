#-------------------------------------
#Summarise IS annotations (best_alignment.tab) of every barcode folder into one .rcp file.
#
#.rcp format (parsed by i_heatmap.py):
#   >condition
#   ?<number of IS_DB insertions>            (denominator of the percentages)
#   #total_insertions=<all annotated insertions>
#   #IS_DB=<n>  #Genome=<n>  #Contamination=<n>  (one per line; other Note values too)
#   #percentage_denominator=IS_DB
#   IS_family, IS_name, percentage, count
#   ---
#Lines starting with '#' are metadata and are ignored by the heatmap module.
#
#Counting unit (B41): one insertion = one unique `Query ID`. If the table has several
#rows for the same Query ID, one row is chosen: IS_DB rows win over Genome over
#Contamination, then the highest Score.
#-------------------------------------
import logging
import os

import pandas as pd

logger = logging.getLogger(__name__)

NOTE_ORDER = ["IS_DB", "Genome", "Contamination"]


def split_subject_id(subject_id):
    """`IS5_ISPs1` -> ('IS5', 'ISPs1'); `IS3_ISPa_41` -> ('IS3', 'ISPa_41');
    `ISPs1` (no underscore) -> ('ISPs1', 'ISPs1') (B40)."""
    subject_id = str(subject_id).strip()
    if "_" in subject_id:
        family, name = subject_id.split("_", 1)
        if family and name:
            return family, name
        subject_id = subject_id.strip("_") or subject_id
    return subject_id, subject_id


def dedupe_insertions(df):
    """Return one row per `Query ID` (see module header for the rule)."""
    if df.empty:
        return df
    df = df.copy()
    rank = {note: i for i, note in enumerate(NOTE_ORDER)}
    df["_note_rank"] = df["Note"].map(rank).fillna(len(NOTE_ORDER))
    score = pd.to_numeric(df["Score"], errors="coerce") if "Score" in df.columns else 0
    df["_score"] = score
    df = df.sort_values(["Query ID", "_note_rank", "_score"], ascending=[True, True, False], kind="mergesort")
    df = df.drop_duplicates(subset="Query ID", keep="first")
    return df.drop(columns=["_note_rank", "_score"])


def summarize_condition(tab_file):
    """Summarise one best_alignment.tab.

    Returns a dict with keys: total, note_counts (dict), is_count, rows
    (list of (family, name, percentage, count) sorted by count desc).
    """
    df_raw = pd.read_csv(tab_file, sep='\t')
    missing = {"Query ID", "Subject ID", "Note"} - set(df_raw.columns)
    if missing:
        raise ValueError(f"{tab_file}: missing columns {sorted(missing)}")

    df = dedupe_insertions(df_raw)
    total = int(df.shape[0])
    note_counts = df["Note"].value_counts().to_dict() if total else {}
    note_counts = {str(k): int(v) for k, v in note_counts.items()}

    is_df = df[df["Note"] == "IS_DB"]
    is_count = int(is_df.shape[0])
    rows = []
    if is_count:
        pairs = is_df["Subject ID"].map(split_subject_id)
        counts = pairs.value_counts()
        for (family, name), count in counts.items():
            rows.append((family, name, 100.0 * count / is_count, int(count)))
        rows.sort(key=lambda r: (-r[3], r[0], r[1]))
    return {"total": total, "note_counts": note_counts, "is_count": is_count, "rows": rows}


def format_condition(condition, summary):
    lines = [f">{condition}", f"?{summary['is_count']}",
             f"#total_insertions={summary['total']}"]
    notes = summary["note_counts"]
    for note in NOTE_ORDER + sorted(set(notes) - set(NOTE_ORDER)):
        lines.append(f"#{note}={notes.get(note, 0)}")
    lines.append("#percentage_denominator=IS_DB")
    for family, name, pct, count in summary["rows"]:
        lines.append(f"{family}, {name}, {pct:.2f}, {count}")
    return "\n".join(lines) + "\n"


def main_is_stat(data_dir, out_file):
    folders = sorted(f for f in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, f)))
    segments = []

    for folder in folders:
        tab_file = os.path.join(data_dir, folder, "best_alignment.tab")
        if not os.path.exists(tab_file):
            logger.warning("'%s' not found. Skipping.", tab_file)
            continue

        summary = summarize_condition(tab_file)
        notes = summary["note_counts"]
        breakdown = ", ".join(f"{k}: {notes.get(k, 0)}" for k in NOTE_ORDER)
        print(f"{folder}: {summary['total']} insertions ({breakdown}); "
              f"IS percentages are relative to the {summary['is_count']} IS_DB insertions")
        for family, name, pct, count in summary["rows"]:
            logger.info("  %s %s: %.2f%% (%d)", family, name, pct, count)
        if summary["is_count"] == 0:
            logger.warning("%s: no IS_DB insertions.", folder)
        segments.append(format_condition(folder, summary))

    if not segments:
        logger.warning("No best_alignment.tab found under %s; writing an empty stats file.", data_dir)

    with open(out_file, "w") as f:
        f.write("---\n".join(segments))
    return out_file
