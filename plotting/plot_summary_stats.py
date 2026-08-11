#!/usr/bin/env python3
import argparse
import base64
from html import escape
import gzip
import os
import re
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns


ASSEMBLER_TO_ASSEMBLY_TYPE = {
    "flye": "long_read",
    "metamdbg": "long_read",
    "metaspades_hybrid": "hybrid_read",
    "opera-ms": "hybrid_read",
    "opera_ms": "hybrid_read",
    "metaspades": "short_read",
    "megahit": "short_read",
    "spades": "short_read",
    "idbaud": "short_read",
    "final": "short_read",
    "unknown": "short_read",
}

ASSEMBLY_TYPE_ORDER = ["long_read", "short_read", "hybrid_read"]
DEFAULT_NG50_TARGET_BP = 5_000_000
DEFAULT_TYPE_COLOR_MAP = {
    "short_read": ["#A1D99B", "#74C476", "#238B45"],
    "long_read": ["#C6DBEF", "#41B6C4", "#08519C"],
    "hybrid_read": ["#D0D1E6", "#6A51A3"],
}


ASSEMBLER_LABEL_ALIASES = {
    "flye": "Flye",
    "metamdbg": "metaMDBG",
    "idbaud": "IDBA-UD",
    "megahit": "MEGAHIT",
    "metaspades": "metaSPAdes",
    "metaspades_hybrid": "metaSPAdes Hybrid",
    "opera-ms": "Opera-MS",
    "opera_ms": "Opera-MS",
    "spades": "SPAdes",
    "final": "Final",
}


def assembler_sort_key(assembler: str):
    assembly_type = classify_assembly_type(assembler)
    type_rank = ASSEMBLY_TYPE_ORDER.index(assembly_type) if assembly_type in ASSEMBLY_TYPE_ORDER else len(ASSEMBLY_TYPE_ORDER)
    return (type_rank, str(assembler))


def build_assembler_color_map(assemblers: Sequence[str], color_map_by_type: Optional[Dict[str, Sequence[str]]] = None) -> Dict[str, str]:
    palette_map = color_map_by_type or DEFAULT_TYPE_COLOR_MAP
    assembler_color: Dict[str, str] = {}
    for assembly_type in ASSEMBLY_TYPE_ORDER:
        members = sorted([a for a in assemblers if classify_assembly_type(a) == assembly_type])
        if not members:
            continue
        palette = list(palette_map.get(assembly_type, []))
        if not palette:
            continue
        for idx, assembler in enumerate(members):
            assembler_color[assembler] = palette[idx % len(palette)]
    return assembler_color


def sample_sort_key(sample: object):
    text = str(sample)
    m = re.match(r"^S(\d+)(?:_subsample_(\d+))?$", text)
    if m:
        sample_num = int(m.group(1))
        subsample_num = int(m.group(2)) if m.group(2) else 0
        return (0, sample_num, subsample_num, text)
    return (1, text)


def parse_size_to_mb(text: str) -> Optional[float]:
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(Mbp|Mb|Kbp|Kb|bp)", text)
    if not match:
        return None
    value = float(match.group(1))
    unit = match.group(2)
    if unit in {"Mbp", "Mb"}:
        return value
    if unit in {"Kbp", "Kb"}:
        return value / 1000.0
    if unit == "bp":
        return value / 1_000_000.0
    return None


def parse_int(text: str) -> Optional[int]:
    match = re.search(r"([0-9][0-9,]*)", text)
    if not match:
        return None
    return int(match.group(1).replace(",", ""))


def parse_int_after_colon(text: str) -> Optional[int]:
    if ":" in text:
        text = text.split(":", 1)[1]
    return parse_int(text)


def parse_bp_threshold_to_kbp(label: str) -> Optional[float]:
    cleaned = " ".join(label.strip().split())
    if cleaned.lower() == "all":
        return 0.0

    m = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*(Kbp|Kb|Mbp|Mb|bp)?$", cleaned)
    if not m:
        return None

    value = float(m.group(1))
    unit = m.group(2)
    if unit in {"Kbp", "Kb"}:
        return value
    if unit in {"Mbp", "Mb"}:
        return value * 1000.0
    if unit == "bp":
        return value / 1000.0
    # No unit: bare integers like 50, 100, 250, 500 in the log are base-pair
    # minimum lengths, not kilobases.
    return value / 1000.0


def parse_assembly_meta(path: str) -> Dict[str, Optional[str]]:
    m = re.search(r"/assemblies/(S\d+)/assembly\.([^/]+)/", path)
    if m:
        return {"sample": m.group(1), "assembler": m.group(2)}

    sample_match = re.search(r"/assemblies/(S\d+)/", path)
    sample = sample_match.group(1) if sample_match else None

    filename = os.path.basename(path)
    assembler_map = {
        "megahit.final.contigs.fa": "megahit",
        "final.contigs.fa": "megahit",
        "metaspades.contigs.fasta": "metaspades",
        "metaspades_hybrid.contigs.fasta": "metaspades_hybrid",
        "metaspades_hybrid.assembly.fasta": "metaspades_hybrid",
        "flye.assembly.fasta": "flye",
        "metamdbg.contigs.fasta": "metamdbg",
        "idbaud.assembly.fasta": "idbaud",
    }
    assembler = assembler_map.get(filename)
    if assembler is None and "." in filename:
        assembler = filename.split(".", 1)[0]

    return {"sample": sample, "assembler": assembler}


def classify_assembly_type(assembler: Optional[str]) -> Optional[str]:
    if assembler is None:
        return None
    return ASSEMBLER_TO_ASSEMBLY_TYPE.get(str(assembler).strip().lower(), "short_read")


def _display_assembler_name(assembler: str) -> str:
    return ASSEMBLER_LABEL_ALIASES.get(str(assembler), str(assembler))


def extract_gc(lines: List[str]) -> Optional[float]:
    for idx, line in enumerate(lines):
        if line.strip().startswith("A\tC\tG\tT\tN"):
            if idx + 1 < len(lines):
                parts = lines[idx + 1].strip().split("\t")
                if len(parts) >= 8:
                    try:
                        return float(parts[7])
                    except ValueError:
                        return None
    return None


def parse_block(block: str) -> Optional[Dict[str, object]]:
    lines = [line.rstrip("\n") for line in block.splitlines() if line.strip()]
    if not lines:
        return None

    assembly_line = next((line for line in lines if line.startswith("Assembly:")), None)
    if not assembly_line:
        return None

    sample_line = next((line for line in lines if line.startswith("Sample:")), None)
    assembler_line = next((line for line in lines if line.startswith("Assembler:")), None)

    assembly_path = assembly_line.split("Assembly:", 1)[1].strip()
    meta = parse_assembly_meta(assembly_path)
    sample = sample_line.split("Sample:", 1)[1].strip() if sample_line else meta["sample"]
    assembler = assembler_line.split("Assembler:", 1)[1].strip() if assembler_line else meta["assembler"]
    assembly_type = classify_assembly_type(assembler)

    values: Dict[str, object] = {
        "assembly_path": assembly_path,
        "sample": sample,
        "assembler": assembler,
        "assembly_type": assembly_type,
        "gc_fraction": extract_gc(lines),
        "scaffold_total": None,
        "sequence_total_mb": None,
        "n50_kbp": None,
        "l50": None,
        "max_scaffold_kbp": None,
        "max_contig_kbp": None,
        "scaffolds_gt_50kb": None,
    }

    for line in lines:
        if line.startswith("Main genome scaffold total:"):
            values["scaffold_total"] = parse_int(line)
        elif line.startswith("Main genome scaffold sequence total:"):
            values["sequence_total_mb"] = parse_size_to_mb(line)
        elif line.startswith("Main genome scaffold N/L50:"):
            m = re.search(
                r"([0-9][0-9,]*)\s*/\s*([0-9][0-9,]*(?:\.[0-9]+)?)\s*(Kbp|Kb|Mbp|Mb|bp)?",
                line,
            )
            if m:
                values["l50"] = int(m.group(1).replace(",", ""))
                n50_val = float(m.group(2).replace(",", ""))
                n50_unit = m.group(3) or "bp"
                if n50_unit in {"Kbp", "Kb"}:
                    values["n50_kbp"] = n50_val
                elif n50_unit in {"Mbp", "Mb"}:
                    values["n50_kbp"] = n50_val * 1000.0
                elif n50_unit == "bp":
                    values["n50_kbp"] = n50_val / 1000.0
        elif line.startswith("Max scaffold length:"):
            max_mb = parse_size_to_mb(line)
            if max_mb is not None:
                values["max_scaffold_kbp"] = max_mb * 1000.0
        elif line.startswith("Max contig length:"):
            max_contig_mb = parse_size_to_mb(line)
            if max_contig_mb is not None:
                values["max_contig_kbp"] = max_contig_mb * 1000.0
        elif line.startswith("Number of scaffolds > 50 KB:"):
            values["scaffolds_gt_50kb"] = parse_int_after_colon(line)

    return values


def parse_threshold_rows(block: str) -> List[Dict[str, object]]:
    lines = [line.rstrip("\n") for line in block.splitlines()]
    sample_line = next((line for line in lines if line.strip().startswith("Sample:")), None)
    assembler_line = next((line for line in lines if line.strip().startswith("Assembler:")), None)
    assembly_line = next((line for line in lines if line.strip().startswith("Assembly:")), None)

    if assembly_line is None:
        return []

    assembly_path = assembly_line.split("Assembly:", 1)[1].strip()
    meta = parse_assembly_meta(assembly_path)
    sample = sample_line.split("Sample:", 1)[1].strip() if sample_line else meta["sample"]
    assembler = assembler_line.split("Assembler:", 1)[1].strip() if assembler_line else meta["assembler"]
    assembly_type = classify_assembly_type(assembler)

    header_sep_idx = None
    for idx, line in enumerate(lines):
        if line.strip().startswith("--------"):
            header_sep_idx = idx
            break

    if header_sep_idx is None:
        return []

    table_rows: List[Dict[str, object]] = []
    row_pattern = re.compile(
        r"^\s*(All|[0-9]+(?:\.[0-9]+)?(?:\s*(?:Kbp|Kb|Mbp|Mb|bp))?)\s+"
        r"([0-9][0-9,]*)\s+"
        r"([0-9][0-9,]*)\s+"
        r"([0-9][0-9,]*)\s+"
        r"([0-9][0-9,]*)\s+"
        r"([0-9]+(?:\.[0-9]+)?)%\s*$"
    )

    for line in lines[header_sep_idx + 1 :]:
        if not line.strip():
            break
        m = row_pattern.match(" ".join(line.strip().split()))
        if not m:
            continue

        threshold_label = m.group(1)
        threshold_kbp = parse_bp_threshold_to_kbp(threshold_label)
        table_rows.append(
            {
                "assembly_path": assembly_path,
                "sample": sample,
                "assembler": assembler,
                "assembly_type": assembly_type,
                "threshold_label": threshold_label,
                "threshold_kbp": threshold_kbp,
                "num_scaffolds": int(m.group(2).replace(",", "")),
                "num_contigs": int(m.group(3).replace(",", "")),
                "scaffold_length_bp": int(m.group(4).replace(",", "")),
                "contig_length_bp": int(m.group(5).replace(",", "")),
                "contig_coverage_pct": float(m.group(6)),
            }
        )

    return table_rows


def split_summary_blocks(log_path: str) -> List[str]:
    with open(log_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    blocks: List[str] = []
    current_block: List[str] = []

    for line in lines:
        if line.startswith("Timestamp:"):
            if current_block:
                blocks.append("".join(current_block))
                current_block = []
            current_block.append(line)
            continue

        if line.startswith("Assembly:") and not current_block:
            current_block.append(line)
            continue

        if current_block:
            current_block.append(line)

    if current_block:
        blocks.append("".join(current_block))

    return blocks


def parse_summary_stats_log(log_path: str) -> pd.DataFrame:
    blocks = split_summary_blocks(log_path)

    parsed = [parse_block(block) for block in blocks]
    parsed = [row for row in parsed if row is not None]
    return pd.DataFrame(parsed)


def parse_threshold_stats_log(log_path: str) -> pd.DataFrame:
    blocks = split_summary_blocks(log_path)
    rows: List[Dict[str, object]] = []
    for block in blocks:
        rows.extend(parse_threshold_rows(block))
    return pd.DataFrame(rows)


def add_percent_of_best(
    df: pd.DataFrame,
    value_col: str,
    group_cols: List[str],
    higher_is_better: bool,
    out_col: str,
) -> pd.DataFrame:
    if df.empty or value_col not in df.columns:
        return df

    result = df.copy()

    if higher_is_better:
        best = result.groupby(group_cols)[value_col].transform("max")
        result[out_col] = (result[value_col] / best) * 100.0
    else:
        best = result.groupby(group_cols)[value_col].transform("min")
        result[out_col] = (best / result[value_col]) * 100.0

    result.loc[best <= 0, out_col] = pd.NA
    result.loc[result[value_col] <= 0, out_col] = pd.NA
    return result


def _open_text_maybe_gzip(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _fasta_lengths(path: str) -> List[int]:
    lengths: List[int] = []
    current = 0

    with _open_text_maybe_gzip(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                current = 0
                continue
            current += len(line)

    if current > 0:
        lengths.append(current)
    return lengths


def _nx(lengths: Sequence[int], cutoff_bp: int) -> Tuple[Optional[int], Optional[int]]:
    if not lengths or cutoff_bp <= 0:
        return None, None

    sorted_lengths = sorted(lengths, reverse=True)
    cumulative = 0
    for idx, length in enumerate(sorted_lengths, start=1):
        cumulative += length
        if cumulative >= cutoff_bp:
            return length, idx
    return None, None


def add_recalculated_ng50(df: pd.DataFrame, target_genome_bp: int = DEFAULT_NG50_TARGET_BP) -> pd.DataFrame:
    result = df.copy()
    if result.empty or "assembly_path" not in result.columns or target_genome_bp <= 0:
        return result

    cache: Dict[str, Tuple[Optional[float], Optional[int], Optional[float], Optional[int]]] = {}
    missing_paths = set()

    def _compute_for_path(path: object) -> Tuple[Optional[float], Optional[int], Optional[float], Optional[int]]:
        if pd.isna(path):
            return None, None, None, None
        apath = os.path.abspath(str(path))
        if apath in cache:
            return cache[apath]
        if not os.path.isfile(apath):
            missing_paths.add(apath)
            cache[apath] = (None, None, None, None)
            return cache[apath]

        lengths = _fasta_lengths(apath)
        total_bp = sum(lengths)
        true_n50_bp, true_l50 = _nx(lengths, (total_bp + 1) // 2)
        ng50_bp, lg50 = _nx(lengths, (target_genome_bp + 1) // 2)

        cache[apath] = (
            (true_n50_bp / 1000.0) if true_n50_bp is not None else None,
            true_l50,
            (ng50_bp / 1000.0) if ng50_bp is not None else None,
            lg50,
        )
        return cache[apath]

    computed = result["assembly_path"].apply(_compute_for_path)
    computed_df = pd.DataFrame(
        computed.tolist(),
        columns=["true_n50_kbp", "true_l50", "ng50_kbp", "lg50"],
        index=result.index,
    )
    for col in computed_df.columns:
        result[col] = computed_df[col]
    result["ng50_target_bp"] = int(target_genome_bp)

    if missing_paths:
        print(f"[WARN] NG50 skipped for {len(missing_paths)} missing assembly file(s)")

    return result


def _contig_lengths_from_fasta(path: str) -> List[int]:
    lengths: List[int] = []
    current = 0

    with _open_text_maybe_gzip(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                current = 0
                continue
            current += len(line)

    if current > 0:
        lengths.append(current)
    return lengths


def _combined_subsample_curve_specs() -> List[Dict[str, str]]:
    return [
        {"label": "S1", "sample": "S1", "log_path": "summary_stats_log.txt"},
        *[
            {
                "label": f"Subsample {idx}",
                "sample": f"S1_subsample_{idx}",
                "log_path": "subsample_assembly_summary_stats_log.txt",
            }
            for idx in range(1, 12)
        ],
    ]


def _parse_contig_curve_log_rows(log_path: str, sample_names: Sequence[str]) -> Dict[Tuple[str, str], Dict[str, object]]:
    if not os.path.isfile(log_path):
        return {}

    wanted = set(sample_names)
    rows: Dict[Tuple[str, str], Dict[str, object]] = {}

    for block in split_summary_blocks(log_path):
        parsed = parse_block(block)
        if not parsed:
            continue
        sample = str(parsed.get("sample") or "")
        if sample not in wanted:
            continue
        assembly_path = parsed.get("assembly_path")
        if not assembly_path or not os.path.isfile(str(assembly_path)):
            continue

        lengths = sorted(_contig_lengths_from_fasta(str(assembly_path)), reverse=True)
        if not lengths:
            continue

        cumulative_bp = []
        running = 0
        for length in lengths:
            running += length
            cumulative_bp.append(running)

        assembler = str(parsed.get("assembler") or "unknown")
        rows[(sample, assembler)] = {
            "sample": sample,
            "assembler": assembler,
            "assembly_type": classify_assembly_type(assembler),
            "cumulative_mbp": [bp / 1_000_000.0 for bp in cumulative_bp],
            "contig_count": list(range(1, len(lengths) + 1)),
            "total_mbp": running / 1_000_000.0,
        }

    return rows


def make_subsample_cumulative_contig_plot(outdir: str, reference_log: str = "summary_stats_log.txt", subsample_log: str = "subsample_assembly_summary_stats_log.txt") -> Optional[str]:
    os.makedirs(outdir, exist_ok=True)

    specs = _combined_subsample_curve_specs()
    sample_names = [spec["sample"] for spec in specs]
    ref_rows = _parse_contig_curve_log_rows(reference_log, ["S1"])
    sub_rows = _parse_contig_curve_log_rows(subsample_log, [name for name in sample_names if name != "S1"])

    rows_by_key = {**ref_rows, **sub_rows}
    if not rows_by_key:
        print("[WARN] No cumulative-length data found for combined subsample plot")
        return None

    assemblers = [row["assembler"] for row in rows_by_key.values() if row.get("assembler")]
    assembler_order = sorted(dict.fromkeys(assemblers), key=assembler_sort_key)
    if not assembler_order:
        print("[WARN] No assemblers found for combined subsample plot")
        return None

    colors = build_assembler_color_map(assembler_order)
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharex=False, sharey=False)
    axes = axes.flatten()

    pane_order = ["flye", "metamdbg", "metaspades_hybrid", "opera_ms"]
    for ax, assembler in zip(axes, pane_order):
        pane_rows = []
        for spec in specs:
            key = (spec["sample"], assembler)
            row = rows_by_key.get(key)
            if row is None:
                continue
            pane_rows.append({"label": spec["label"], **row})
        if not pane_rows:
            ax.set_title(f"{assembler} (no data)")
            ax.axis("off")
            continue

        for row in pane_rows:
            label = row["label"]
            x_vals = row["contig_count"]
            y_vals = row["cumulative_mbp"]
            color = colors.get(row["assembler"], "#4c72b0")
            ax.plot(x_vals, y_vals, linewidth=2, color=color, label=label)

        ax.set_title(_display_assembler_name(assembler))
        ax.set_xlabel("Contig count (longest to shortest)")
        ax.set_ylabel("Cumulative contig length (Mb)")
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.legend(title="Sample", fontsize=7, title_fontsize=8)

    fig.suptitle("S1 True Assembly + Subsample 1-11 Cumulative Contig Length", y=1.02)
    fig.tight_layout()
    outpath = os.path.join(outdir, "subsample_s1_cumulative_contig_length_four_panel.png")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def make_plots(df: pd.DataFrame, threshold_df: pd.DataFrame, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    sns.set_context("notebook")
    sns.set_style("darkgrid")
    sns.set_palette("muted")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "legend.title_fontsize": 11,
        }
    )

    combined_plot = make_subsample_cumulative_contig_plot(outdir)
    if combined_plot:
        print(f"Saved combined subsample plot: {combined_plot}")

    plot_df = df.dropna(subset=["sample", "assembler"]).copy()
    sample_order = sorted(plot_df["sample"].unique(), key=sample_sort_key)
    assembler_order = sorted(plot_df["assembler"].unique(), key=assembler_sort_key)
    plot_df["sample"] = pd.Categorical(plot_df["sample"], categories=sample_order, ordered=True)
    plot_df["assembler"] = pd.Categorical(plot_df["assembler"], categories=assembler_order, ordered=True)
    palette = build_assembler_color_map(assembler_order)

    for metric, higher_is_better in [
        ("sequence_total_mb", True),
        ("n50_kbp", True),
        ("ng50_kbp", True),
        ("l50", False),
        ("max_scaffold_kbp", True),
        ("max_contig_kbp", True),
        ("scaffold_total", False),
        ("scaffolds_gt_50kb", True),
    ]:
        if metric in plot_df.columns:
            plot_df = add_percent_of_best(
                plot_df,
                value_col=metric,
                group_cols=["sample"],
                higher_is_better=higher_is_better,
                out_col=f"{metric}_pct_best",
            )

    # 1) Total assembly size
    size_df = plot_df.dropna(subset=["sequence_total_mb"])
    if not size_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=size_df, x="sample", y="sequence_total_mb", hue="assembler", palette=palette, errorbar=None)
        plt.title("Assembly Size by Sample and Assembler")
        plt.ylabel("Total Scaffold Sequence (Mb)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "assembly_size_mb.png"), dpi=300)
        plt.close()

    # 2) N50
    n50_df = plot_df.dropna(subset=["n50_kbp"])
    if not n50_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=n50_df, x="sample", y="n50_kbp", hue="assembler", palette=palette, errorbar=None)
        plt.title("Scaffold N50 by Sample and Assembler")
        plt.ylabel("N50 (Kbp)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "assembly_n50_kbp.png"), dpi=300)
        plt.close()

    # 2a) NG50
    ng50_df = plot_df.dropna(subset=["ng50_kbp"])
    if not ng50_df.empty:
        target_bp = int(ng50_df["ng50_target_bp"].dropna().iloc[0]) if "ng50_target_bp" in ng50_df.columns else DEFAULT_NG50_TARGET_BP
        plt.figure(figsize=(14, 7))
        sns.barplot(data=ng50_df, x="sample", y="ng50_kbp", hue="assembler", palette=palette, errorbar=None)
        plt.title(f"Scaffold NG50 by Sample and Assembler (Target: {target_bp:,} bp)")
        plt.ylabel("NG50 (Kbp)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "assembly_ng50_kbp.png"), dpi=300)
        plt.close()

    # 2b) L50
    l50_df = plot_df.dropna(subset=["l50"])
    if not l50_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=l50_df, x="sample", y="l50", hue="assembler", palette=palette, errorbar=None)
        plt.title("Scaffold L50 by Sample and Assembler")
        plt.ylabel("L50")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "assembly_l50.png"), dpi=300)
        plt.close()

    # 2c) Main genome scaffold total
    scaffold_total_df = plot_df.dropna(subset=["scaffold_total"])
    if not scaffold_total_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=scaffold_total_df, x="sample", y="scaffold_total", hue="assembler", palette=palette, errorbar=None)
        plt.title("Main Genome Scaffold Total by Sample and Assembler")
        plt.ylabel("Scaffold Total")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "main_genome_scaffold_total.png"), dpi=300)
        plt.close()

    # 3) Max scaffold length
    max_df = plot_df.dropna(subset=["max_scaffold_kbp"])
    if not max_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=max_df, x="sample", y="max_scaffold_kbp", hue="assembler", palette=palette, errorbar=None)
        plt.title("Max Scaffold Length by Sample and Assembler")
        plt.ylabel("Max Scaffold Length (Kbp)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "max_scaffold_kbp.png"), dpi=300)
        plt.close()

    # 3b) Max contig length
    max_contig_df = plot_df.dropna(subset=["max_contig_kbp"])
    if not max_contig_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=max_contig_df, x="sample", y="max_contig_kbp", hue="assembler", palette=palette, errorbar=None)
        plt.title("Max Contig Length by Sample and Assembler")
        plt.ylabel("Max Contig Length (Kbp)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "max_contig_kbp.png"), dpi=300)
        plt.close()

    # 4) GC vs size scatter
    plt.figure(figsize=(10, 8))
    gc_df = plot_df.dropna(subset=["gc_fraction", "sequence_total_mb"])
    if not gc_df.empty:
        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            data=gc_df,
            x="gc_fraction",
            y="sequence_total_mb",
            hue="assembler",
            palette=palette,
            style="sample",
            s=140,
        )
        plt.title("GC Fraction vs Assembly Size")
        plt.xlabel("GC Fraction")
        plt.ylabel("Total Scaffold Sequence (Mb)")
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "gc_vs_size.png"), dpi=300)
        plt.close()

    # 5) Within-sample comparison (size)
    if not size_df.empty:
        g = sns.catplot(
            data=size_df,
            x="assembler",
            y="sequence_total_mb",
            col="sample",
            col_wrap=4,
            kind="bar",
            order=assembler_order,
            sharey=False,
            height=3.4,
            aspect=1.1,
            errorbar=None,
            palette=palette,
        )
        g.set_axis_labels("Assembler", "Total Scaffold Sequence (Mb)")
        g.set_titles("{col_name}")
        for ax in g.axes.flatten():
            ax.tick_params(axis="x", rotation=35)
        g.fig.suptitle("Within-Sample Assembler Comparison: Assembly Size", y=1.03)
        g.fig.tight_layout()
        g.fig.savefig(os.path.join(outdir, "within_sample_size_facets.png"), dpi=300)
        plt.close(g.fig)

    # 6) Within-sample comparison (N50)
    if not n50_df.empty:
        g = sns.catplot(
            data=n50_df,
            x="assembler",
            y="n50_kbp",
            col="sample",
            col_wrap=4,
            kind="bar",
            order=assembler_order,
            sharey=False,
            height=3.4,
            aspect=1.1,
            errorbar=None,
            palette=palette,
        )
        g.set_axis_labels("Assembler", "N50 (Kbp)")
        g.set_titles("{col_name}")
        for ax in g.axes.flatten():
            ax.tick_params(axis="x", rotation=35)
        g.fig.suptitle("Within-Sample Assembler Comparison: N50", y=1.03)
        g.fig.tight_layout()
        g.fig.savefig(os.path.join(outdir, "within_sample_n50_facets.png"), dpi=300)
        plt.close(g.fig)

    # 6b) Within-sample comparison (NG50)
    if not ng50_df.empty:
        g = sns.catplot(
            data=ng50_df,
            x="assembler",
            y="ng50_kbp",
            col="sample",
            col_wrap=4,
            kind="bar",
            order=assembler_order,
            sharey=False,
            height=3.4,
            aspect=1.1,
            errorbar=None,
            palette=palette,
        )
        g.set_axis_labels("Assembler", "NG50 (Kbp)")
        g.set_titles("{col_name}")
        for ax in g.axes.flatten():
            ax.tick_params(axis="x", rotation=35)
        g.fig.suptitle("Within-Sample Assembler Comparison: NG50", y=1.03)
        g.fig.tight_layout()
        g.fig.savefig(os.path.join(outdir, "within_sample_ng50_facets.png"), dpi=300)
        plt.close(g.fig)

    # 7) Heatmap view for quick within-sample comparisons
    size_matrix = size_df.pivot_table(index="assembler", columns="sample", values="sequence_total_mb", aggfunc="first")
    if not size_matrix.empty and size_matrix.notna().any().any():
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(size_matrix, annot=True, fmt=".1f", cmap="viridis")
        plt.title("Assembly Size (Mb): Assembler vs Sample")
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "within_sample_size_heatmap.png"), dpi=300)
        plt.close()

    n50_matrix = n50_df.pivot_table(index="assembler", columns="sample", values="n50_kbp", aggfunc="first")
    if not n50_matrix.empty and n50_matrix.notna().any().any():
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(n50_matrix, annot=True, fmt=".1f", cmap="magma")
        plt.title("N50 (Kbp): Assembler vs Sample")
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "within_sample_n50_heatmap.png"), dpi=300)
        plt.close()

    ng50_matrix = ng50_df.pivot_table(index="assembler", columns="sample", values="ng50_kbp", aggfunc="first")
    if not ng50_matrix.empty and ng50_matrix.notna().any().any():
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(ng50_matrix, annot=True, fmt=".1f", cmap="rocket")
        plt.title("NG50 (Kbp): Assembler vs Sample")
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "within_sample_ng50_heatmap.png"), dpi=300)
        plt.close()

    l50_matrix = l50_df.pivot_table(index="assembler", columns="sample", values="l50", aggfunc="first")
    if not l50_matrix.empty and l50_matrix.notna().any().any():
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(l50_matrix, annot=True, fmt=".0f", cmap="cividis")
        plt.title("L50: Assembler vs Sample")
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "within_sample_l50_heatmap.png"), dpi=300)
        plt.close()

    for pct_col, title, outfile in [
        ("n50_kbp_pct_best", "N50 (% of Best Assembler in Sample)", "n50_pct_best_heatmap.png"),
        ("ng50_kbp_pct_best", "NG50 (% of Best Assembler in Sample)", "ng50_pct_best_heatmap.png"),
        (
            "scaffold_total_pct_best",
            "Main Genome Scaffold Total (% of Best; lower is better)",
            "main_genome_scaffold_total_pct_best_heatmap.png",
        ),
        (
            "max_contig_kbp_pct_best",
            "Max Contig Length (% of Best Assembler in Sample)",
            "max_contig_pct_best_heatmap.png",
        ),
    ]:
        if pct_col not in plot_df.columns:
            continue
        pct_matrix = plot_df.pivot_table(index="assembler", columns="sample", values=pct_col, aggfunc="first")
        if pct_matrix.empty or not pct_matrix.notna().any().any():
            continue
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(pct_matrix, annot=True, fmt=".1f", cmap="Blues", vmin=0, vmax=100)
        plt.title(title)
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, outfile), dpi=300)
        plt.close()

    if not threshold_df.empty:
        threshold_plot_df = threshold_df.dropna(subset=["sample", "assembler", "threshold_kbp"]).copy()
        if not threshold_plot_df.empty:
            threshold_plot_df["scaffold_length_mb"] = threshold_plot_df["scaffold_length_bp"] / 1_000_000.0
            threshold_plot_df["threshold_label"] = threshold_plot_df["threshold_label"].astype(str)
            threshold_plot_df = add_percent_of_best(
                threshold_plot_df,
                value_col="num_contigs",
                group_cols=["sample", "threshold_kbp"],
                higher_is_better=True,
                out_col="num_contigs_pct_best",
            )
            threshold_plot_df = add_percent_of_best(
                threshold_plot_df,
                value_col="scaffold_length_mb",
                group_cols=["sample", "threshold_kbp"],
                higher_is_better=True,
                out_col="scaffold_length_mb_pct_best",
            )

            threshold_order = (
                threshold_plot_df[["threshold_label", "threshold_kbp"]]
                .drop_duplicates()
                .sort_values("threshold_kbp")
                ["threshold_label"]
                .tolist()
            )

            threshold_plot_df["sample"] = pd.Categorical(
                threshold_plot_df["sample"], categories=sample_order, ordered=True
            )
            threshold_plot_df["assembler"] = pd.Categorical(
                threshold_plot_df["assembler"], categories=assembler_order, ordered=True
            )
            threshold_plot_df["threshold_label"] = pd.Categorical(
                threshold_plot_df["threshold_label"], categories=threshold_order, ordered=True
            )

            key_thresholds = [250.0, 500.0, 1000.0]
            key_df = threshold_plot_df[threshold_plot_df["threshold_kbp"].isin(key_thresholds)].copy()
            if not key_df.empty:
                g = sns.catplot(
                    data=key_df,
                    x="assembler",
                    y="num_contigs",
                    hue="threshold_label",
                    col="sample",
                    col_wrap=4,
                    kind="bar",
                    order=assembler_order,
                    sharey=False,
                    height=3.5,
                    aspect=1.15,
                    errorbar=None,
                )
                g.set_axis_labels("Assembler", "Number of Contigs")
                g.set_titles("{col_name}")
                for ax in g.axes.flatten():
                    ax.tick_params(axis="x", rotation=35)
                g.fig.suptitle(
                    "Contigs Above Key Thresholds (250/500 Kbp, 1 Mbp): Within-Sample Comparison", y=1.03
                )
                g.fig.tight_layout()
                g.fig.savefig(os.path.join(outdir, "within_sample_contigs_key_thresholds.png"), dpi=300)
                plt.close(g.fig)

                for threshold_kbp, label_slug in [(250.0, "250kbp"), (500.0, "500kbp"), (1000.0, "1mbp")]:
                    threshold_slice = key_df[key_df["threshold_kbp"] == threshold_kbp]
                    if threshold_slice.empty:
                        continue
                    matrix = threshold_slice.pivot_table(
                        index="assembler", columns="sample", values="num_contigs", aggfunc="first"
                    )
                    if matrix.empty or not matrix.notna().any().any():
                        continue
                    plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
                    sns.heatmap(matrix, annot=True, fmt=".0f", cmap="YlGnBu")
                    plt.title(f"Number of Contigs >= {threshold_kbp:.0f} Kbp")
                    plt.xlabel("Sample")
                    plt.ylabel("Assembler")
                    plt.tight_layout()
                    plt.savefig(os.path.join(outdir, f"contigs_ge_{label_slug}_heatmap.png"), dpi=300)
                    plt.close()

                for threshold_kbp, label_slug in [(250.0, "250kbp"), (500.0, "500kbp"), (1000.0, "1mbp")]:
                    threshold_slice = key_df[key_df["threshold_kbp"] == threshold_kbp]
                    if threshold_slice.empty:
                        continue
                    pct_matrix = threshold_slice.pivot_table(
                        index="assembler", columns="sample", values="num_contigs_pct_best", aggfunc="first"
                    )
                    if pct_matrix.empty or not pct_matrix.notna().any().any():
                        continue
                    plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
                    sns.heatmap(pct_matrix, annot=True, fmt=".1f", cmap="Blues", vmin=0, vmax=100)
                    plt.title(f"Contigs >= {threshold_kbp:.0f} Kbp (% of Best in Sample)")
                    plt.xlabel("Sample")
                    plt.ylabel("Assembler")
                    plt.tight_layout()
                    plt.savefig(os.path.join(outdir, f"contigs_ge_{label_slug}_pct_best_heatmap.png"), dpi=300)
                    plt.close()

            g = sns.catplot(
                data=threshold_plot_df,
                x="threshold_label",
                y="scaffold_length_mb",
                hue="assembler",
                col="sample",
                col_wrap=4,
                kind="point",
                sharey=False,
                height=3.6,
                aspect=1.2,
                palette=palette,
            )
            g.set_axis_labels("Minimum Scaffold Length", "Total Scaffold Length (Mb)")
            g.set_titles("{col_name}")
            for ax in g.axes.flatten():
                ax.tick_params(axis="x", rotation=45)
            g.fig.suptitle("Scaffold Length Across Minimum Scaffold Length Thresholds", y=1.03)
            g.fig.tight_layout()
            g.fig.savefig(os.path.join(outdir, "within_sample_scaffold_length_across_thresholds.png"), dpi=300)
            plt.close(g.fig)

            g = sns.catplot(
                data=threshold_plot_df,
                x="threshold_label",
                y="scaffold_length_mb_pct_best",
                hue="assembler",
                col="sample",
                col_wrap=4,
                kind="point",
                sharey=True,
                height=3.6,
                aspect=1.2,
                palette=palette,
            )
            g.set_axis_labels("Minimum Scaffold Length", "Scaffold Length (% of Best)")
            g.set_titles("{col_name}")
            for ax in g.axes.flatten():
                ax.tick_params(axis="x", rotation=45)
                ax.set_ylim(0, 105)
            g.fig.suptitle("Scaffold Length Across Thresholds (% of Best in Sample)", y=1.03)
            g.fig.tight_layout()
            g.fig.savefig(os.path.join(outdir, "within_sample_scaffold_length_pct_best_across_thresholds.png"), dpi=300)
            plt.close(g.fig)

    focus_samples = [sample for sample in ["S1", "S2", "S5"] if sample in sample_order]
    focus_df = plot_df[plot_df["sample"].isin(focus_samples)].copy()
    if not focus_df.empty:
        focus_assembler_order = sorted(focus_df["assembler"].dropna().unique())

        pct_metric_map = {
            "sequence_total_mb_pct_best": "Assembly Size\n(% best)",
            "n50_kbp_pct_best": "N50\n(% best)",
            "ng50_kbp_pct_best": "NG50\n(% best)",
            "l50_pct_best": "L50\n(% best, lower=better)",
            "max_contig_kbp_pct_best": "Max Contig\n(% best)",
            "scaffold_total_pct_best": "Scaffold Total\n(% best, lower=better)",
            "scaffolds_gt_50kb_pct_best": ">50 Kb Scaffolds\n(% best)",
        }

        pct_cols = [col for col in pct_metric_map if col in focus_df.columns]
        if pct_cols:
            zoom_pct = focus_df[["sample", "assembler"] + pct_cols].copy()
            zoom_pct["sample_assembler"] = zoom_pct["sample"].astype(str) + " | " + zoom_pct["assembler"].astype(str)
            metric_renamed = zoom_pct.rename(columns=pct_metric_map)
            heatmap_data = metric_renamed.set_index("sample_assembler")[list(pct_metric_map[col] for col in pct_cols)]

            row_order = []
            for sample in focus_samples:
                for assembler in focus_assembler_order:
                    label = f"{sample} | {assembler}"
                    if label in heatmap_data.index:
                        row_order.append(label)
            heatmap_data = heatmap_data.loc[row_order]

            if not heatmap_data.empty:
                plt.figure(figsize=(max(10, heatmap_data.shape[1] * 1.8), max(6, heatmap_data.shape[0] * 0.45)))
                sns.heatmap(heatmap_data, annot=True, fmt=".1f", cmap="Blues", vmin=0, vmax=100)
                plt.title("Zoom: S1/S2/S5 Key Metrics (% of Best Assembler per Sample)")
                plt.xlabel("Metric")
                plt.ylabel("Sample | Assembler")
                plt.tight_layout()
                plt.savefig(os.path.join(outdir, "zoom_s1_s2_s5_key_metrics_pct_best_heatmap.png"), dpi=300)
                plt.close()

        raw_metric_map = {
            "sequence_total_mb": "Assembly Size (Mb)",
            "n50_kbp": "N50 (Kbp)",
            "ng50_kbp": "NG50 (Kbp)",
            "l50": "L50",
            "max_contig_kbp": "Max Contig (Kbp)",
            "scaffold_total": "Main Genome Scaffold Total",
            "scaffolds_gt_50kb": "Scaffolds >50 Kb",
        }
        raw_cols = [col for col in raw_metric_map if col in focus_df.columns]
        if raw_cols:
            long_raw = focus_df[["sample", "assembler"] + raw_cols].melt(
                id_vars=["sample", "assembler"], value_vars=raw_cols, var_name="metric", value_name="value"
            )
            long_raw = long_raw.dropna(subset=["value"])
            long_raw["metric"] = long_raw["metric"].map(raw_metric_map)
            long_raw["sample"] = pd.Categorical(long_raw["sample"], categories=focus_samples, ordered=True)
            long_raw["assembler"] = pd.Categorical(
                long_raw["assembler"], categories=focus_assembler_order, ordered=True
            )

            if not long_raw.empty:
                g = sns.catplot(
                    data=long_raw,
                    x="assembler",
                    y="value",
                    hue="sample",
                    col="metric",
                    col_wrap=3,
                    kind="bar",
                    sharey=False,
                    height=3.8,
                    aspect=1.15,
                    errorbar=None,
                )
                g.set_axis_labels("Assembler", "Value")
                g.set_titles("{col_name}")
                for ax in g.axes.flatten():
                    ax.tick_params(axis="x", rotation=35)
                g.fig.suptitle("Zoom: S1/S2/S5 Key Raw Metrics by Assembler", y=1.03)
                g.fig.tight_layout()
                g.fig.savefig(os.path.join(outdir, "zoom_s1_s2_s5_key_metrics_raw_facets.png"), dpi=300)
                plt.close(g.fig)

    if not threshold_df.empty:
        focus_samples = [sample for sample in ["S1", "S2", "S5"] if sample in sample_order]
        threshold_focus = threshold_df[
            threshold_df["sample"].isin(focus_samples) & threshold_df["threshold_kbp"].isin([250.0, 500.0, 1000.0])
        ].copy()
        if not threshold_focus.empty:
            threshold_focus = add_percent_of_best(
                threshold_focus,
                value_col="num_contigs",
                group_cols=["sample", "threshold_kbp"],
                higher_is_better=True,
                out_col="num_contigs_pct_best",
            )
            threshold_focus["threshold_label"] = threshold_focus["threshold_kbp"].map(
                {250.0: "250 Kbp", 500.0: "500 Kbp", 1000.0: "1 Mbp"}
            )

            assembler_focus_order = sorted(threshold_focus["assembler"].dropna().unique())
            threshold_focus["sample"] = pd.Categorical(
                threshold_focus["sample"], categories=focus_samples, ordered=True
            )
            threshold_focus["assembler"] = pd.Categorical(
                threshold_focus["assembler"], categories=assembler_focus_order, ordered=True
            )

            g = sns.catplot(
                data=threshold_focus,
                x="assembler",
                y="num_contigs_pct_best",
                hue="threshold_label",
                col="sample",
                col_wrap=3,
                kind="bar",
                sharey=True,
                height=4.0,
                aspect=1.15,
                errorbar=None,
            )
            g.set_axis_labels("Assembler", "Contigs (% of Best)")
            g.set_titles("{col_name}")
            for ax in g.axes.flatten():
                ax.tick_params(axis="x", rotation=35)
                ax.set_ylim(0, 105)
            g.fig.suptitle("Zoom: S1/S2/S5 Contigs at 250/500 Kbp and 1 Mbp (% of Best)", y=1.03)
            g.fig.tight_layout()
            g.fig.savefig(os.path.join(outdir, "zoom_s1_s2_s5_contigs_key_thresholds_pct_best.png"), dpi=300)
            plt.close(g.fig)

    plt.close()


def write_embedded_html_report(
    outdir: str,
    report_filename: str = "report.html",
    include_dirs: Optional[List[str]] = None,
    extra_png_paths: Optional[List[str]] = None,
) -> str:
    scan_dirs: List[str] = [outdir]
    if include_dirs:
        for path in include_dirs:
            if path and os.path.isdir(path):
                scan_dirs.append(path)

    image_paths: List[str] = []
    seen_paths = set()
    for scan_dir in scan_dirs:
        for root, _, files in os.walk(scan_dir):
            for name in files:
                if not name.lower().endswith(".png"):
                    continue
                path = os.path.join(root, name)
                if path in seen_paths:
                    continue
                seen_paths.add(path)
                image_paths.append(path)

    image_paths = sorted(image_paths)
    if extra_png_paths:
        for path in extra_png_paths:
            if path and os.path.isfile(path) and path not in seen_paths:
                seen_paths.add(path)
                image_paths.append(path)
    report_path = os.path.join(outdir, report_filename)

    sections: List[str] = []
    for image_path in image_paths:
        image_name = os.path.relpath(image_path, start=outdir)
        with open(image_path, "rb") as handle:
            encoded = base64.b64encode(handle.read()).decode("ascii")
        title = image_name.rsplit(".", 1)[0].replace("_", " ").title()
        sections.append(
            "\n".join(
                [
                    "<section class='plot'>",
                    f"  <h2>{escape(title)}</h2>",
                    f"  <img alt='{escape(image_name)}' src='data:image/png;base64,{encoded}' />",
                    "</section>",
                ]
            )
        )

    html = "\n".join(
        [
            "<!doctype html>",
            "<html lang='en'>",
            "<head>",
            "  <meta charset='utf-8' />",
            "  <meta name='viewport' content='width=device-width, initial-scale=1' />",
            "  <title>Summary Stats Report</title>",
            "  <style>",
            "    body { font-family: Arial, sans-serif; margin: 24px; background: #f8fafc; color: #0f172a; }",
            "    h1 { margin: 0 0 8px; }",
            "    p { margin: 0 0 18px; color: #334155; }",
            "    .plot { margin: 24px 0 32px; padding: 16px; border: 1px solid #cbd5e1; border-radius: 8px; background: #ffffff; }",
            "    .plot h2 { margin: 0 0 12px; font-size: 1.05rem; }",
            "    .plot img { max-width: 100%; height: auto; display: block; }",
            "  </style>",
            "</head>",
            "<body>",
            "  <h1>Summary Statistics Graph Report</h1>",
            f"  <p>Output directory: {escape(outdir)}<br />Embedded plots: {len(image_paths)}</p>",
            *(sections if sections else ["  <p>No PNG plots found to embed.</p>"]),
            "</body>",
            "</html>",
        ]
    )

    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(html)

    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot assembly stats from summary_stats_log.txt")
    parser.add_argument("--log", default="summary_stats_log.txt", help="Path to summary stats log")
    parser.add_argument("--outdir", default="plots/summary_stats", help="Directory to write plots and parsed CSV")
    parser.add_argument("--data-dir", default=os.path.join("data"), help="Directory to write parsed CSV files")
    parser.add_argument("--plots-root", default=None, help="Root directory containing plot subdirectories")
    parser.add_argument(
        "--ng50-target-bp",
        type=int,
        default=DEFAULT_NG50_TARGET_BP,
        help="Expected genome size in bp for NG50 (default: 5,000,000)",
    )
    parser.add_argument("--skip-report", action="store_true", help="Skip writing HTML report")
    parser.add_argument("--report-only", action="store_true", help="Only write HTML report from existing PNG outputs")
    args = parser.parse_args()

    plots_root_dir = args.plots_root if args.plots_root else os.path.dirname(args.outdir)
    all_plot_dirs = [
        os.path.join(plots_root_dir, "summary_stats"),
        os.path.join(plots_root_dir, "fragmentation"),
        os.path.join(plots_root_dir, "rarefaction"),
        os.path.join(plots_root_dir, "threshold_curves"),
        os.path.join(plots_root_dir, "yield_complexity"),
    ]
    extra_png_paths = [
        os.path.join(
            plots_root_dir,
            "threshold_curves",
            "assembly_type_total_length_vs_min_contig_length.png",
        )
    ]

    if args.report_only:
        os.makedirs(args.outdir, exist_ok=True)
        report_path = write_embedded_html_report(
            args.outdir,
            include_dirs=all_plot_dirs,
            extra_png_paths=extra_png_paths,
        )
        print(f"Wrote HTML report: {report_path}")
        return

    df = parse_summary_stats_log(args.log)
    if df.empty:
        raise SystemExit("No assembly blocks parsed from log. Check --log path and format.")

    threshold_df = parse_threshold_stats_log(args.log)
    df = add_recalculated_ng50(df, args.ng50_target_bp)

    if not df.empty:
        for col, digits in {
            "gc_fraction": 4,
            "sequence_total_mb": 3,
            "n50_kbp": 3,
            "ng50_kbp": 3,
            "max_scaffold_kbp": 3,
            "max_contig_kbp": 3,
        }.items():
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").round(digits)

    summary_dedup_cols = [
        "sample",
        "assembler",
        "gc_fraction",
        "scaffold_total",
        "sequence_total_mb",
        "n50_kbp",
        "ng50_kbp",
        "lg50",
        "l50",
        "max_scaffold_kbp",
        "max_contig_kbp",
        "scaffolds_gt_50kb",
    ]
    summary_subset = [col for col in summary_dedup_cols if col in df.columns]
    summary_rows_before = len(df)
    if summary_subset:
        df = df.drop_duplicates(subset=summary_subset, keep="first").copy()
    summary_rows_removed = summary_rows_before - len(df)

    threshold_dedup_cols = [
        "sample",
        "assembler",
        "threshold_label",
        "threshold_kbp",
        "num_scaffolds",
        "num_contigs",
        "scaffold_length_bp",
        "contig_length_bp",
        "contig_coverage_pct",
    ]
    threshold_subset = [col for col in threshold_dedup_cols if col in threshold_df.columns]
    threshold_rows_before = len(threshold_df)
    if threshold_subset:
        threshold_df = threshold_df.drop_duplicates(subset=threshold_subset, keep="first").copy()
    threshold_rows_removed = threshold_rows_before - len(threshold_df)

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.data_dir, exist_ok=True)
    df.to_csv(os.path.join(args.data_dir, "parsed_summary_stats.csv"), index=False)
    if not threshold_df.empty:
        threshold_df.to_csv(os.path.join(args.data_dir, "parsed_threshold_stats.csv"), index=False)
    make_plots(df, threshold_df, args.outdir)

    report_path = ""
    if not args.skip_report:
        report_path = write_embedded_html_report(
            args.outdir,
            include_dirs=all_plot_dirs,
            extra_png_paths=extra_png_paths,
        )

    print(f"Parsed {len(df)} assembly entries")
    print(f"Parsed {len(threshold_df)} threshold-table rows")
    print(f"Removed {summary_rows_removed} duplicate assembly rows")
    print(f"Removed {threshold_rows_removed} duplicate threshold-table rows")
    print(f"Wrote CSV: {os.path.join(args.data_dir, 'parsed_summary_stats.csv')}")
    if not threshold_df.empty:
        print(f"Wrote CSV: {os.path.join(args.data_dir, 'parsed_threshold_stats.csv')}")
    print(f"Wrote plots to: {args.outdir}")
    if report_path:
        print(f"Wrote HTML report: {report_path}")


if __name__ == "__main__":
    main()
