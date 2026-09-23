#!/usr/bin/env python3
"""
validate_minimap2.py

Diagnostic tool to check whether a minimap2 run used sensible parameters,
by inspecting the resulting (sorted, indexed) BAM file(s).

This does NOT re-run alignment. It inspects what minimap2 already produced
and flags patterns that usually mean the preset/parameters were wrong for
your data type. Pass multiple BAMs to get a side-by-side comparison table
(e.g. one BAM per parameter variant you're testing).

Usage:
    pip install pysam --break-system-packages   # if not already installed

    # single file, full report
    python validate_minimap2.py aligned.sorted.bam --preset map-ont

    # multiple files, comparison table across runs
    python validate_minimap2.py run1.bam run2.bam run3.bam \\
        --preset map-ont \\
        --labels baseline k13 no_secondary \\
        --csv comparison.csv

Requires: pysam
"""

import argparse
import sys
import statistics
from collections import Counter, defaultdict

try:
    import pysam
except ImportError:
    sys.exit("This script requires pysam. Install with:\n"
              "  pip install pysam --break-system-packages")

# Rough expectation table per preset. These are heuristics from common
# practice, not hard thresholds -- use them to spot outliers, not as
# pass/fail gates.
PRESET_EXPECTATIONS = {
    "map-ont":  {"typical_identity": (0.85, 0.98), "clip_warn": 0.20},
    "map-hifi": {"typical_identity": (0.98, 1.00), "clip_warn": 0.10},
    "map-pb":   {"typical_identity": (0.80, 0.95), "clip_warn": 0.20},
    "sr":       {"typical_identity": (0.95, 1.00), "clip_warn": 0.05},
    "asm5":     {"typical_identity": (0.99, 1.00), "clip_warn": 0.05},
    "asm10":    {"typical_identity": (0.95, 1.00), "clip_warn": 0.05},
    "asm20":    {"typical_identity": (0.90, 1.00), "clip_warn": 0.08},
}


def safe_div(a, b):
    return a / b if b else float("nan")


# ---------------------------------------------------------------------------
# Stats computation (pure -- no printing) so it can feed either the single
# -file report or the multi-file comparison table.
# ---------------------------------------------------------------------------

def compute_stats(bam_path, preset=None):
    bam = pysam.AlignmentFile(bam_path, "rb")

    total = 0
    unmapped = 0
    mapped = 0
    secondary = 0
    supplementary = 0
    primary_mapped = 0
    duplicate = 0

    mapq_counter = Counter()
    clip_fractions = []
    identities = []
    read_lengths = []
    insert_sizes = []
    chimeric_reads = set()
    per_ref_coverage_reads = defaultdict(int)

    for read in bam.fetch(until_eof=True):
        total += 1

        if read.is_unmapped:
            unmapped += 1
            continue

        mapped += 1
        per_ref_coverage_reads[read.reference_name] += 1

        if read.is_secondary:
            secondary += 1
        if read.is_supplementary:
            supplementary += 1
            chimeric_reads.add(read.query_name)
        if read.is_duplicate:
            duplicate += 1

        if not read.is_secondary and not read.is_supplementary:
            primary_mapped += 1
            mapq_counter[read.mapping_quality] += 1

            qlen = read.query_length or (read.infer_query_length() or 0)
            if qlen:
                read_lengths.append(qlen)

            if read.cigartuples:
                clipped = sum(length for op, length in read.cigartuples if op in (4, 5))
                total_len = sum(length for op, length in read.cigartuples if op in (0, 1, 4, 5, 7, 8))
                if total_len:
                    clip_fractions.append(clipped / total_len)

            try:
                nm = read.get_tag("NM")
                aligned_len = read.query_alignment_length
                if aligned_len:
                    identities.append(1 - (nm / aligned_len))
            except KeyError:
                pass

            if read.is_paired and read.is_proper_pair and read.template_length:
                tlen = abs(read.template_length)
                if 0 < tlen < 100000:
                    insert_sizes.append(tlen)

    bam.close()

    n_primary = sum(mapq_counter.values())
    mapq0 = mapq_counter.get(0, 0)
    mapq_low = sum(v for k, v in mapq_counter.items() if k < 10)

    stats = {
        "bam_path": bam_path,
        "preset": preset,
        "total": total,
        "mapped": mapped,
        "unmapped": unmapped,
        "mapped_pct": safe_div(mapped, total),
        "unmapped_pct": safe_div(unmapped, total),
        "primary_mapped": primary_mapped,
        "secondary": secondary,
        "supplementary": supplementary,
        "duplicate": duplicate,
        "n_chimeric_reads": len(chimeric_reads),
        "mapq_counter": mapq_counter,
        "n_primary_for_mapq": n_primary,
        "mapq0_pct": safe_div(mapq0, n_primary),
        "mapq_lt10_pct": safe_div(mapq_low, n_primary),
        "mean_clip": statistics.mean(clip_fractions) if clip_fractions else None,
        "median_clip": statistics.median(clip_fractions) if clip_fractions else None,
        "high_clip_pct": (sum(1 for c in clip_fractions if c > 0.3) / len(clip_fractions)) if clip_fractions else None,
        "mean_identity": statistics.mean(identities) if identities else None,
        "median_identity": statistics.median(identities) if identities else None,
        "min_identity": min(identities) if identities else None,
        "max_identity": max(identities) if identities else None,
        "mean_read_len": statistics.mean(read_lengths) if read_lengths else None,
        "median_read_len": statistics.median(read_lengths) if read_lengths else None,
        "n_reads_with_length": len(read_lengths),
        "mean_insert": statistics.mean(insert_sizes) if insert_sizes else None,
        "median_insert": statistics.median(insert_sizes) if insert_sizes else None,
        "stdev_insert": statistics.pstdev(insert_sizes) if len(insert_sizes) > 1 else None,
        "per_ref_coverage_reads": per_ref_coverage_reads,
    }
    return stats


def get_flags(stats):
    """Return a list of (severity, message) flags for a single stats dict."""
    flags = []
    preset = stats["preset"]

    if stats["unmapped_pct"] is not None and stats["unmapped_pct"] > 0.10:
        flags.append(("WARN", f">10% unmapped ({stats['unmapped_pct']:.1%}). Check reference, "
                               "contamination, or seeding parameters (-k/-w)."))

    if stats["mapq0_pct"] is not None and stats["mapq0_pct"] > 0.30:
        flags.append(("WARN", f"High MAPQ=0 fraction ({stats['mapq0_pct']:.1%}). Repetitive "
                               "region, -N too high, or reads too short/divergent."))

    if stats["mean_clip"] is not None:
        clip_warn = PRESET_EXPECTATIONS.get(preset, {}).get("clip_warn", 0.15)
        if stats["mean_clip"] > clip_warn:
            flags.append(("WARN", f"Mean clip fraction {stats['mean_clip']:.1%} exceeds "
                                   f"typical {clip_warn:.0%} for preset '{preset}'. Adapter "
                                   "contamination, chimeric reads, or wrong preset/scoring."))

    if stats["mean_identity"] is not None and preset in PRESET_EXPECTATIONS:
        lo, hi = PRESET_EXPECTATIONS[preset]["typical_identity"]
        if not (lo <= stats["mean_identity"] <= hi + 0.02):
            flags.append(("WARN", f"Mean identity {stats['mean_identity']:.2%} outside typical "
                                   f"{lo:.0%}-{hi:.0%} range for preset '{preset}'."))

    return flags


# ---------------------------------------------------------------------------
# Single-file report (verbose, human-readable)
# ---------------------------------------------------------------------------

def print_single_report(stats):
    print("=" * 70)
    print(f"BAM: {stats['bam_path']}")
    if stats["preset"]:
        print(f"Declared preset: {stats['preset']}")
    print("=" * 70)

    print(f"\nTotal records:            {stats['total']:,}")
    print(f"Mapped records:           {stats['mapped']:,} ({stats['mapped_pct']:.1%})")
    print(f"Unmapped records:         {stats['unmapped']:,} ({stats['unmapped_pct']:.1%})")
    print(f"Primary mapped reads:     {stats['primary_mapped']:,}")
    print(f"Secondary alignments:     {stats['secondary']:,}")
    print(f"Supplementary alignments: {stats['supplementary']:,}")
    print(f"Duplicate-flagged:        {stats['duplicate']:,}")
    print(f"Reads w/ supplementary (possible chimera/SV): {stats['n_chimeric_reads']:,}")

    mapq_counter = stats["mapq_counter"]
    if mapq_counter:
        n_primary = stats["n_primary_for_mapq"]
        print(f"\nMAPQ distribution (primary alignments, n={n_primary:,}):")
        for q in sorted(mapq_counter):
            bar = "#" * int(50 * mapq_counter[q] / n_primary)
            print(f"  MAPQ {q:>3}: {mapq_counter[q]:>8,} {bar}")
        print(f"\n  MAPQ==0 fraction:  {stats['mapq0_pct']:.1%}")
        print(f"  MAPQ<10 fraction:  {stats['mapq_lt10_pct']:.1%}")

    if stats["mean_clip"] is not None:
        print(f"\nSoft/hard-clip fraction of aligned read length:")
        print(f"  mean:   {stats['mean_clip']:.1%}")
        print(f"  median: {stats['median_clip']:.1%}")
        print(f"  % of reads with >30% clipped: {stats['high_clip_pct']:.1%}")

    if stats["mean_identity"] is not None:
        print(f"\nAlignment identity (from NM tag):")
        print(f"  mean:   {stats['mean_identity']:.2%}")
        print(f"  median: {stats['median_identity']:.2%}")
        print(f"  min/max: {stats['min_identity']:.2%} / {stats['max_identity']:.2%}")

    if stats["mean_read_len"] is not None:
        print(f"\nAligned read length: mean={stats['mean_read_len']:.0f}, "
              f"median={stats['median_read_len']:.0f}, n={stats['n_reads_with_length']:,}")

    if stats["mean_insert"] is not None:
        print(f"\nInsert size (proper pairs): mean={stats['mean_insert']:.0f}, "
              f"median={stats['median_insert']:.0f}, stdev={stats['stdev_insert'] or 0:.0f}")

    per_ref = stats["per_ref_coverage_reads"]
    if len(per_ref) > 1:
        print(f"\nReads per reference contig (top 10 of {len(per_ref)}):")
        for ref, cnt in sorted(per_ref.items(), key=lambda x: -x[1])[:10]:
            print(f"  {ref}: {cnt:,}")

    flags = get_flags(stats)
    if flags:
        print("\n" + "-" * 70)
        print("Flags:")
        for severity, msg in flags:
            print(f"  [{severity}] {msg}")

    print("\n" + "=" * 70)
    print("Interpretation notes:")
    print("- No single metric proves parameters are 'right' -- look for internal")
    print("  consistency: does identity/clip/MAPQ match what's expected for your")
    print("  read type and preset? Big mismatches are the signal, not the absolute values.")
    print("=" * 70)


# ---------------------------------------------------------------------------
# Multi-file comparison
# ---------------------------------------------------------------------------

COMPARISON_ROWS = [
    ("Mapped %",        lambda s: s["mapped_pct"],       "{:.1%}"),
    ("Unmapped %",      lambda s: s["unmapped_pct"],      "{:.1%}"),
    ("Primary mapped",  lambda s: s["primary_mapped"],    "{:,}"),
    ("Secondary aln",   lambda s: s["secondary"],         "{:,}"),
    ("Supplementary",   lambda s: s["supplementary"],     "{:,}"),
    ("Chimeric reads",  lambda s: s["n_chimeric_reads"],  "{:,}"),
    ("MAPQ==0 %",       lambda s: s["mapq0_pct"],         "{:.1%}"),
    ("MAPQ<10 %",       lambda s: s["mapq_lt10_pct"],     "{:.1%}"),
    ("Mean clip %",     lambda s: s["mean_clip"],         "{:.1%}"),
    ("Mean identity",   lambda s: s["mean_identity"],     "{:.2%}"),
    ("Median identity", lambda s: s["median_identity"],   "{:.2%}"),
    ("Mean read len",   lambda s: s["mean_read_len"],     "{:.0f}"),
    ("Mean insert",     lambda s: s["mean_insert"],       "{:.0f}"),
]


def fmt(value, fmt_str):
    if value is None:
        return "-"
    try:
        return fmt_str.format(value)
    except (ValueError, TypeError):
        return str(value)


def print_comparison_table(all_stats, labels):
    col_width = max(14, max(len(l) for l in labels) + 2)
    label_col = max(16, max(len(name) for name, _, _ in COMPARISON_ROWS) + 2)

    print("=" * 70)
    print("COMPARISON ACROSS RUNS")
    print("=" * 70)

    header = " " * label_col + "".join(f"{l:>{col_width}}" for l in labels)
    print(header)
    print("-" * len(header))

    for name, getter, fmt_str in COMPARISON_ROWS:
        row_vals = [getter(s) for s in all_stats]
        cells = "".join(f"{fmt(v, fmt_str):>{col_width}}" for v in row_vals)
        print(f"{name:<{label_col}}{cells}")

    # Flag the metrics with the largest relative spread across runs --
    # these are the params most likely "doing something" to your alignment.
    print("\n" + "-" * 70)
    print("Biggest differences across runs (by relative spread):")
    spreads = []
    for name, getter, fmt_str in COMPARISON_ROWS:
        vals = [getter(s) for s in all_stats if getter(s) is not None]
        if len(vals) >= 2 and max(vals) != 0:
            spread = (max(vals) - min(vals)) / (abs(max(vals)) if max(vals) != 0 else 1)
            spreads.append((spread, name, min(vals), max(vals), fmt_str))
    spreads.sort(reverse=True)
    for spread, name, lo, hi, fmt_str in spreads[:6]:
        if spread > 0.01:
            print(f"  {name}: {fmt(lo, fmt_str)} to {fmt(hi, fmt_str)}  (Δ {spread:.0%})")

    print("\n" + "-" * 70)
    print("Per-run flags:")
    for label, s in zip(labels, all_stats):
        flags = get_flags(s)
        if flags:
            print(f"  {label}:")
            for severity, msg in flags:
                print(f"    [{severity}] {msg}")
        else:
            print(f"  {label}: no flags raised")

    print("=" * 70)


def write_csv(all_stats, labels, csv_path):
    import csv as csv_mod
    with open(csv_path, "w", newline="") as fh:
        writer = csv_mod.writer(fh)
        writer.writerow(["metric"] + labels)
        for name, getter, _ in COMPARISON_ROWS:
            writer.writerow([name] + [getter(s) for s in all_stats])
    print(f"\nComparison table written to {csv_path}")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bams", nargs="+", help="One or more sorted, indexed BAM files from minimap2")
    ap.add_argument("--preset", nargs="+",
                     help="Preset(s) used. Give one value to apply to all BAMs, or one per BAM "
                          "in the same order (e.g. --preset map-ont map-ont sr)")
    ap.add_argument("--labels", nargs="+",
                     help="Short names for each BAM, used as column headers in comparison mode "
                          "(defaults to filenames)")
    ap.add_argument("--csv", help="Write the comparison table to this CSV path (multi-BAM mode only)")
    args = ap.parse_args()

    n = len(args.bams)

    if args.preset and len(args.preset) not in (1, n):
        sys.exit(f"--preset given {len(args.preset)} values but {n} BAM(s) provided; "
                  "give exactly one (applied to all) or one per BAM.")
    presets = (args.preset * n) if (args.preset and len(args.preset) == 1) else \
              (args.preset if args.preset else [None] * n)

    if args.labels and len(args.labels) != n:
        sys.exit(f"--labels given {len(args.labels)} values but {n} BAM(s) provided.")
    labels = args.labels if args.labels else [b.rsplit("/", 1)[-1] for b in args.bams]

    all_stats = [compute_stats(bam, preset=preset) for bam, preset in zip(args.bams, presets)]

    if n == 1:
        print_single_report(all_stats[0])
    else:
        print_comparison_table(all_stats, labels)
        if args.csv:
            write_csv(all_stats, labels, args.csv)


if __name__ == "__main__":
    main()