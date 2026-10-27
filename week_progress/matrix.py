# Precision / recall of eris loci against the truth BED.

"""
# Usage
import matrix as mt

# Complete genome (chromosome + plasmid)
gt = mt.read_truth(mt.GT_COMPLETE)                    # or read_truth(..., exclude=["ISSen10"])
fa = mt.read_report("../result/fa_check/fa_check_report.tsv")
summary, detail = mt.precision_recall(gt, fa, buffer=0)

# Simulated GFAs were assembled from the chromosome only
gt = mt.read_truth(mt.GT_CHROMOSOME)
gfa = mt.read_report("../result/new_merge_id/new_merge_id_report.tsv")
summary, detail = mt.precision_recall_gfa(gt, gfa, count_copies=True)
"""

import csv
from collections import defaultdict
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

import IS_dict as isd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "data"

GT_COMPLETE = DATA / "complete_ec958_IS.reanotated.sorted.bed"  # chromosome + plasmid (HG941719.1)
GT_CHROMOSOME = DATA / "ec958_IS.reannotated.sorted.bed"         # chromosome only
FA = DATA / "EC958_complete_genome.fasta"
GFA = DATA / "contigs_shovill_simulated.gfa"
GFA_FASTQ = DATA / "contigs_shovill_realfastq.gfa"

# Same ISfinder table eris --merge-synonyms uses, so both see the same synonyms
IS_TABLE = Path.home() / ".cache" / "eris" / "IS.csv"


# ---- Reading files ---------------------------------------------------------------------------------------------------

@cache
def _synonym_to_name():
    if IS_TABLE.exists():
        with open(IS_TABLE, newline="") as handle:
            rows = list(csv.reader(handle))
    else:
        rows = isd.load_is_table()
    return isd.build_synonym_dicts(rows)[1]


def is_name(target):
    # ISEc23_IS66_unknown_stitched / ISEc23_1 -> ISEc23
    return target.split("_")[0]


def is_group(name):
    return _synonym_to_name().get(name, name)


def read_truth(bed, exclude=()):
    truth = pd.read_table(bed, header=None, usecols=[0, 1, 2, 3],
                          names=["contig", "start", "end", "target"])
    truth["name"] = truth["target"].map(is_name)
    truth = truth[~truth["name"].isin(exclude)].copy()
    truth["group"] = truth["name"].map(is_group)
    return truth[["contig", "start", "end", "name", "group"]].reset_index(drop=True)


def read_report(report):
    loci = pd.read_csv(report, sep="\t").drop_duplicates("locus_id")
    if "locus_start" not in loci:
        raise ValueError(f"{report} has no locus coordinates; rerun eris to get them")

    # collapse mode calls the depth column copy_number, variant mode fractional_depth
    copies = "copy_number" if "copy_number" in loci else "fractional_depth"

    # A locus can hold several IS (comma separated): one row each
    loci["target"] = loci["target"].str.split(",")
    loci = loci.explode("target")

    pred = pd.DataFrame({
        "locus_id": loci["locus_id"],
        "contig": loci["locus_contig"].astype(str),
        "start": loci["locus_start"],
        "end": loci["locus_end"],
        "name": loci["target"].map(is_name),
        "copies": loci[copies].astype(float).round().astype(int),
    })
    pred["group"] = pred["name"].map(is_group)
    return pred.reset_index(drop=True)


# ---- Precision / recall ----------------------------------------------------------------------------------------------

_NO_EDGE = 1e9  # cost of an impossible pair; never beats a real match


def _scores(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"TP": int(tp), "FP": int(fp), "FN": int(fn),
            "precision": round(float(precision), 3),
            "recall": round(float(recall), 3),
            "f1": round(float(f1), 3)}


def precision_recall(truth, pred, buffer=0):
    # Complete genome: a prediction is TP when it overlaps a truth on the same contig
    # (within `buffer` bp) AND has the same IS group. Pairs are one-to-one, chosen by
    # the Hungarian algorithm: most pairs first, then most total overlap.
    truth = truth.reset_index(drop=True)
    pred = pred.reset_index(drop=True)

    cost = np.full((len(pred), len(truth)), _NO_EDGE)
    reason = {}
    for i, p in pred.iterrows():
        overlapping = truth[(truth["contig"] == p["contig"])
                            & (truth["start"] < p["end"] + buffer)
                            & (truth["end"] > p["start"] - buffer)]
        same_group = overlapping[overlapping["group"] == p["group"]]
        if overlapping.empty:
            reason[i] = "no_overlap"
        elif same_group.empty:
            reason[i] = "wrong_name"
        for j, t in same_group.iterrows():
            cost[i, j] = -(min(p["end"], t["end"]) - max(p["start"], t["start"]))

    matched = defaultdict(set)
    if cost.size:
        for i, j in zip(*linear_sum_assignment(cost)):
            if cost[i, j] < _NO_EDGE:
                matched["pred"].add(i)
                matched["truth"].add(j)

    for i in pred.index:
        if i not in matched["pred"] and i not in reason:
            reason[i] = "duplicate"  # right place and name, but another prediction took that truth

    tp, fp, fn = len(matched["pred"]), len(reason), len(truth) - len(matched["truth"])
    summary = {"buffer": buffer, **_scores(tp, fp, fn)}

    pred = pred.assign(result="FP", fp_reason=pd.Series(reason, dtype=object))
    pred.loc[list(matched["pred"]), "result"] = "TP"
    missed = truth.drop(index=list(matched["truth"])).assign(result="FN")
    detail = pd.concat([pred, missed], ignore_index=True)

    return summary, detail


def precision_recall_gfa(truth, pred, count_copies=True):
    # Assembly graph: loci have no genome coordinates, so match on IS group only.
    # count_copies=True  -> truth: one copy per BED line, eris: sum of copy numbers
    # count_copies=False -> presence only, each IS counted once on each side
    counts = pd.concat([
        truth["group"].value_counts().rename("truth"),
        pred.groupby("group")["copies"].sum().rename("eris"),
    ], axis=1).fillna(0).astype(int)
    if not count_copies:
        counts = counts.clip(upper=1)

    counts["TP"] = counts[["truth", "eris"]].min(axis=1)
    counts["FP"] = (counts["eris"] - counts["truth"]).clip(lower=0)
    counts["FN"] = (counts["truth"] - counts["eris"]).clip(lower=0)
    counts["names"] = pd.concat([truth, pred]).groupby("group")["name"].agg(lambda n: ",".join(sorted(set(n))))

    summary = {"Copies": "Counted Copies" if count_copies else "Only by name",
               **_scores(counts["TP"].sum(), counts["FP"].sum(), counts["FN"].sum())}
    detail = counts.rename_axis("IS").reset_index().sort_values("IS", ignore_index=True)

    return summary, detail
