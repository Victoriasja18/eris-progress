
# This is a script to calculate the precision and recall
import ast
import pandas as pd
from collections import defaultdict
import IS_dict as isd
from scipy.optimize import linear_sum_assignment
import numpy as np


# PATH 
FA = "../../data/EC958_complete_genome.fasta"
GT = "../../data/complete_ec958_IS.reanotated.sorted.bed"
GFA = "../../data/contigs_shovill_simulated.gfa"

def format_fa(df, context: str):

    """
    Context whether you want inside, downstream, upstream
    """
    clean_fa = pd.read_csv(df, sep="\t")

    clean_fa = clean_fa.iloc[:, :4]
    # Extract the locations
    clean_fa["transcript"] = clean_fa['gene_id'].str.split('_').str[0]
    clean_fa["chromStart"] = clean_fa['gene_id'].str.split('_').str[1].astype(int)
    clean_fa["chromEnd"] = clean_fa['gene_id'].str.split('_').str[2].astype(int)

    # Extract the name
    clean_fa["name"] = clean_fa["target"].str.split("_").str[0]
    clean_fa["repeats"] = clean_fa.groupby("name")["locus_id"].transform(lambda x: pd.factorize(x)[0] + 1)

    # Exclude unnecessary columns
    clean_fa = clean_fa.drop(columns=["gene_id", "locus_id", "target"])

    # Just collect what is inside
    clean_fa = clean_fa[clean_fa["context"] == context]

    # Sort the location
    clean_fa = clean_fa.sort_values(["transcript", "name", "repeats", "chromStart"]).reset_index(drop=True)

    # Collapsing the location
    gap = 1000  # max distance between rows to be considered the same instance, adjust as needed

    group_keys = ["context", "transcript", "name", "repeats"]
    new_group = (
        (clean_fa[group_keys] != clean_fa[group_keys].shift()).any(axis=1)
        | (clean_fa["chromStart"] > clean_fa["chromEnd"].shift() + gap)
    )
    clean_fa["cluster_id"] = new_group.cumsum()

    agg = {"chromStart": "min", "chromEnd": "max"}
    for c in clean_fa.columns:
        if c not in group_keys and c not in agg and c != "cluster_id":
            agg[c] = "first"

    clean_fa = clean_fa.groupby(["cluster_id"] + group_keys, dropna=False, as_index=False).agg(agg)
    clean_fa = clean_fa.drop(columns="cluster_id")

    # Synonyms
    rows = isd.load_is_table()
    name_to_synonyms, synonym_to_name = isd.build_synonym_dicts(rows)
    clean_fa["synonyms"] = clean_fa["name"].apply(
    lambda n: isd.get_synonyms(n, name_to_synonyms, synonym_to_name) or pd.NA)

    return clean_fa

def format_gt(df):
    # Build the lookup once
    rows = isd.load_is_table()
    name_to_synonyms, synonym_to_name = isd.build_synonym_dicts(rows)

    bed_columns = ['transcript', 'chromStart', 'chromEnd', 'name', 'score', 'strand']

    df = pd.read_table(df, header=None, names=bed_columns[:6])

    gt = df.copy()
    gt["repeats"] = gt["name"].str.split("_").str[1]
    gt["name"] = gt["name"].str.split("_").str[0]

    gt["synonyms"] = gt["name"].apply(
        lambda n: isd.get_synonyms(n, name_to_synonyms, synonym_to_name) or pd.NA)

    return gt


def _parse_synonyms(val):
    """Coerce a synonyms field (list, stringified list, or NaN) into a list."""
    if isinstance(val, str):
        try:
            val = ast.literal_eval(val)
        except (ValueError, SyntaxError):
            return []
    if isinstance(val, list):
        return val
    return []


def _names_match(pred, truth):
    """True if pred and truth refer to the same IS element by name or synonym,
    checked in both directions."""
    pred_name = pred["name"]
    truth_name = truth["name"]

    if pred_name == truth_name:
        return True

    pred_syns = _parse_synonyms(pred.get("synonyms", pd.NA))
    truth_syns = _parse_synonyms(truth.get("synonyms", pd.NA))

    if truth_name in pred_syns:
        return True
    if pred_name in truth_syns:
        return True
    if pred_syns and truth_syns and set(pred_syns) & set(truth_syns):
        return True

    return False


_PENALTY = 1e9  # cost for a non-edge; large enough to never beat a real match


def precision_recall(df1, df2, buffer):
    """
    df1 = ground truth, df2 = predictions.
    A prediction is TP only if its location overlaps a truth (within buffer)
    AND its name matches the truth's name or either side's synonym list.
    Matching is a global one-to-one optimal assignment (Hungarian algorithm):
    among all ways to pair eligible predictions to truths, it picks the one
    with the most pairs, breaking ties by total overlap. A prediction that
    loses that assignment (or never had an eligible truth) is FP.
    """
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    # --- Stage 1: location candidates ---
    by_pred = defaultdict(list)
    for i, pred in df2.iterrows():
        for j, truth in df1.iterrows():
            if pred["transcript"] != truth["transcript"]:
                continue
            if (pred["chromStart"] - buffer) <= truth["chromEnd"] and \
               (pred["chromEnd"] + buffer) >= truth["chromStart"]:
                by_pred[i].append(j)

    def overlap_size(pred, truth):
        return min(pred["chromEnd"], truth["chromEnd"]) - max(pred["chromStart"], truth["chromStart"])

    # --- Stage 2: name-eligible edges among location candidates ---
    pred_ids = sorted(by_pred.keys())
    truth_ids = sorted({j for js in by_pred.values() for j in js})
    pred_pos = {i: r for r, i in enumerate(pred_ids)}
    truth_pos = {j: c for c, j in enumerate(truth_ids)}

    fp_reason = {}  # i -> "wrong_name" once we know it has no name-matching candidate
    cost = np.full((len(pred_ids), len(truth_ids)), _PENALTY)

    for i in pred_ids:
        pred = df2.loc[i]
        any_name_match = False
        for j in by_pred[i]:
            truth = df1.loc[j]
            if _names_match(pred, truth):
                any_name_match = True
                cost[pred_pos[i], truth_pos[j]] = -overlap_size(pred, truth)
        if not any_name_match:
            fp_reason[i] = "wrong_name"

    # --- Stage 3: global optimal one-to-one matching (max pairs, then max total overlap) ---
    matched_df1 = set()
    matched_df2 = set()
    if cost.size:
        row_idx, col_idx = linear_sum_assignment(cost)
        for r, c in zip(row_idx, col_idx):
            if cost[r, c] >= _PENALTY:
                continue  # this slot had no real eligible edge, ignore it
            matched_df2.add(pred_ids[r])
            matched_df1.add(truth_ids[c])

    for i in pred_ids:
        if i not in matched_df2 and i not in fp_reason:
            fp_reason[i] = "duplicate"  # had a name match, but lost the optimal assignment to another prediction

    # --- metrics ---
    tp = len(matched_df2)
    fn = len(df1) - len(matched_df1)

    no_candidate_idx = set(df2.index) - set(by_pred.keys())
    fp = len(no_candidate_idx) + len(fp_reason)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    summary = {
        "buffer": buffer, "TP": tp, "FP": fp, "FN": fn,
        "precision": round(precision, 3), "recall": round(recall, 3), "f1": round(f1, 3),
    }

    # --- detail dataframe ---
    tp_df = df2.loc[list(matched_df2)].copy()
    tp_df[f"buffer_{buffer}"] = "TP"

    fp_no_candidate_df = df2.loc[list(no_candidate_idx)].copy()
    fp_no_candidate_df[f"buffer_{buffer}"] = "FP"

    fp_named_records = []
    for i, reason in fp_reason.items():
        record = df2.loc[i].to_dict()
        record[f"buffer_{buffer}"] = "FP"
        record["fp_reason"] = reason
        fp_named_records.append(record)
    fp_named_df = pd.DataFrame(fp_named_records)

    fn_df = df1.loc[df1.index.difference(matched_df1)].copy()
    fn_df[f"buffer_{buffer}"] = "FN"

    detail_df = pd.concat(
        [tp_df, fp_no_candidate_df, fp_named_df, fn_df],
        ignore_index=True,
    )

    return summary, detail_df



