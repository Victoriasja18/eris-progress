
# This is a script to calculate the precision and recall
import ast
import pandas as pd
from collections import defaultdict
import IS_dict as isd

def format_fa(df, context: str):
    clean_fa = pd.read_csv(df, sep="\t")

    clean_fa = clean_fa.iloc[:,:4]

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


def precision_recall(df1, df2, buffer):
    """
    df1 = ground truth, df2 = predictions.
    A prediction is TP only if its location overlaps a truth (within buffer)
    AND its name matches the truth's name or either side's synonym list.
    Any prediction whose location matches but whose name doesn't is scored as FP.
    """
    # --- Stage 1: location candidates ---
    by_pred = defaultdict(list)
    for i, pred in df2.iterrows():
        for j, truth in df1.iterrows():
            if pred["transcript"] != truth["transcript"]:
                continue
            if (pred["chromStart"] - buffer) <= truth["chromEnd"] and \
               (pred["chromEnd"] + buffer) >= truth["chromStart"]:
                by_pred[i].append(j)

    # --- Stage 2: name check within each prediction's candidates ---
    matched_df1 = set()
    matched_df2 = set()
    fp_wrong_name_records = []

    for i, truth_idxs in by_pred.items():
        pred = df2.loc[i]
        tp_found = False

        for j in truth_idxs:
            truth = df1.loc[j]

            if _names_match(pred, truth):
                matched_df1.add(j)
                matched_df2.add(i)
                tp_found = True
                break  # one correct match is enough for this prediction

        if not tp_found:
            record = pred.to_dict()
            record[f"buffer_{buffer}"] = "FP"
            record["fp_reason"] = "wrong_name"
            fp_wrong_name_records.append(record)

    # --- metrics ---
    tp = len(matched_df2)
    fn = len(df1) - len(matched_df1)

    no_candidate_idx = set(df2.index) - set(by_pred.keys())      # no location overlap at all
    wrong_name_idx   = set(by_pred.keys()) - matched_df2          # location matched, name didn't
    fp_idx = no_candidate_idx | wrong_name_idx
    fp = len(fp_idx)

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

    fp_wrong_name_df = pd.DataFrame(fp_wrong_name_records)

    fn_df = df1.loc[df1.index.difference(matched_df1)].copy()
    fn_df[f"buffer_{buffer}"] = "FN"

    detail_df = pd.concat(
        [tp_df, fp_no_candidate_df, fp_wrong_name_df, fn_df],
        ignore_index=True,
    )

    return summary, detail_df





