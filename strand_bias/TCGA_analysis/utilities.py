
from collections import Counter
import numpy as np


_COMP_DICT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _normalize_minimum(minimum_raw_count, normalize_by_nt_content=False,
                       normalization_factor=1, normalization_data=None):
    if minimum_raw_count is None:
        minimums = {nt: 0 for nt in "ACGT"}
    elif normalize_by_nt_content:
        minimums = minimum_raw_count / normalization_data.sum(axis=0) * normalization_factor
    else:
        minimums = {nt: minimum_raw_count for nt in "ACGT"}
    return minimums


def rectify_ratios(data):
    return np.exp(np.abs(np.log(data)))


def make_complement(seq):
    comp = "".join([_COMP_DICT[x] for x in seq])
    return comp


def reverse_complement(seq):
    return make_complement(seq)[::-1]


_CHANGES = ["TC", "CT", "TA", "TG", "CG", "GT"]
_COMPS = [make_complement(change) for change in _CHANGES]


def calculate_tukey_values(data, columns=None, iqr_factor=1.5):
    if columns is not None:
        data = data[columns]
    metrics = data.describe()
    fences = {col: (metrics[col]["25%"] - iqr_factor * (iqr := metrics[col]["75%"] - metrics[col]["25%"]),
                    metrics[col]["75%"] + iqr_factor * iqr) for col in metrics.columns}
    return fences


def filter_by_fences(table, columns=None, iqr_factor=1.5, return_outliers=False):
    if columns is None:
        columns = table.columns
    fences = calculate_tukey_values(table, columns, iqr_factor)
    if not return_outliers:
        for col in columns:
            table = table.loc[(fences[col][0] <= table[col]) & (table[col] <= fences[col][1])]
    else:
        for col in columns:
            table = table.loc[(table[col] < fences[col][0]) | (fences[col][1] < table[col])]
    return table


def filter_raw_call_data_dups(call_table):
    calls = (call_table.groupby(["project_id", "case_id", "file_id"]).sum())
    id_counts = Counter([x[1] for x in calls.index.to_list()])
    dups = [x for x, y in id_counts.items() if y > 1]
    calls = (call_table
             .reset_index()
             .set_index("case_id")
             .drop(dups)
             .reset_index()
             .set_index(["project_id", "case_id", "file_id", "contig", "pos",
                         "gene_orientation", "ref", "alt"]))
    return calls
