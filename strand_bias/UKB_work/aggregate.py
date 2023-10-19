
import pandas as pd
from pandas import IndexSlice as idx


VCOUNT_HEADER = ["PSC", "id", "sample_id", "nRefHom", "nNonRefHom",
                 "nHets", "nTransitions", "nTransversions", "nIndels",
                 "averagedepth", "nSingletons", "nHapRef", "nHapAlt",
                 "nMissing"]

COMPLEMENT_DICT = {"A": "T", "C": "G", "G": "C", "T": "A"}
VARIANT_COMPLEMENT_DICT = {'AtoC': 'TtoG',
                           'AtoG': 'TtoC',
                           'AtoT': 'TtoA',
                           'CtoA': 'GtoT',
                           'CtoG': 'GtoC',
                           'CtoT': 'GtoA',
                           'GtoA': 'CtoT',
                           'GtoC': 'CtoG',
                           'GtoT': 'CtoA',
                           'TtoA': 'AtoT',
                           'TtoC': 'AtoG',
                           'TtoG': 'AtoC'}
VARIANTS = list(VARIANT_COMPLEMENT_DICT.keys())


def read_variant_counts(file, preprocessed=False):
    if preprocessed:
        header = VCOUNT_HEADER[2:9] + VCOUNT_HEADER[10:]
        cols = None
    else:
        header = VCOUNT_HEADER
        cols = set(range(14))-{0, 1, 9}
    data = pd.read_csv(file, sep="\t", header=0, names=header,
                       usecols=cols, index_col=0)
    return data


def read_sample_mapping(file):
    data = pd.read_csv(file, sep="\t", header=0, index_col=0)
    return data


def aggregate_chr_variant_counts(files, preprocessed=False):
    data = read_variant_counts(files[0], preprocessed)
    for file in files[1:]:
        data += read_variant_counts(file)
    return data


def encode_sample_ids(variant_count_file, mapping_file, preprocessed=True):
    counts = read_variant_counts(variant_count_file, preprocessed)
    mapping = read_sample_mapping(mapping_file)
    merged = pd.merge(counts, mapping, left_index=True, right_index=True)
    merged = merged.set_index("key")
    return merged


def read_multidex_vcounts(file):
    data = pd.read_csv(file, sep="\t", header=0, index_col=[0, 1, 2])
    return data


def reverse_complement_variants(variants):
    rev_variants = [VARIANT_COMPLEMENT_DICT[x] for x in variants]
    return rev_variants


def sort_variant_index(data):
    data.sort_index(key=lambda x: ["".join({y[0], y[-1]}) for y in x])
    return data


def calculate_complemented_totals(data):
    coding = (data.loc[idx[:, "forward", "reference", :]]
              + data.loc[idx[:, "reverse", "complement", :]])
    template = (data.loc[idx[:, "reverse", "reference", :]]
                + data.loc[idx[:, "forward", "complement", :]])
    return coding, template
