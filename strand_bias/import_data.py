
import pandas as pd
from pandas.api.types import CategoricalDtype
import plotly.io as pio
from numpy import log2

pio.renderers.default = "browser"


def read_asymmetry_data(asymmetry_file):
    data = pd.read_csv(asymmetry_file, sep="\t", index_col=0)
    data["Burden_per_1k_reads_log2"] = log2(1000*data.Nbdn)
    return data


def read_mapped_singleton_counts(mapped_singleton_file):
    data = pd.read_csv(mapped_singleton_file, sep="\t", index_col=0)
    return data


def read_mapped_genotype_counts(geno_counts):
    data = pd.read_csv(geno_counts, sep="\t", index_col=0)
    return data


def read_ancestry_file(ancestry_file):
    data = pd.read_csv(ancestry_file, sep="\t", index_col=0)
    return data


def read_reference_nuc_file(nuc_file):
    data = pd.read_csv(nuc_file, sep="\t", header=0,
                       names=["chr", "start", "stop", "targets", "score", "strand",
                              "pct_at", "pct_gc", "A_count", "C_count", "G_count",
                              "T_count", "N_count", "other_count", "length"])
    return data


def normalize_mismatch_counts(asymmetry_data, nuc_counts):
    asymmetry_data["CodingGT_norm"] = asymmetry_data.CodingGT / nuc_counts.G_count.sum()
    asymmetry_data["CodingCA_norm"] = asymmetry_data.CodingCA / nuc_counts.C_count.sum()

    asymmetry_data["TemplateGT_norm"] = asymmetry_data.TemplateGT / nuc_counts.G_count.sum()
    asymmetry_data["TemplateCA_norm"] = asymmetry_data.TemplateCA / nuc_counts.C_count.sum()

    # asymmetry_data["PosCA_norm"] = asymmetry_data.PosCA / nuc_counts.C_count.sum()
    # asymmetry_data["NegCA_norm"] = asymmetry_data.NegCA / nuc_counts.G_count.sum()
    # asymmetry_data["CA_ratio_norm"] = (asymmetry_data.PosCA_norm
    #                                    / asymmetry_data.NegCA_norm)
    return asymmetry_data


def calculate_asym_ratios(asymmetry_data, normalized=True):
    gt, ca = [f"{mut}_norm" if normalized else mut for mut in ["GT", "CA"]]
    asym = "Asymmetry_norm" if normalized else "Asymmetry"
    # for mut in [gt, ca]:
    #     asymmetry_data[f"{mut}_ratio"] = (asymmetry_data[f"Coding{mut}"]
    #                                       / asymmetry_data[f"Template{mut}"])
    #     asymmetry_data[f"{mut}_ratio_log"] = log2(asymmetry_data[f"{mut}_ratio"])
    asymmetry_data[asym] = (
        (asymmetry_data[f"Coding{ca}"] + asymmetry_data[f"Template{gt}"])
        / (asymmetry_data[f"Coding{gt}"] + asymmetry_data[f"Template{ca}"])
        )
    asymmetry_data[f"{asym}_log2"] = log2(asymmetry_data[asym])
    return asymmetry_data


def adjust_batches_in_asymmetry_data(data):
    batch_burden_median = dict()
    batch_asym_median = dict()
    batches = set(data.batch)
    for batch in batches:
        batch_data = data.loc[data.batch == batch]
        batch_burden_median[batch] = batch_data.Nbdn.median()
        batch_asym_median[batch] = batch_data.Asymmetry_norm.median()

    data["bpkr_batch_adjusted"] = data.apply(func=lambda x: x["Nbdn"] / batch_burden_median[x["batch"]],
                                             axis=1)
    data["bpkr_batch_adjusted_log2"] = log2(data["bpkr_batch_adjusted"])

    data["asym_batch_adjusted"] = data.apply(func=lambda x: x["Asymmetry_norm"] / batch_asym_median[x["batch"]],
                                             axis=1)
    data["asym_batch_adjusted_log2"] = log2(data["asym_batch_adjusted"])
    return data


def read_all_asymmetry_data(asym_file, singleton_file, geno_counts, nuc_file,
                            ancestry_file):
    asym = read_asymmetry_data(asym_file)

    scounts = read_mapped_singleton_counts(singleton_file)
    geno_counts = read_mapped_genotype_counts(geno_counts)
    ancestry = read_ancestry_file(ancestry_file)

    data = pd.merge(scounts, ancestry, left_index=True, right_index=True)
    data = pd.merge(data, asym, left_index=True, right_index=True)
    data = pd.merge(data, geno_counts, left_index=True, right_index=True)

    nuc_counts = read_reference_nuc_file(nuc_file)

    data = normalize_mismatch_counts(data, nuc_counts)
    data = calculate_asym_ratios(data)
    data = adjust_batches_in_asymmetry_data(data)

    return data


def read_vcounts(vcounts_file, master_file=True):
    data = pd.read_csv(vcounts_file, sep="\t", header=0)
    if not master_file:
        data["orientation"] = ["complement" if "complement" in x else "reference"
                               for x in data["coding_strand"]]
        data["coding_strand"] = [x.split("_")[0] for x in data["coding_strand"]]

    cat_cols = ["coding_strand", "orientation", "variant"]
    variants = ['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT',
                'GtoA', 'GtoC', 'GtoT', 'TtoA', 'TtoC', 'TtoG']
    dtypes = [["forward", "reverse"], ["reference", "complement"], variants]
    for col, dtype in zip(cat_cols, dtypes):
        data[col] = data[col].astype(CategoricalDtype(dtype))

    data = data.set_index(["sample_id"] + cat_cols)
    # data["nNonRefAlleles"]
    return data


def read_default_asym_data():
    wd = "~/Documents/Projects/StrandBias/"
    data = read_all_asymmetry_data(
        # asym_file=wd+"asymmetry.declan.na_filtered.tsv",
        asym_file=wd+"asymmetry.declan.na_filtered.renamed.tsv",
        singleton_file=wd+"singleton_counts.mapped.tsv",
        geno_counts=wd+"genotype_counts.mapped.tsv",
        nuc_file=wd+"capture_kit_info/xgen-exome-hyb-panel-v2/stranded/xgen.merged.nuc_file",
        ancestry_file=wd+"ancestry.mapped.tsv"
        )
    return data


if __name__ == "__main__":
    asymmetry_data = read_default_asym_data()
