
import numpy as np
import pandas as pd
from pandas import IndexSlice as idx
import plotly.graph_objects as go
from scipy.stats import spearmanr

from import_data import (
    read_default_asym_data,
    read_reference_nuc_file,
    read_vcounts,
    )
from plotting import plot_hexbin, plot_per_variant_ratios, plot_by_batch_with_trendline
import aggregate
from aggregate import VARIANTS


# %% Data preparation.
# Read asym data.
print("Reading data...")
asym_data = read_default_asym_data()
nuc_counts = read_reference_nuc_file("./StrandBias/capture_kit_info/"
                                     "xgen-exome-hyb-panel-v2/stranded/xgen.merged.nuc_file")
# Read genotype counts.
vcounts = read_vcounts("./StrandBias/genotype_counts/vcount_master.tsv",
                       master_file=True)

# Subset samples to only Brits and only those in the asym data file.
print("Organizing data...")
asym_data = asym_data.loc[(asym_data.British == True) & (asym_data.ancestry_codes == "1001")]
vcounts = vcounts.loc[idx[list(asym_data.index), :, :]]

# Pull out variant counts, listed as variant with respect to strand, i.e. forward coding
# strand variants are reference direction nucleotides, while reverse coding strand
# variants are complementary to the reference nucleotides.
# Note: Table is multi-indexed with 4 levels - sample, coding strand, variant orientation,
# and variant.
vcounts_oriented = (vcounts.loc[idx[:, "forward", "reference", :]]
                    + vcounts.loc[idx[:, "reverse", "complement", :]])

# Calculate asymmetry of genotype calls.
# Ratio of G nucleotides in forward vs. reverse coding regions:
g_ratio = nuc_counts.G_count[0] / nuc_counts.G_count[1]
# Non-captured-region variant counts:
vcounts_rev_oriented = (vcounts.loc[idx[:, "forward", "complement", :]]
                        + vcounts.loc[idx[:, "reverse", "reference", :]])
# Ratio of captured-region variants vs. non-captured-region variants:
vcount_ratios = vcounts_oriented / vcounts_rev_oriented
# Ratio normalized by the number of Gs in each region:
# FIX: THIS IS WRONG. FIX ME.
vcount_ratios_g_normalized = vcount_ratios.loc[idx[:, "GtoT"], ] / g_ratio


# %% Colinearity tests: Asymmetry vs. genotypes and singletons.
print("Running colinearity tests...")
cols = ["nHets", "nSingletons", "nNonRefHom", "nNonRefAlleles"]
spearman_results = dict(zip(cols, [dict(), dict(), dict(), dict()]))

for variant in VARIANTS:
    data = vcounts_oriented.loc[idx[:, variant], ]
    for col in cols:
        spearman_results[col][variant] = spearmanr(asym_data.Asymmetry_norm, data[col])
spearman_results = pd.DataFrame.from_dict(spearman_results)
spearman_results = spearman_results.applymap(lambda x: x.statistic if x.pvalue < 0.05
                                             else np.nan)

# %% Colinearity tests: Asymmetry batch-adjusted vs. genotypes and singletons.
print("Running batch-adjusted colinearity tests...")
spearman_results_adjusted = dict(zip(cols, [dict(), dict(), dict(), dict()]))

for variant in VARIANTS:
    data = vcounts_oriented.loc[idx[:, variant], ]
    for col in cols:
        spearman_results_adjusted[col][variant] = spearmanr(asym_data.asym_batch_adjusted, data[col])
spearman_results_adjusted = pd.DataFrame.from_dict(spearman_results_adjusted)
spearman_results_adjusted = spearman_results_adjusted.applymap(
    lambda x: x.statistic if x.pvalue < 0.05
    else np.nan
    )


# %% Filtered colinearity tests.
print("Filtering results...")
asym_data_filtered = asym_data.loc[asym_data.Asymmetry_norm_log2 > 1]
vcounts_oriented_filtered = vcounts_oriented.loc[asym_data_filtered.index]
spearman_results_filtered = dict(zip(cols, [dict(), dict(), dict(), dict()]))

for variant in VARIANTS:
    data = vcounts_oriented_filtered.loc[idx[:, variant], ]
    for col in cols:
        spearman_results_filtered[col][variant] = spearmanr(asym_data_filtered.Asymmetry_norm,
                                                            data[col])
spearman_results_filtered = pd.DataFrame.from_dict(spearman_results_filtered)
spearman_results_filtered = spearman_results_filtered.applymap(lambda x: x.statistic if x.pvalue < 0.05
                                                               else np.nan)

# %% Batch analysis.
print("Running colinearity tests by batch...")
batches = set(asym_data.batch)
spearmanr_batches = dict()
for batch in batches:
    filter_mask = asym_data.batch == batch
    spearmanr_batches[batch] = spearmanr(
        asym_data.Asymmetry_norm_log2[filter_mask],
        asym_data.Burden_per_1k_reads_log2[filter_mask]
        )
spearmanr_batches = {key: (val.statistic, val.pvalue)
                     for key, val in spearmanr_batches.items()}
spearmanr_batches_table = pd.DataFrame.from_dict(spearmanr_batches, orient="index")
spearmanr_batches_table.columns = ["statistic", "pvalue"]


# %% Hexbin plots.
print("Generating plots...")
hexbin_asym_vs_singletons = plot_hexbin(asym_data.Asymmetry_norm, asym_data.HET_SNP_CT, "inferno")
hexbin_asym_vs_singletons_bed = plot_hexbin(
    asym_data.Asymmetry_norm,
    vcounts.loc[idx[:, :, "reference", :]].groupby(["sample_id"]).agg(sum).nSingletons,
    "inferno"
    )
hexbin_ca_ratio_vs_genotypes = plot_hexbin(asym_data.Asymmetry_norm, asym_data.nHets, "inferno")

# batch_plot = plot_asym_vs_burden_by_batch(asym_data)
