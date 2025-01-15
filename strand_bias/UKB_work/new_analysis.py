
from itertools import product
from pathlib import Path

import pandas as pd
from pandas import IndexSlice as idx
from plotly.subplots import make_subplots
from scipy import stats

from plot_help import VarCounts, PlotDirectory
from hexbin import plot_hexbin
from ukb_mismatch_table import UkbMismatchTable
from ukb_variant_table import UkbVariantsUnstacked
from strand_bias.TCGA_analysis.capture_kit_counts import CaptureKit

# %% Read data.
data_dir = Path("~/StrandBias/UKB_analysis/")
var_dir = data_dir/"variant_reanalysis"

mismatches = UkbMismatchTable.read_csv(data_dir/"asymmetry.tsv")
oxog_variants = VarCounts(
    filtered=UkbVariantsUnstacked.read_csv(var_dir/"oxog_vars.snp_filtered.stats.project_ids.tsv"),
    unfiltered=UkbVariantsUnstacked.read_csv(var_dir/"oxog_vars.unfiltered.stats.project_ids.tsv")
    )
xgen = CaptureKit.read_capture_kit_nucleotide_summary(
    "/home/tyler/Documents/Resource_Data/capture_kits/IDT_xGen_Exome_Hyb_Panel/nt_counts.tsv"
    )

# %% Plot everything.

plots = PlotDirectory(mismatches, oxog_variants,
                      normalization_data=xgen,
                      mismatch_normalization_factor=1000,
                      variant_normalization_factor=1000000,
                      ratio_by_template=True)
# plots2 = PlotDirectory(mismatches, oxog_variants)


# %% Correlation testing.

def test_correlation(mismatch_table, variant_tables: VarCounts, bias_type):
    mismatch_data = mismatch_table.ukb_mismatches.calculate_bias(bias_type, True, False)[["ratio"]].sort_index()
    cols = product(["spearman", "pearson"], ["stat", "pvalue"])
    stat_table = pd.DataFrame(index=list(variant_tables._fields),
                              columns=pd.MultiIndex.from_tuples(cols))
    for name, variant_table in variant_tables._asdict().items():
        variant_data = variant_table.ukb_variants.calculate_bias(bias_type, True).Hets.sort_index()
        for test_name, test in [("spearman", stats.spearmanr), ("pearson", stats.pearsonr)]:
            result = test(mismatch_data.ratio, variant_data.ratio)
            stat_table.loc[name, (test_name, "stat")] = result.statistic
            stat_table.loc[name, (test_name, "pvalue")] = result.pvalue
    return stat_table


def test_subset_correlation(mismatch_table, variant_tables: VarCounts, bias_type,
                            normalization_data=None, ratio_by_template=True):
    data = mismatch_table.ukb_mismatches.calculate_bias(bias_type, True, False,
                                                           normalization_data,
                                                           ratio_by_template=ratio_by_template)
    var_data = variant_tables.filtered.ukb_variants.calculate_bias(bias_type, True,
                                                                   normalization_data,
                                                                   ratio_by_template=ratio_by_template)
    data.sort_values("ratio", inplace=True)
    cutoff = data.shape[0]//10
    upper_samples = list(data.index.unique("sample_id"))[-cutoff:]
    lower_samples = list(data.index.unique("sample_id"))[:cutoff]
    upper_samples = data.loc[idx[:, upper_samples],]
    lower_samples = data.loc[idx[:, lower_samples],]
    for group, table in (("all", data), ("upper", upper_samples), ("lower", lower_samples)):
        combo = table[["ratio"]].join(var_data.Hets, lsuffix="_m", rsuffix="_v")
        for test in (stats.spearmanr, stats.pearsonr):
            print(group, test(combo.ratio_m, combo.ratio_v))
    # return upper_samples, lower_samples


reference_bias_correlation = test_correlation(mismatches, oxog_variants, "reference")
transcription_bias_correlation = test_correlation(mismatches, oxog_variants, "transcription")

# %% Singleton Check.

counts = pd.read_csv("/home/tyler/Documents/Projects/StrandBias/UKB_analysis/singleton_checking/top10singletons"
                     ".unfiltered.counts.tsv", sep="\t", index_col=list(range(5)), header=0)

# %% Paper figures

def make_sfig_7(variant_table):
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=[f"<b>{label+')':<120}</b>" for label in "ab"])
    for i, bias_type in enumerate(["reference", "transcription"]):
        sub_objs = variant_table.ukb_variants._plot_bias_hexbin(
            bias_type=bias_type,
            normalization_data=xgen,
            normalization_factor=1000000,
            ratio_by_template=True,
            showlegend=bool(i),
            return_graph_objects=True
            )
        fig.layout[f"xaxis{i+1}"]["title"]["text"] = "C➔A variant calls"
        fig.layout[f"yaxis{i+1}"]["title"]["text"] = "G➔T variant calls"
        for datum in sub_objs:
            fig.add_trace(datum, row=1, col=i+1)
    return fig


fig5 = (mismatches.ukb_mismatches
        .plot_combined_bias_figure(xgen, 1000, True, manual_yx_line=[0, 30],
                                   highlight_top10=True, add_subplot_titles=False)
        .update_layout(title=None, height=1500, width=1500, legend={'itemsizing': 'constant'})
        )

fig7 = make_sfig_7(oxog_variants.filtered).update_layout(width=1500, height=750, font=dict(size=20))

