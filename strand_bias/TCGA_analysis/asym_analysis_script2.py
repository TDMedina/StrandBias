
from collections import Counter, defaultdict
from datetime import datetime
from itertools import product
from operator import invert

from pybedtools import BedTool
import pandas as pd
from pandas import IndexSlice as idx
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import numpy as np
from numpy import inf
from scipy import stats

from strand_bias.pileup_parser.pileup_table import ConcatenatedPileupTable
from strand_bias.pileup_parser.call_table import CallTallyTable
from strand_bias.TCGA_analysis.plotting import plot_single_call_asymmetry_box, plot_heatmap_of_unique_site_counts
from strand_bias.pileup_parser.aggregate_asym_tables import read_capture_kit_nucleotide_summary
from strand_bias.misc import gdc_api

from capture_kit_counts import CaptureKit
from utilities import (
    _CHANGES,
    _COMPS,
    _normalize_minimum,
    make_complement,
    reverse_complement,
    filter_by_fences,
    filter_raw_call_data_dups
    )
from plotting import (
    plot_mismatch_asymmetry,
    # plot_mismatch_asym_per_change_per_project,
    plot_call_asymmetry,
    plot_oxog_calls_vs_mismatches_per_project
    )
pio.renderers.default = "browser"


# %% Read data.

# VCROME = read_capture_kit_nucleotide_summary("/home/tyler/Documents/Projects/StrandBias/VCRome.hg38.nt_counts.tsv")
VCROME = CaptureKit.read_capture_kit_nucleotide_summary("/home/tyler/Documents/Projects/StrandBias/VCRome.hg38.nt_counts.tsv")
mismatch_data = ConcatenatedPileupTable.read_csv("/home/tyler/StrandBias/Analysis/concatenated_asym_tables.rename.tsv")
call_data = CallTallyTable.read_csv("/home/tyler/StrandBias/VCF_analysis/variant_filter_tally3.tsv")
frequencies = pd.read_csv("/home/tyler/Documents/Projects/StrandBias/VCF_analysis/variant_frequencies.tsv",
                          sep="\t", index_col=list(range(8)))
call_data = call_data.join(frequencies)

# call_data = CallTable.read_csv("/home/tyler/StrandBias/Analysis/vcf_counts.tsv")
# call_reference_asymmetry = call_data.call_tools.simplify_reference_asymmetry()

_proj_numbers = {y: x for x, y in enumerate(mismatch_data.index.levels[0], start=-4)}


# %% Make mismatch and call reference asymmetry plots.

mismatch_ref_plots = plot_mismatch_asymmetry(
    bias_type="reference",
    mismatch_table=mismatch_data,
    log_transform_ratio=False,
    normalize_by_nt_content=True,
    normalization_data=VCROME,
    normalization_factor=1000,
    shared_yscale="box")

mismatch_trans_plots = plot_mismatch_asymmetry(
    bias_type="transcription",
    mismatch_table=mismatch_data,
    log_transform_ratio=False,
    normalize_by_nt_content=True,
    normalization_data=VCROME,
    normalization_factor=1000,
    shared_yscale="box")

call_ref_plots = plot_call_asymmetry(
    call_table=call_data,
    bias_type="reference",
    filter_outliers=True,
    normalize_by_nt_content=True,
    normalization_data=VCROME,
    # normalization_factor=1000000,
    log_transform_ratio=False,
    shared_yscale="box")

call_trans_plots = plot_call_asymmetry(
    call_table=call_data,
    bias_type="transcription",
    filter_outliers=True,
    normalize_by_nt_content=True,
    normalization_data=VCROME,
    # normalization_factor=1000000,
    log_transform_ratio=False,
    shared_yscale="box")

call_ref_gt_plot = plot_single_call_asymmetry_box(
    call_table=call_data,
    bias_type="reference",
    numerator="GT",
    denominator="CA",
    filter_outliers=True,
    iqr_factor=1.5,
    normalize_by_nt_content=True,
    normalization_factor=1000000,
    normalization_data=VCROME,
    )

call_trans_gt_plot = plot_single_call_asymmetry_box(
    call_table=call_data,
    bias_type="transcription",
    numerator="GT",
    denominator="CA",
    filter_outliers=True,
    iqr_factor=1.5,
    normalize_by_nt_content=True,
    normalization_factor=1000000,
    normalization_data=VCROME,
    )

# %% Binomial tests comparing incidence of each nt change pair.

def calculate_binomtest_for_mismatches(mismatch_table, filter_outliers=False, iqr_factor=1.5,
                                       normalize_by_nt_content=False, normalization_factor=1,
                                       minimum_raw_count=None, changes=None):
    normalization_counts = VCROME.sum(axis=0) if normalize_by_nt_content else None
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    if changes is None:
        changes = _CHANGES
    binom_results = dict()

    # Binomial test per nucleotide change pair.
    for i, change in enumerate(changes):
        binom_change = dict()

        comp = make_complement(change)
        change_comparison_data = (mismatch_table.mismatch_tools.calculate_reference_asymmetry(
            change_numerator=change,
            change_denominator=comp,
            add_ratio_column=False,
            normalization_counts=normalization_counts,
            normalization_factor=normalization_factor
            ))

        # Filtering for total binomial test.
        if filter_outliers:
            filtered = filter_by_fences(change_comparison_data, [change, comp], iqr_factor)
        elif minimum_raw_count is not None:
            filtered = change_comparison_data.loc[(minimums[change[0]] <= change_comparison_data[change])
                                                  & (minimums[comp[0]] <= change_comparison_data[comp])]
        else:
            filtered = change_comparison_data

        # Binomial test on totals across TCGA.
        binom_change["total"] = stats.binomtest(
            k=int(filtered[change].sum()),
            n=int(filtered[[change, comp]].sum().sum())
            ).pvalue

        # Binomial tests per project.
        for proj in _proj_numbers.keys():

            # Filtering per project.
            filtered = change_comparison_data.loc[idx[proj, :],]
            if filter_outliers:
                filtered = filter_by_fences(filtered, [change, comp], iqr_factor)
            elif minimum_raw_count is not None:
                filtered = filtered.loc[(minimums[change[0]] <= filtered[change])
                                        & (minimums[comp[0]] <= filtered[comp])]

            # Binomial test for project.
            binom_change[proj] = stats.binomtest(
                k=int(filtered[change].sum()),
                n=int(filtered[[change, comp]].sum().sum())
                ).pvalue

        binom_results[change] = binom_change
    binom_results = pd.DataFrame(binom_results)
    return binom_results


def calculate_binomtest_for_calls(call_table, filter_outliers=False, iqr_factor=1.5,
                                  normalize_by_nt_content=False, normalization_factor=1,
                                  minimum_raw_count=None, changes=None):
    normalization_counts = VCROME.sum(axis=0) if normalize_by_nt_content else None
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    if changes is None:
        changes = _CHANGES
    binom_results = dict()

    # Binomial test per nucleotide change pair.
    for i, change in enumerate(changes):
        for filter_status in ["TOTAL", "PASS", "FAIL"]:
            binom_change = dict()

            comp = make_complement(change)
            change_comparison_data = (
                call_table.call_tools.calculate_reference_asymmetry(change_numerator=change,
                                                                    change_denominator=comp,
                                                                    add_ratio_column=False,
                                                                    normalization_counts=normalization_counts,
                                                                    normalization_factor=normalization_factor)
                )[filter_status]

            # Filtering for total binomial test.
            if filter_outliers:
                filtered = filter_by_fences(change_comparison_data, [change, comp], iqr_factor)
            elif minimum_raw_count is not None:
                filtered = change_comparison_data.loc[(minimums[change[0]] <= change_comparison_data[change])
                                                      & (minimums[comp[0]] <= change_comparison_data[comp])]
            else:
                filtered = change_comparison_data

            # Binomial test on totals across TCGA.
            binom_change["total"] = stats.binomtest(
                k=int(filtered[change].sum()),
                n=int(filtered[[change, comp]].sum().sum())
                ).pvalue

            # Binomial tests per project.
            for proj in _proj_numbers.keys():

                # Filtering per project.
                filtered = change_comparison_data.loc[idx[proj, :],]
                if filter_outliers:
                    filtered = filter_by_fences(filtered, [change, comp], iqr_factor)
                elif minimum_raw_count is not None:
                    filtered = filtered.loc[(minimums[change[0]] <= filtered[change])
                                            & (minimums[comp[0]] <= filtered[comp])]

                # Binomial test for project.
                binom_change[proj] = stats.binomtest(
                    k=int(filtered[change].sum()),
                    n=int(filtered[[change, comp]].sum().sum())
                    ).pvalue

            binom_results[(filter_status, change)] = binom_change
    binom_results = pd.DataFrame(binom_results)
    binom_results.columns = pd.MultiIndex.from_tuples(binom_results.columns)
    binom_results = binom_results[sorted(binom_results.columns)]
    return binom_results


mismatch_binom_results_min = calculate_binomtest_for_mismatches(
    mismatch_table=mismatch_data,
    minimum_raw_count=100
    )

# mismatch_binom_results = calculate_binomtest_for_mismatches(
#     mismatch_table=mismatch_data,
#     )

call_binom_results = calculate_binomtest_for_calls(
    call_table=call_data
    )


# %% Calculate reference asymmetry mismatch vs. call correlation.

def filter_call_dups(call_change_table):
    id_counts = Counter([x[1] for x in call_change_table.index.to_list()])
    dups = [x for x, y in id_counts.items() if y > 1]
    table = (call_change_table
             .reset_index(level=["project_id", "file_id"])
             .drop(dups)
             .reset_index().set_index(["project_id", "case_id", "file_id"]))
    return table


def make_corr_table(
        mismatch_table, call_table, bias_type,
        filter_outliers=False, iqr_factor=1.5,
        normalization_data=None,
        # log_transform_ratio=False,
        minimum_mismatch_count=None,
        corr_test="spearman",
        filters=None):
    if filters is None:
        filters = ["TOTAL", "PASS"]
    corr_test = {"spearman": stats.spearmanr,
                 "pearson": stats.pearsonr}[corr_test]
    if normalization_data is not None:
        normalization_data = normalization_data.capkit.calculate_counts(bias_type)
    minimums = _normalize_minimum(minimum_mismatch_count, True, 1, normalization_data)

    corr_results = defaultdict(dict)

    call_table = filter_raw_call_data_dups(call_table)

    # Correlation test per nucleotide change pair.
    for i, change in enumerate(_CHANGES):

        comp = make_complement(change)
        ratio = f"{change}{comp}_ratio".lower()
        mismatch_subset = mismatch_table.mismatch_tools.calculate_asymmetry(
            bias_type=bias_type,
            change_numerator=change,
            change_denominator=comp,
            add_ratio_column=True,
            log_transform_ratio=False,
            add_fraction_column=False,
            normalization_counts=normalization_data,
            normalization_factor=1).droplevel("file_id")

        call_subset = call_table.call_tools.calculate_asymmetry(
            bias_type=bias_type,
            change_numerator=change,
            change_denominator=comp,
            add_ratio_column=True,
            log_transform_ratio=False,
            add_fraction_column=False,
            normalization_counts=normalization_data,
            normalization_factor=1)

        for filter_status in filters:
            combo = call_subset[filter_status].join(mismatch_subset,
                                                    lsuffix="_calls",
                                                    rsuffix="_mismatches")
            if minimum_mismatch_count is not None:
                combo = combo.loc[(minimums[change[0]] <= combo[f"{change}_mismatches"])
                                  & (minimums[comp[0]] <= combo[f"{comp}_mismatches"])]
            for field in [f"{ratio}_calls", f"{ratio}_mismatches"]:
                combo = combo.loc[~ pd.isna(combo[field])]
                combo = combo.loc[combo[field].abs() != inf]
            # Binomial test on totals across TCGA.
            result = corr_test(combo[f"{ratio}_calls"],
                               combo[f"{ratio}_mismatches"])
            corr_results["total"][(filter_status, change, "statistic")] = result.statistic
            corr_results["total"][(filter_status, change, "pvalue")] = result.pvalue

            # Binomial tests per project.
            for proj in _proj_numbers.keys():
                # Filtering per project.
                proj_combo = combo.loc[idx[proj, :],]

                # Binomial test for project.
                result = corr_test(proj_combo[f"{ratio}_calls"],
                                   proj_combo[f"{ratio}_mismatches"])
                corr_results[proj][(filter_status, change, "statistic")] = result.statistic
                corr_results[proj][(filter_status, change, "pvalue")] = result.pvalue

    corr_results = pd.DataFrame.from_dict(corr_results, orient="index")
    corr_results = corr_results[sorted(corr_results.columns)]
    return corr_results


def calculate_call_mismatch_correlation(call_table, mismatch_table,
                                        minimum_mismatch_count=None,
                                        corr_test="spearman"):
    corr_test = {"spearman": stats.spearmanr,
                 "pearson": stats.pearsonr}[corr_test]
    minimums = _normalize_minimum(minimum_mismatch_count)

    corr_results = defaultdict(dict)

    call_table = filter_raw_call_data_dups(call_table)

    # Correlation test per nucleotide change pair.
    for i, change in enumerate(_CHANGES):

        comp = make_complement(change)
        ratio = f"{change}{comp}_ratio".lower()
        mismatch_subset = (mismatch_table.mismatch_tools.calculate_reference_asymmetry(change, comp)
                           .droplevel("file_id"))

        call_subset = call_table.call_tools.calculate_reference_asymmetry(change, comp)
        # call_subset = filter_call_dups(call_subset).droplevel("file_id")

        for filter_status in ["TOTAL", "PASS", "FAIL"]:
            combo = call_subset[filter_status].join(mismatch_subset,
                                                    lsuffix="_calls",
                                                    rsuffix="_mismatches")
            if minimum_mismatch_count is not None:
                combo = combo.loc[(minimums[change[0]] <= combo[f"{change}_mismatches"])
                                  & (minimums[comp[0]] <= combo[f"{comp}_mismatches"])]
            for field in [f"{ratio}_calls", f"{ratio}_mismatches"]:
                combo = combo.loc[~ pd.isna(combo[field])]
                combo = combo.loc[combo[field].abs() != inf]
            # Binomial test on totals across TCGA.
            result = corr_test(combo[f"{ratio}_calls"],
                               combo[f"{ratio}_mismatches"])
            corr_results["total"][(filter_status, change, "statistic")] = result.statistic
            corr_results["total"][(filter_status, change, "pvalue")] = result.pvalue

            # Binomial tests per project.
            for proj in _proj_numbers.keys():

                # Filtering per project.
                proj_combo = combo.loc[idx[proj, :],]

                # Binomial test for project.
                result = corr_test(proj_combo[f"{ratio}_calls"],
                                   proj_combo[f"{ratio}_mismatches"])
                corr_results[proj][(filter_status, change, "statistic")] = result.statistic
                corr_results[proj][(filter_status, change, "pvalue")] = result.pvalue

    corr_results = pd.DataFrame.from_dict(corr_results, orient="index")
    corr_results = corr_results[sorted(corr_results.columns)]
    return corr_results


corr = calculate_call_mismatch_correlation(call_data, mismatch_data)
corr_pearson = calculate_call_mismatch_correlation(call_data, mismatch_data,
                                                   corr_test="pearson")


# %% Correlation by allele filtering.

def calculate_af_correlation(mismatch_table, call_table, projects=None):
    if projects is not None:
        if isinstance(projects, str):
            projects = [projects]
    else:
        projects = ["total"] + list(_proj_numbers.keys())
    results = dict()
    afs = [round(n, 3) for n in np.linspace(0, 1, 101)]
    dedup_calls = filter_raw_call_data_dups(call_table)
    steps = len(_CHANGES)*len(afs)*len(projects)*2
    i = 1
    for change, comp in zip(_CHANGES, _COMPS):
        mismatches = mismatch_table
        for pass_fail in ["TOTAL", "PASS"]:
            for proj in projects:
                for af in afs:
                    print(f"Step {i}/{steps}.\r", end="")
                    calls = (dedup_calls.loc[dedup_calls.frequency >= af])[pass_fail]
                    proj_slice = idx[proj if not proj == "total" else slice(None), :, :]
                    calls = calls.loc[proj_slice].droplevel("file_id")
                    combo = (calls
                             .join(mismatches.droplevel("file_id"),
                                   lsuffix="_calls", rsuffix="_mismatches"))
                    fields = [f"{change}{comp}_ratio_mismatches".lower(),
                              f"{change}{comp}_ratio_calls".lower()]
                    for field in fields:
                        combo = combo.loc[~ pd.isna(combo[field])]
                        combo = combo.loc[combo[field].abs() != inf]
                    corr_test = stats.spearmanr(combo[fields[0]],
                                                combo[fields[1]])
                    label = (proj, pass_fail, change, comp, af)
                    results[label] = [corr_test.statistic, corr_test.pvalue, combo.shape[0]]
                    i += 1
    results = pd.DataFrame.from_dict(results, orient="index")
    results.index = pd.MultiIndex.from_tuples(results.index)
    results.index.names = ["project_id", "filtered", "change", "comp", "af"]
    results.columns = ["stat", "pvalue", "sample_count"]
    return results


def calculate_depth_correlation(mismatch_table, call_table, projects=None):
    if projects is not None:
        if isinstance(projects, str):
            projects = [projects]
    else:
        projects = ["total"] + list(_proj_numbers.keys())
    results = dict()
    depths = list(range(41))
    depths = [2**n for n in range(7)]
    dedup_calls = filter_raw_call_data_dups(call_table)
    steps = len(_CHANGES)*len(depths)*len(projects)*2
    i = 1
    for change, comp in zip(_CHANGES, _COMPS):
        mismatches = mismatch_table
        for pass_fail in ["TOTAL", "PASS"]:
            for proj in projects:
                for depth in depths:
                    print(f"Step {i}/{steps}.\r", end="")
                    try:
                        calls = (dedup_calls.loc[dedup_calls.alt_depth >= depth])[pass_fail]
                        proj_slice = idx[proj if not proj == "total" else slice(None), :, :]
                        calls = calls.loc[proj_slice].droplevel("file_id")
                    except KeyError:
                        continue
                    combo = (calls
                             .join(mismatches.droplevel("file_id"),
                                   lsuffix="_calls", rsuffix="_mismatches"))
                    fields = [f"{change}{comp}_ratio_mismatches".lower(),
                              f"{change}{comp}_ratio_calls".lower()]
                    for field in fields:
                        combo = combo.loc[~ pd.isna(combo[field])]
                        combo = combo.loc[combo[field].abs() != inf]
                    corr_test = stats.spearmanr(combo[fields[0]],
                                                combo[fields[1]])
                    label = (proj, pass_fail, change, comp, depth)
                    results[label] = [corr_test.statistic, corr_test.pvalue, combo.shape[0]]
                    i += 1
    results = pd.DataFrame.from_dict(results, orient="index")
    results.index = pd.MultiIndex.from_tuples(results.index)
    results.index.names = ["project_id", "filtered", "change", "comp", "depth"]
    results.columns = ["stat", "pvalue", "sample_count"]
    return results


# %% Binomial test of reference call asymmetry in aggregate, counting each site once across
# all samples.
def calculate_binomtest_for_unique_calls(call_table, changes=None):
    if changes is None:
        changes = _CHANGES
    filters = ["TOTAL", "PASS", "FAIL"]
    binom_results = dict()

    for proj in list(_proj_numbers) + [slice(None)]:
        proj_dict = dict()
        uniques = ((call_table.loc[proj].groupby(["contig", "pos", "ref", "alt"]).sum() > 0)
                   .astype(int)
                   .groupby(["ref", "alt"])
                   .sum())
        for pass_fail in filters:
            filtered = uniques[pass_fail]
            for change in changes:
                comp = make_complement(change)
                change_count = filtered.loc[idx[tuple(change)]]
                comp_count = filtered.loc[idx[tuple(comp)]]
                test = stats.binomtest(change_count, change_count + comp_count).pvalue
                proj_dict[(pass_fail, change)] = test
        binom_results[proj if isinstance(proj, str) else "total"] = proj_dict
    binom_results = pd.DataFrame.from_dict(binom_results, orient="index")
    return binom_results


def calculate_binomtest_for_mismatches(mismatch_table, changes=None, p=None):
    if changes is None:
        changes = _CHANGES
    if p is None:
        p = {change: 0.5 for change in changes}
    binom_results = dict()

    for change in changes:
        comp = make_complement(change)
        subset = (mismatch_table.mismatch_tools
                  .calculate_reference_asymmetry(change, comp,
                                                 add_ratio_column=False,
                                                 add_fraction_column=False))
        subset = subset.groupby("project_id").sum()
        subset = subset.astype(int)
        change_dict = dict()
        for proj in list(_proj_numbers):
            test = stats.binomtest(subset.loc[proj, change],
                                   subset.loc[proj, change]+subset.loc[proj, comp],
                                   p=p[change])
            change_dict[proj if isinstance(proj, str) else "total"] = test.pvalue
        change_dict["total"] = stats.binomtest(subset[change].sum(),
                                               subset.sum().sum(), p=p[change]).pvalue
        binom_results[change] = change_dict
    binom_results = pd.DataFrame(binom_results)
    return binom_results


unique_binom = calculate_binomtest_for_unique_calls(call_data)


# %% Count unique PASS sites per project.
all_uniques = ((call_data
                .groupby(["project_id", "contig", "pos", "ref", "alt"])
                .sum() > 0).astype(int))["PASS"]
all_unique_counts = all_uniques.loc[all_uniques > 0].groupby(["project_id", "ref", "alt"]).sum()


def plot_aggregate_oxog_frequencies(call_table):
    all_uniques = ((call_table
                    .groupby(["project_id", "contig", "pos", "ref", "alt"])
                    .sum() > 0).astype(int))["PASS"]
    all_unique_counts = all_uniques.loc[all_uniques > 0].groupby(["project_id", "ref", "alt"]).sum()
    gt = all_unique_counts.loc[idx[:, "G", "T"]].reset_index()
    gt["change"] = "GT"
    ca = all_unique_counts.loc[idx[:, "C", "A"]].reset_index()
    ca["change"] = "CA"
    plot_data = pd.concat([gt, ca])
    plot_data = plot_data.set_index("project_id").sort_index()
    for proj in _proj_numbers:
        plot_data.loc[proj, "PASS"] = plot_data.loc[proj, "PASS"] / plot_data.loc[proj, "PASS"].sum()
    fig = px.bar(plot_data.reset_index(), x="project_id", y="PASS", color="change")
    return fig


# %% Check variant contexts of all oxo-G variants.

def count_oxog_contexts(call_table, ref_genome):
    all_uniques = ((call_table
                    .groupby(["project_id", "contig", "pos", "ref", "alt"])
                    .sum() > 0).astype(int))["PASS"]

    all_contexts = defaultdict(dict)
    for proj in _proj_numbers:
        for change in ["GT", "CA"]:
            strings = [f"{x.contig}\t{x.pos - 2}\t{x.pos+1}"
                       for x in (all_uniques
                                 .loc[all_uniques > 0]
                                 .loc[idx[proj, :, :, change[0], change[1]]]
                                 .reset_index()
                                 .itertuples())]
            bedtool = BedTool("\n".join(strings), True)
            bedtool.sequence(fi=ref_genome)
            with open(bedtool.seqfn) as infile:
                context_results = infile.readlines()
            context_results = [x.rstrip() for x in context_results if not x.startswith(">")]
            all_contexts[proj][change] = Counter(context_results)
    return all_contexts


def make_oxog_context_table(oxog_context_counts, complement_C_contexts=True):
    table = pd.DataFrame()
    # contexts = product("ACGT", "CG", "ACGT")
    for cohort, change_dict in oxog_context_counts.items():
        for change, counter in change_dict.items():
            for context, count in counter.items():
                table.loc[idx[cohort], context] = count
    cols = [(context[1], context) for context in table.columns]
    if complement_C_contexts:
        for i, (nt, context) in enumerate(cols):
            if nt == "C":
                cols[i] = (nt, reverse_complement(context))
    table.columns = pd.MultiIndex.from_tuples(cols, names=["base", "context"])
    table = table[sorted(table.columns)]
    return table

oxog_contexts = count_oxog_contexts(call_data, "/home/tyler/Documents/Resource_Data/reference_genomes/GRCh38.d1.vd1.fa")


# %% Binomial simulation of TGCT.

def plot_tgct_oxog_simulation(call_table, iterations):
    prob = VCROME.G.sum() / (VCROME.G.sum() + VCROME.C.sum())
    tgct_oxog = (call_table.call_tools.calculate_reference_asymmetry("GT", "CA")
                 .loc[idx["TCGA-TGCT", :, :], "PASS"])
    tgct_oxog["total"] = (tgct_oxog.GT + tgct_oxog.CA).astype(int)
    total_g = tgct_oxog.GT.sum()
    simulations = [stats.binom.rvs(tgct_oxog.total, prob).sum()
                   for _ in range(iterations)]
    fig = make_subplots(1, 2)
    fig.add_trace(go.Histogram(x=simulations, histnorm="probability density"))
    y_limit = max(Counter(simulations).values()) / iterations
    fig.add_trace(go.Scatter(x=[total_g]*2, y=[0, y_limit], mode="lines"))
    fig.add_trace(go.Box(y=tgct_oxog.GT / tgct_oxog.total, name="TCGA-TGCT"), row=1, col=2)
    fig.add_trace(go.Box(y=stats.binom.rvs(tgct_oxog.total, prob) / tgct_oxog.total, name="Simulated Example"),
                  row=1, col=2)
    fig.update_layout(showlegend=False)
    gap = (max(simulations) - min(simulations)) / 5
    fig.layout.xaxis.range = [min(simulations) - gap, total_g + gap]
    return fig


# %% Read GDC file data.

bam_data = gdc_api.retrieve_all_data()
bam_data = gdc_api.filter_for_paired_sequences(bam_data, true_pairs=True)
bam_data = gdc_api.filter_by_capture_kit_number(bam_data, "06 465 668 001", include_wgs=False)
bam_data = bam_data.reset_index().set_index([("cases", "project.project_id"), "case_id", "file_id"])
bam_data.index.names = ["project_id", "case_id", "file_id"]

oxo_by_date = mismatch_data.mismatch_tools.calculate_reference_asymmetry("GT", "CA")
oxo_by_date = oxo_by_date.join(bam_data.analysis.sequencing_date)
oxo_by_date.sequencing_date = [datetime.strptime(x, "%Y-%m-%dT%H") if isinstance(x, str) else x
                               for x in oxo_by_date.sequencing_date]


def plot_seq_date(oxo_by_date_table):
    fig = px.scatter(oxo_by_date_table.reset_index(), x="sequencing_date", y="gtca_ratio", color="project_id")
    fig.update_layout(legend_title="TCGA Project:", xaxis_title="Sequencing Date", yaxis_title="G➔T/C➔A ratio")
    return fig

# oxo_date_plot = px.scatter(oxo_by_date.reset_index(), x="sequencing_date", y="gtca_ratio", color="project_id")


# %% Paper figures.
fig1 = mismatch_ref_plots[2].update_layout(title=None, height=1000, width=1000)
fig2 = plot_call_asymmetry(call_table=call_data, bias_type="reference", filter_outliers=True,
                           normalize_by_nt_content=True, normalization_data=VCROME,
                           shared_yscale="box", filters=["TOTAL"])[1].update_layout(title=None, height=1000, width=1000)
fig3 = call_ref_gt_plot.update_layout(title=None, height=600, width=1000)
fig4 = plot_heatmap_of_unique_site_counts(call_data).update_layout(title=None, height=800, width=900)


# %% Supplement figures.
sfig1 = mismatch_ref_plots[0].update_layout(title=None, height=1000, width=1000)
sfig2 = plot_seq_date(oxo_by_date).update_layout(height=500, width=1000)
sfig3 = mismatch_trans_plots[0].update_layout(title=None, height=1000, width=1000)
# sfig4 is a WGS figure.
# sfig5 is a WGS figure.
sfig6 = call_trans_gt_plot.update_layout(title=None, height=600, width=1000)
