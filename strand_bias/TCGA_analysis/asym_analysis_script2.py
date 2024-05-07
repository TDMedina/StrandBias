
from collections import Counter, defaultdict
from datetime import datetime
from itertools import product, combinations, permutations

from pybedtools import BedTool
import pandas as pd
from pandas import IndexSlice as idx
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import numpy as np
from numpy import inf, log2
from scipy import stats

from strand_bias.pileup_parser.pileup_table import ConcatenatedPileupTable
from strand_bias.pileup_parser.call_table import CallTallyTable
from strand_bias.pileup_parser.aggregate_asym_tables import read_capture_kit_nucleotide_summary
from strand_bias.misc import gdc_api

pio.renderers.default = "browser"

_REVCOMP_DICT = {"A": "T", "C": "G"}
_REVCOMP_DICT |= {value: key for key, value in _REVCOMP_DICT.items()}

VCROME = read_capture_kit_nucleotide_summary("/home/tyler/Documents/Projects/StrandBias/VCRome.hg38.nt_counts.tsv")


# %% Functions.

# def plot_per_project_ref_asym(asym_table, x, y, title=""):
#     fig = px.scatter(asym_table.reset_index(),
#                      x=x, y=y, color="project_id",
#                      marginal_x="box", marginal_y="box")
#     fig.update_traces(notched=False, selector=dict(type="box"))
#
#     ratio_corner = px.box(asym_table.reset_index(), y="gtca_ratio", color="project_id")
#     fig.add_traces(ratio_corner.data, rows=2, cols=2)
#     fig.update_traces(showlegend=False, row=2, col=2)
#     fig.for_each_trace(
#         lambda trace: trace.update(x0=_proj_numbers[trace.name]) if trace.xaxis == "x4" else ())
#
#     fig.update_layout(dict(
#         title=title,
#         yaxis4_showticklabels=True,
#         yaxis4_dtick=1,
#         yaxis4_title="G➔T / C➔A Ratio"
#         ))
#     return fig


# %% Read data.

mismatch_data = ConcatenatedPileupTable.read_csv("/home/tyler/StrandBias/Analysis/concatenated_asym_tables.rename.tsv")
call_data = CallTallyTable.read_csv("/home/tyler/StrandBias/VCF_analysis/variant_filter_tally3.tsv")
frequencies = pd.read_csv("/home/tyler/Documents/Projects/StrandBias/VCF_analysis/variant_frequencies.tsv",
                          sep="\t", index_col=list(range(8)))
call_data = call_data.join(frequencies)

# call_data = CallTable.read_csv("/home/tyler/StrandBias/Analysis/vcf_counts.tsv")
# call_reference_asymmetry = call_data.call_tools.simplify_reference_asymmetry()

_proj_numbers = {y: x for x, y in enumerate(mismatch_data.index.levels[0], start=-4)}
CHANGES = ["TC", "CT", "TA", "TG", "CG", "GT"]


# %% Per-project reference asymmetry by nucleotide change.

def reverse_complement(seq):
    revcomp = "".join([_REVCOMP_DICT[x] for x in seq])
    return revcomp


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


def make_change_labels(change_numerators):
    change_labels = [[change, reverse_complement(change)] for change in change_numerators]
    change_labels = [["➔".join(list(change)) for change in changes]
                     for changes in change_labels]
    change_labels = ["<br>".join([f"  {changes[0]} /", changes[1]]) for changes in change_labels]
    return change_labels


def make_change_label(change_string):
    change_string = "➔".join(list(change_string))
    return change_string


SCALE_LABELS = {1: "", 1000: "1K ", 1000000: "1M"}


def make_scaled_label(change_label, scale_factor):
    label = f"{change_label} per {SCALE_LABELS[scale_factor]}{change_label[0]} nts"
    return label


def _normalize_minimum(minimum_raw_count, normalize_by_nt_content=False,
                       normalization_factor=1):
    if minimum_raw_count is None:
        minimums = {nt: 0 for nt in "ACGT"}
    elif normalize_by_nt_content:
        minimums = minimum_raw_count / VCROME.sum(axis=0) * normalization_factor
    else:
        minimums = {nt: minimum_raw_count for nt in "ACGT"}
    return minimums


def _make_axis_range(mins, maxes):
    axis_min, axis_max = min(mins), max(maxes)
    inter_range = axis_max - axis_min
    axis_range = [axis_min - inter_range/5, axis_max+inter_range/5]
    return axis_range


def plot_mismatch_ref_asym_per_change_per_project(mismatch_table, filter_outliers=False, iqr_factor=1.5,
                                                  normalize_by_nt_content=False, normalization_factor=1,
                                                  use_scaled_labels=True, log_transform_ratio=False,
                                                  shared_yscale=False,
                                                  minimum_raw_count=None):
    normalization_counts = VCROME.sum(axis=0) if normalize_by_nt_content else None
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    plot_box_ratio_by_nt_change = make_subplots(3, 2, start_cell="bottom-left")
    plot_scatter_by_nt_change = make_subplots(3, 2, start_cell="bottom-left")
    box_range = [[], []]
    scatter_range = [[], []]

    comps = [reverse_complement(change) for change in CHANGES]
    for i, (change, comp) in enumerate(zip(CHANGES, comps)):
        ratio = f"{change}{comp}_ratio".lower()

        change_comparison_data = (mismatch_table.mismatch_tools.simplify_for_reference_asymmetry(
            change_numerator=change,
            change_denominator=comp,
            log_transform_ratio=log_transform_ratio,
            normalization_counts=normalization_counts,
            normalization_factor=normalization_factor
            ))

        for j, proj in enumerate(_proj_numbers.keys()):

            # Pre-subset and treat data.
            proj_data = change_comparison_data.loc[idx[proj, :],]
            if filter_outliers:
                proj_data = filter_by_fences(proj_data, [change, comp], iqr_factor)
            if minimum_raw_count is not None:
                proj_data = proj_data.loc[(minimums[change[0]] <= proj_data[change])
                                          & (minimums[comp[0]] <= proj_data[comp])]

            # Make scatter plot of change for project.
            plot_scatter_by_nt_change.add_trace(go.Scatter(x=proj_data[comp], y=proj_data[change],
                                                           legendgroup=proj, showlegend=False,
                                                           name=proj, mode="markers"),
                                                row=i // 2 + 1, col=i % 2 + 1)
            scatter_range[0] += [proj_data[change].min(), proj_data[comp].min()]
            scatter_range[1] += [proj_data[change].max(), proj_data[comp].max()]

            # Make box plot of change for project.
            proj_data = proj_data.loc[~ pd.isna(proj_data[ratio])]
            proj_data = proj_data.loc[proj_data[ratio] != inf]
            plot_box_ratio_by_nt_change.add_trace(go.Box(y=proj_data[ratio], name=proj,
                                                         legendgroup=proj, showlegend=False),
                                                  row=i // 2 + 1, col=i % 2 + 1)
            box_range[0].append(proj_data[ratio].min())
            box_range[1].append(proj_data[ratio].max())

            # Fix colors for consistency across traces for both scatter and box.
            for plot in (plot_scatter_by_nt_change, plot_box_ratio_by_nt_change):
                plot.update_traces(
                    marker=dict(color=plot.layout["template"]["layout"]["colorway"][j]),
                    selector=dict(name=proj)
                    )

        # Add y=x line to scatter plot.
        plot_scatter_by_nt_change.add_trace(go.Scatter(x=[0, change_comparison_data[comp].max()],
                                                       y=[0, change_comparison_data[comp].max()],
                                                       mode="lines", marker_color="rgba(0, 0, 0, .25)",
                                                       name="y=x", legendgroup="y=x", showlegend=False),
                                            row=i // 2 + 1, col=i % 2 + 1)

    # Adjust axes and add titles.
    box_range = _make_axis_range(*box_range)
    scatter_range = _make_axis_range(*scatter_range)
    for i, (change, comp) in enumerate(zip(CHANGES, comps)):
        change, comp = make_change_label(change), make_change_label(comp)
        box_y_label = f"{change} / {comp}"
        if log_transform_ratio:
            box_y_label = f"log2({box_y_label})"
        plot_box_ratio_by_nt_change.layout[f"yaxis{i+1}"]["title"]["text"] = box_y_label
        if shared_yscale in ["box", "both", True]:
            plot_box_ratio_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = box_range
            # plot_box_ratio_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = [0, 1]
        if shared_yscale in ["scatter", "both", True]:
            plot_scatter_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = scatter_range
        if normalization_counts is not None:
            if use_scaled_labels:
                change, comp = (make_scaled_label(x, normalization_factor) for x in (change, comp))
        plot_scatter_by_nt_change.layout[f"xaxis{i+1}"]["title"]["text"] = comp
        plot_scatter_by_nt_change.layout[f"yaxis{i+1}"]["title"]["text"] = change
    # if shared_yscale in ["box", "both", True]:
    #     plot_box_ratio_by_nt_change.update_yaxes(nticks=20)

    # Show one trace legend in each plot.
    plot_box_ratio_by_nt_change.update_traces(showlegend=True, row=1, col=1)
    plot_scatter_by_nt_change.update_traces(showlegend=True, row=1, col=1)

    # Move boxplot legend.
    plot_box_ratio_by_nt_change.update_layout(legend=dict(
        orientation="h", yanchor="bottom", title_text="TCGA Project:", x=0.5, xanchor="center",
        ))

    # Add titles.
    plot_box_ratio_by_nt_change.update_layout(
        title=f"Mismatch / complement {'log2 ' if log_transform_ratio else ''}ratios per TCGA project"
        )
    plot_scatter_by_nt_change.update_layout(title="Mismatch vs. complement counts per TCGA project")

    # Remove box plot tick labels.
    # plot_box_ratio_by_nt_change.update_xaxes(showticklabels=False)

    return plot_box_ratio_by_nt_change, plot_scatter_by_nt_change


def plot_horizon(mismatch_table, filter_outliers=False, iqr_factor=1.5,
                 normalize_by_nt_content=False, normalization_factor=1,
                 use_scaled_labels=True, log_transform_ratio=False, shared_yscale=False,
                 minimum_raw_count=None):
    normalization_counts = VCROME.sum(axis=0) if normalize_by_nt_content else None
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    plot_scatter_by_nt_change = make_subplots(1, 6, start_cell="bottom-left")
    scatter_range = [0, 0]

    comps = [reverse_complement(change) for change in CHANGES]
    for i, (change, comp) in enumerate(zip(CHANGES, comps)):

        change_comparison_data = (mismatch_table.mismatch_tools.simplify_for_reference_asymmetry(
            change_numerator=change,
            change_denominator=comp,
            log_transform_ratio=log_transform_ratio,
            normalization_counts=normalization_counts,
            normalization_factor=normalization_factor
            ))

        for j, proj in enumerate(_proj_numbers.keys()):

            # Pre-subset and treat data.
            proj_data = change_comparison_data.loc[idx[proj, :],]
            if filter_outliers:
                proj_data = filter_by_fences(proj_data, [change, comp], iqr_factor)
            if minimum_raw_count is not None:
                proj_data = proj_data.loc[(minimums[change[0]] <= proj_data[change])
                                          & (minimums[comp[0]] <= proj_data[comp])]

            # Make scatter plot of change for project.
            plot_scatter_by_nt_change.add_trace(go.Scatter(x=proj_data[comp], y=proj_data[change],
                                                           legendgroup=proj, showlegend=False,
                                                           name=proj, mode="markers"),
                                                row=1, col=i+1)
            scatter_range = [min(scatter_range[0], proj_data[change].min(), proj_data[comp].min()),
                             max(scatter_range[1], proj_data[change].max(), proj_data[comp].max())]

            # Fix colors for consistency across traces for both scatter and box.
            plot_scatter_by_nt_change.update_traces(
                marker=dict(color=plot_scatter_by_nt_change.layout["template"]["layout"]["colorway"][j]),
                selector=dict(name=proj)
                )

        # Add y=x line to scatter plot.
        plot_scatter_by_nt_change.add_trace(go.Scatter(x=[0, change_comparison_data[comp].max()],
                                                       y=[0, change_comparison_data[comp].max()],
                                                       mode="lines", marker_color="rgba(0, 0, 0, .25)",
                                                       name="y=x", legendgroup="y=x", showlegend=False),
                                            row=1, col=i + 1)

    # Adjust axes and add titles.
    for i, (change, comp) in enumerate(zip(CHANGES, comps)):
        change, comp = make_change_label(change), make_change_label(comp)
        if shared_yscale in ["scatter", "both", True]:
            plot_scatter_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = scatter_range
        if normalization_counts is not None:
            if use_scaled_labels:
                change, comp = (make_scaled_label(x, normalization_factor) for x in (change, comp))
        plot_scatter_by_nt_change.layout[f"xaxis{i+1}"]["title"]["text"] = comp
        plot_scatter_by_nt_change.layout[f"yaxis{i+1}"]["title"]["text"] = change

    # Show one trace legend in each plot.
    plot_scatter_by_nt_change.update_traces(showlegend=True, row=1, col=1)

    # Add titles.
    plot_scatter_by_nt_change.update_layout(title="Mismatch vs. complement counts per 1k nts per TCGA project")

    # Remove box plot tick labels.
    # plot_box_ratio_by_nt_change.update_xaxes(showticklabels=False)

    return plot_scatter_by_nt_change


def plot_call_ref_asym_per_change_per_project(call_table,
                                              filter_outliers=False, iqr_factor=1.5,
                                              normalize_by_nt_content=False, normalization_factor=1,
                                              log_transform_ratio=False, shared_yscale=False,
                                              minimum_raw_count=None):
    normalization_counts = VCROME.sum(axis=0) if normalize_by_nt_content else None
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    box = make_subplots(3, 2, start_cell="bottom-left")
    scatter = make_subplots(3, 4, start_cell="bottom-left")
    box_range = [[], []]
    scatter_range = [[], []]

    comps = [reverse_complement(change) for change in CHANGES]
    for i, (change, comp) in enumerate(zip(CHANGES, comps)):
        ratio = f"{change}{comp}_ratio".lower()

        change_comparison_data = (
            call_table.call_tools.simplify_for_reference_asymmetry(change_numerator=change,
                                                                   change_denominator=comp,
                                                                   log_transform_ratio=log_transform_ratio,
                                                                   normalization_counts=normalization_counts,
                                                                   normalization_factor=normalization_factor)
            )

        for j, proj in enumerate(_proj_numbers.keys()):
            for k, filtered in enumerate(["PASS"], start=1):
            # for k, filtered in enumerate(["PASS", "FAIL"], start=1):
            # for k, filtered in enumerate(["FAIL"], start=1):

                # Pre-subset and filter data.
                proj_data = change_comparison_data.loc[idx[proj, :, :],][filtered]
                if filter_outliers:
                    proj_data = filter_by_fences(proj_data, [change, comp], iqr_factor)
                if minimum_raw_count is not None:
                    proj_data = proj_data.loc[(minimums[change[0]] <= proj_data[change])
                                              & (minimums[comp[0]] <= proj_data[comp])]

                # Make scatter plot of change for project.
                scatter.add_trace(go.Scatter(x=proj_data[comp], y=proj_data[change],
                                             legendgroup=proj, showlegend=False,
                                             name=proj, mode="markers"),
                                  row=i // 2 + 1, col=i % 2 * 2 + k)
                scatter_range[0] += [proj_data[change].min(), proj_data[comp].min()]
                scatter_range[1] += [proj_data[change].max(), proj_data[comp].max()]

                # Make box plot of change for project.
                proj_data = proj_data.loc[~ pd.isna(proj_data[ratio])]
                proj_data = proj_data.loc[proj_data[ratio].abs() != inf]
                box.add_trace(go.Box(y=proj_data[ratio], name=f"{proj} {filtered}",
                                     legendgroup=proj, showlegend=False),
                              row=i // 2 + 1, col=i % 2 + 1)
                box_range[0].append(proj_data[ratio].min())
                box_range[1].append(proj_data[ratio].max())
                # Fix colors for consistency across traces for both scatter and box.
                for plot in (scatter, box):
                    plot.update_traces(
                        marker=dict(color=plot.layout["template"]["layout"]["colorway"][j]),
                        selector=dict(legendgroup=proj)
                        )

        # Add y=x line to scatter plot.
        for k, filtered in enumerate(["PASS", "FAIL"], start=1):
            xy_line = go.Scatter(
                x=[0, change_comparison_data[(filtered, change)].max()],
                y=[0, change_comparison_data[(filtered, comp)].max()],
                mode="lines", marker_color="rgba(0, 0, 0, .25)",
                name="y = x", legendgroup="y = x", showlegend=False
                )
            scatter.add_trace(xy_line, row=i // 2 + 1, col=i % 2 * 2 + k)

    # Adjust axes and add titles.
    box_range = _make_axis_range(*box_range)
    scatter_range = _make_axis_range(*scatter_range)
    for i, (change, comp) in enumerate(zip(CHANGES, comps), start=1):
        change, comp = make_change_label(change), make_change_label(comp)
        box_y_label = f"{change} / {comp}"
        if log_transform_ratio:
            box_y_label = f"log2({box_y_label})"
        box.layout[f"yaxis{i}"]["title"]["text"] = box_y_label
        if shared_yscale in ["box", "both", True]:
            box.layout[f"yaxis{i}"]["range"] = box_range

    if shared_yscale in ["scatter", "both", True]:
        for k in range(1, 13):
            scatter.layout[f"yaxis{k}"]["range"] = scatter_range

    for i, (change, comp) in enumerate(zip(CHANGES, comps)):
        if normalization_counts is not None:
            change, comp = (make_scaled_label(x, normalization_factor) for x in (change, comp))
        for k in range(1, 3):
            scatter.layout[f"xaxis{i*2+k}"]["title"]["text"] = comp
            scatter.layout[f"yaxis{i*2+k}"]["title"]["text"] = change

    # Show one trace legend in each plot.
    box.update_traces(showlegend=True, row=1, col=1)
    scatter.update_traces(showlegend=True, row=1, col=1)

    # Move boxplot legend.
    box.update_layout(legend=dict(
        orientation="h", yanchor="bottom", title_text="TCGA Project:", x=0.5, xanchor="center",
        ))

    # Add titles.
    box.update_layout(
        title=f"Call / complement {'log2 ' if log_transform_ratio else ''}ratios per TCGA project"
        )
    scatter.update_layout(title="Call vs. complement counts per TCGA project")

    # Remove box plot tick labels.
    # plot_box_ratio_by_nt_change.update_xaxes(showticklabels=False)

    return box, scatter


def plot_oxog_calls_vs_mismatches_per_project(mismatch_table, call_table,
                                              minimum_mismatch_count=None):
    minimums = _normalize_minimum(minimum_mismatch_count)

    plots = dict()
    mismatches = mismatch_table.mismatch_tools.simplify_for_reference_asymmetry("GT", "CA").droplevel("file_id")

    calls = filter_raw_call_data_dups(call_table)
    calls = calls.call_tools.simplify_for_reference_asymmetry("GT", "CA")

    for filter_status in ["PASS", "FAIL"]:
        scatter = go.Figure()
        combo = calls[filter_status].join(mismatches,
                                          lsuffix="_calls",
                                          rsuffix="_mismatches")
        if minimum_mismatch_count is not None:
            combo = combo.loc[(minimums["G"] <= combo["GT_mismatches"])
                              & (minimums["C"] <= combo["CA_mismatches"])]
        for field in ["gtca_ratio_calls", "gtca_ratio_mismatches"]:
            combo = combo.loc[~ pd.isna(combo[field])]
            combo = combo.loc[combo[field].abs() != inf]

        for i, proj in enumerate(_proj_numbers.keys()):
            proj_combo = combo.loc[idx[proj, :],]
            scatter.add_trace(go.Scatter(x=proj_combo["gtca_ratio_mismatches"],
                                         y=proj_combo["gtca_ratio_calls"],
                                         mode="markers",
                                         name=proj))
            scatter.update_traces(
                    marker=dict(color=scatter.layout["template"]["layout"]["colorway"][i]),
                    selector=dict(name=proj)
                    )

        scatter.update_layout(title=f"G➔T / C➔A ratio of mismatches vs. "
                                    f"{filter_status.lower()} calls per TCGA project",
                              xaxis_title="G➔T / C➔A mismatch ratio",
                              yaxis_title="G➔T / C➔A call ratio")

        plots[filter_status] = scatter
    return plots


# %% Make mismatch and call reference asymmetry plots.

mismatch_plots = plot_mismatch_ref_asym_per_change_per_project(
    mismatch_table=mismatch_data,
    log_transform_ratio=False,
    normalize_by_nt_content=True,
    normalization_factor=1000,
    shared_yscale="box")


call_plots = plot_call_ref_asym_per_change_per_project(
    call_table=call_data,
    # normalize_by_nt_content=True,
    # normalization_factor=1000000,
    log_transform_ratio=False,
    shared_yscale="box")


# %% Binomial tests comparing incidence of each nt change pair.

def calculate_binomtest_for_mismatches(mismatch_table, filter_outliers=False, iqr_factor=1.5,
                                       normalize_by_nt_content=False, normalization_factor=1,
                                       minimum_raw_count=None, changes=None):
    normalization_counts = VCROME.sum(axis=0) if normalize_by_nt_content else None
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    if changes is None:
        changes = CHANGES
    binom_results = dict()

    # Binomial test per nucleotide change pair.
    for i, change in enumerate(changes):
        binom_change = dict()

        comp = reverse_complement(change)
        change_comparison_data = (mismatch_table.mismatch_tools.simplify_for_reference_asymmetry(
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
        changes = CHANGES
    binom_results = dict()

    # Binomial test per nucleotide change pair.
    for i, change in enumerate(changes):
        for filter_status in ["PASS", "FAIL"]:
            binom_change = dict()

            comp = reverse_complement(change)
            change_comparison_data = (
                call_table.call_tools.simplify_for_reference_asymmetry(change_numerator=change,
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


def filter_raw_call_data_dups(call_table):
    calls = (call_table.groupby(["project_id", "case_id", "file_id"]).agg(sum))
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


def calculate_call_mismatch_correlation(call_table, mismatch_table,
                                        minimum_mismatch_count=None,
                                        corr_test="spearman"):
    corr_test = {"spearman": stats.spearmanr,
                 "pearson": stats.pearsonr}[corr_test]
    minimums = _normalize_minimum(minimum_mismatch_count)

    corr_results = defaultdict(dict)

    call_table = filter_raw_call_data_dups(call_table)

    # Correlation test per nucleotide change pair.
    for i, change in enumerate(CHANGES):

        comp = reverse_complement(change)
        ratio = f"{change}{comp}_ratio".lower()
        mismatch_subset = (mismatch_table.mismatch_tools.simplify_for_reference_asymmetry(change, comp)
                           .droplevel("file_id"))

        call_subset = call_table.call_tools.simplify_for_reference_asymmetry(change, comp)
        # call_subset = filter_call_dups(call_subset).droplevel("file_id")

        for filter_status in ["PASS", "FAIL"]:
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
    comps = [reverse_complement(change) for change in CHANGES]
    steps = len(CHANGES)*len(afs)*len(projects)*2
    i = 1
    for change, comp in zip(CHANGES, comps):
        mismatches = mismatch_table
        for pass_fail in ["PASS", "FAIL"]:
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
    comps = [reverse_complement(change) for change in CHANGES]
    steps = len(CHANGES)*len(depths)*len(projects)*2
    i = 1
    for change, comp in zip(CHANGES, comps):
        mismatches = mismatch_table
        for pass_fail in ["PASS", "FAIL"]:
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

    #             calls = call_table.loc[call_table.frequency]
    # for l1_a, (change, comp) in zip(range(0, 12, 2),
    #                               zip(CHANGES, comps)):
    #     mismatches = (mismatch_table
    #                   .mismatch_tools
    #                   .simplify_for_reference_asymmetry(change, comp))
    #     for l1_b, pass_fail in enumerate(["PASS", "FAIL"]):
    #         l1 = l1_a + l1_b
    #         for l2, proj in enumerate(_proj_numbers.keys()):
    #             for l3, af in enumerate(np.linspace(0, 1, 201)):
    #                 calls = call_table.loc[call_table.frequency >= af]
    #                 calls = (calls
    #                          .call_tools
    #                          .simplify_for_reference_asymmetry(change, comp))[pass_fail]
    #                 combo = calls.join(mismatches, lsuffix="_calls", rsuffix="_mismatches")
    #                 corr_test = stats.spearmanr(combo[f"{change}{comp}_ratio_mismatches"],
    #                                             combo[f"{change}{comp}_ratio_calls"])


# %% Binomial test of reference call asymmetry in aggregate, counting each site once across
# all samples.

def calculate_binomtest_for_unique_calls(call_table, changes=None):
    if changes is None:
        changes = CHANGES
    filters = ["PASS", "FAIL"]
    binom_results = dict()

    for proj in list(_proj_numbers) + [slice(None)]:
        proj_dict = dict()
        uniques = ((call_table.loc[proj].groupby(["contig", "pos", "ref", "alt"]).agg(sum) > 0)
                   .astype(int)
                   .groupby(["ref", "alt"])
                   .agg(sum))
        for pass_fail in filters:
            filtered = uniques[pass_fail]
            for change in changes:
                comp = reverse_complement(change)
                change_count = filtered.loc[idx[tuple(change)]]
                comp_count = filtered.loc[idx[tuple(comp)]]
                test = stats.binomtest(change_count, change_count + comp_count).pvalue
                proj_dict[(pass_fail, change)] = test
        binom_results[proj if isinstance(proj, str) else "total"] = proj_dict
    binom_results = pd.DataFrame.from_dict(binom_results, orient="index")
    return binom_results


def calculate_binomtest_for_mismatches(mismatch_table, changes=None, p=None):
    if changes is None:
        changes = CHANGES
    if p is None:
        p = {change: 0.5 for change in changes}
    binom_results = dict()

    for change in changes:
        comp = reverse_complement(change)
        subset = (mismatch_table.mismatch_tools
                  .simplify_for_reference_asymmetry(change, comp,
                                                    add_ratio_column=False,
                                                    add_fraction_column=False))
        subset = subset.groupby("project_id").agg(sum)
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
                .agg(sum) > 0).astype(int))["PASS"]
all_unique_counts = all_uniques.loc[all_uniques > 0].groupby(["project_id", "ref", "alt"]).agg(sum)


def plot_aggregate_oxog_frequencies(call_table):
    all_uniques = ((call_table
                    .groupby(["project_id", "contig", "pos", "ref", "alt"])
                    .agg(sum) > 0).astype(int))["PASS"]
    all_unique_counts = all_uniques.loc[all_uniques > 0].groupby(["project_id", "ref", "alt"]).agg(sum)
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


# %% Check variant contexts of TGCT oxo-G variants.

# unique_tgct = ((call_data.loc["TCGA-TGCT"]
#                 .groupby(["contig", "pos", "ref", "alt"])
#                 .agg(sum) > 0).astype(int))["PASS"]
#
# contexts = dict()
# for change in ["GT", "CA"]:
#     strings = [f"{x.contig}\t{x.pos - 2}\t{x.pos+1}"
#                for x in (unique_tgct
#                          .loc[unique_tgct > 0]
#                          .loc[idx[:, :, change[0], change[1]]]
#                          .reset_index()
#                          .itertuples())]
#     bedtool = BedTool("\n".join(strings), True)
#     bedtool.sequence(fi="/home/tyler/Documents/Resource_Data/reference_genomes/GRCh38.d1.vd1.fa")
#     with open(bedtool.seqfn) as infile:
#         context_results = infile.readlines()
#     context_results = [x.rstrip() for x in context_results if not x.startswith(">")]
#     contexts[change] = Counter(context_results)


# %% Check variant contexts of all oxo-G variants.

def count_oxog_contexts(call_table, ref_genome):
    all_uniques = ((call_table
                    .groupby(["project_id", "contig", "pos", "ref", "alt"])
                    .agg(sum) > 0).astype(int))["PASS"]

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


oxog_contexts = count_oxog_contexts(call_data, "/home/tyler/Documents/Resource_Data/reference_genomes/GRCh38.d1.vd1.fa")


# %% Binomial simulation of TGCT.

def plot_tgct_oxog_simulation(call_table, iterations):
    prob = VCROME.G.sum() / (VCROME.G.sum() + VCROME.C.sum())
    tgct_oxog = (call_table.call_tools.simplify_for_reference_asymmetry("GT", "CA")
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

# bam_data = gdc_api.retrieve_all_data()
# bam_data = gdc_api.filter_for_paired_sequences(bam_data, true_pairs=True)
# bam_data = gdc_api.filter_by_capture_kit_number(bam_data, "06 465 668 001", include_wgs=False)
# bam_data = bam_data.reset_index().set_index([("cases", "project.project_id"), "case_id", "file_id"])
# bam_data.index.names = ["project_id", "case_id", "file_id"]
#
# oxo_by_date = mismatch_data.mismatch_tools.simplify_for_reference_asymmetry("GT", "CA")
# oxo_by_date = oxo_by_date.join(bam_data.analysis.sequencing_date)
# oxo_by_date.sequencing_date = [datetime.strptime(x, "%Y-%m-%dT%H") if isinstance(x, str) else x
#                                for x in oxo_by_date.sequencing_date]
#
# oxo_date_plot = px.scatter(oxo_by_date.reset_index(), x="sequencing_date", y="gtca_ratio", color="project_id")
