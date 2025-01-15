import pandas as pd
import plotly.graph_objects as go
from numpy import inf, log2
from pandas import IndexSlice as idx
from plotly.subplots import make_subplots

from utilities import (
    _CHANGES, _COMPS,
    _normalize_minimum,
    make_complement,
    filter_by_fences,
    filter_raw_call_data_dups
    )


def make_change_labels(change_numerators):
    change_labels = [[change, make_complement(change)] for change in change_numerators]
    change_labels = [["➔".join(list(change)) for change in changes]
                     for changes in change_labels]
    change_labels = ["<br>".join([f"  {changes[0]} /", changes[1]]) for changes in change_labels]
    return change_labels


def make_change_label(change_string):
    change_string = "➔".join(list(change_string))
    return change_string


SCALE_LABELS = {1: "", 1000: "1K ", 1000000: "1M"}


def make_scaled_label(change_label, scale_factor):
    # label = f"{change_label} per {SCALE_LABELS[scale_factor]}{change_label[0]} nts"
    label = f"{change_label} per {scale_factor} {change_label[0]}"
    return label


def _make_axis_range(mins, maxes):
    axis_min, axis_max = min(mins), max(maxes)
    inter_range = axis_max - axis_min
    axis_range = [axis_min - inter_range/5, axis_max+inter_range/5]
    return axis_range


# def plot_mismatch_asym_per_change_per_project(mismatch_table, bias_type, filter_outliers=False, iqr_factor=1.5,
#                                               normalize_by_nt_content=False, normalization_factor=1,
#                                               normalization_data=None,
#                                               use_scaled_labels=True, log_transform_ratio=False,
#                                               shared_yscale=False,
#                                               minimum_raw_count=None):
#     normalization_counts = normalization_data.sum(axis=0) if normalize_by_nt_content else None
#     minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
#     _proj_numbers = {y: x for x, y in enumerate(mismatch_table.index.levels[0], start=-4)}
#     plot_box_ratio_by_nt_change = make_subplots(3, 2, start_cell="bottom-left")
#     plot_scatter_by_nt_change = make_subplots(3, 2, start_cell="bottom-left")
#     box_range = [[], []]
#     scatter_range = [[], []]
#
#     comps = [reverse_complement(change) for change in _CHANGES]
#     for i, (change, comp) in enumerate(zip(_CHANGES, comps)):
#         ratio = f"{change}{comp}_ratio".lower()
#
#         change_comparison_data = (mismatch_table.mismatch_tools.calculate_asymmetry(
#             bias_type=bias_type,
#             change_numerator=change,
#             change_denominator=comp,
#             log_transform_ratio=log_transform_ratio,
#             normalization_counts=normalization_counts,
#             normalization_factor=normalization_factor
#             ))
#
#         for j, proj in enumerate(_proj_numbers.keys()):
#
#             # Pre-subset and treat data.
#             proj_data = change_comparison_data.loc[idx[proj, :],]
#             if filter_outliers:
#                 proj_data = filter_by_fences(proj_data, [change, comp], iqr_factor)
#             if minimum_raw_count is not None:
#                 proj_data = proj_data.loc[(minimums[change[0]] <= proj_data[change])
#                                           & (minimums[comp[0]] <= proj_data[comp])]
#
#             # Make scatter plot of change for project.
#             plot_scatter_by_nt_change.add_trace(go.Scatter(x=proj_data[comp], y=proj_data[change],
#                                                            legendgroup=proj, showlegend=False,
#                                                            name=proj, mode="markers"),
#                                                 row=i // 2 + 1, col=i % 2 + 1)
#             scatter_range[0] += [proj_data[change].min(), proj_data[comp].min()]
#             scatter_range[1] += [proj_data[change].max(), proj_data[comp].max()]
#
#             # Make box plot of change for project.
#             proj_data = proj_data.loc[~ pd.isna(proj_data[ratio])]
#             proj_data = proj_data.loc[proj_data[ratio] != inf]
#             plot_box_ratio_by_nt_change.add_trace(go.Box(y=proj_data[ratio], name=proj,
#                                                          legendgroup=proj, showlegend=False),
#                                                   row=i // 2 + 1, col=i % 2 + 1)
#             box_range[0].append(proj_data[ratio].min())
#             box_range[1].append(proj_data[ratio].max())
#
#             # Fix colors for consistency across traces for both scatter and box.
#             for plot in (plot_scatter_by_nt_change, plot_box_ratio_by_nt_change):
#                 plot.update_traces(
#                     marker=dict(color=plot.layout["template"]["layout"]["colorway"][j]),
#                     selector=dict(name=proj)
#                     )
#
#         # Add y=x line to scatter plot.
#         plot_scatter_by_nt_change.add_trace(go.Scatter(x=[0, change_comparison_data[comp].max()],
#                                                        y=[0, change_comparison_data[comp].max()],
#                                                        mode="lines", marker_color="rgba(0, 0, 0, .25)",
#                                                        name="y=x", legendgroup="y=x", showlegend=False),
#                                             row=i // 2 + 1, col=i % 2 + 1)
#
#     # Adjust axes and add titles.
#     box_range = _make_axis_range(*box_range)
#     scatter_range = _make_axis_range(*scatter_range)
#     for i, (change, comp) in enumerate(zip(_CHANGES, comps)):
#         change, comp = make_change_label(change), make_change_label(comp)
#         box_y_label = f"{change} / {comp}"
#         if log_transform_ratio:
#             box_y_label = f"log<sub>2</sub>({box_y_label})"
#         plot_box_ratio_by_nt_change.layout[f"yaxis{i+1}"]["title"]["text"] = box_y_label
#         if shared_yscale in ["box", "both", True]:
#             plot_box_ratio_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = box_range
#             # plot_box_ratio_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = [0, 1]
#         if shared_yscale in ["scatter", "both", True]:
#             plot_scatter_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = scatter_range
#         if normalization_counts is not None:
#             if use_scaled_labels:
#                 change, comp = (make_scaled_label(x, normalization_factor) for x in (change, comp))
#         plot_scatter_by_nt_change.layout[f"xaxis{i+1}"]["title"]["text"] = comp
#         plot_scatter_by_nt_change.layout[f"yaxis{i+1}"]["title"]["text"] = change
#     # if shared_yscale in ["box", "both", True]:
#     #     plot_box_ratio_by_nt_change.update_yaxes(nticks=20)
#
#     # Show one trace legend in each plot.
#     plot_box_ratio_by_nt_change.update_traces(showlegend=True, row=1, col=1)
#     plot_scatter_by_nt_change.update_traces(showlegend=True, row=1, col=1)
#
#     # Move boxplot legend.
#     plot_box_ratio_by_nt_change.update_layout(legend=dict(
#         orientation="h", yanchor="bottom", title_text="TCGA Project:", x=0.5, xanchor="center",
#         ))
#
#     # Add titles.
#     plot_box_ratio_by_nt_change.update_layout(
#         title=f"Asymmetry by {bias_type} strand: mismatch vs. complement "
#               f"{'log<sub>2</sub>' if log_transform_ratio else ''}ratio per TCGA project"
#         )
#     plot_scatter_by_nt_change.update_layout(title="Reference strand mismatch vs. complement counts per TCGA project")
#
#     # Remove box plot tick labels.
#     # plot_box_ratio_by_nt_change.update_xaxes(showticklabels=False)
#
#     return plot_box_ratio_by_nt_change, plot_scatter_by_nt_change


def plot_mismatch_asymmetry(mismatch_table, bias_type,
                            filter_outliers=False, iqr_factor=1.5,
                            normalize_by_nt_content=False, normalization_factor=1,
                            normalization_data=None, use_scaled_labels=True,
                            log_transform_ratio=False,
                            shared_yscale=False,
                            minimum_raw_count=None, seq_type="WES"):
    # normalization_counts = normalization_data.sum(axis=0) if normalize_by_nt_content else None
    normalization_counts = None
    if normalize_by_nt_content:
        normalization_counts = normalization_data.capkit.calculate_counts(bias_type)
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    _proj_numbers = {y: x for x, y in enumerate(mismatch_table.index.levels[0], start=-4)}
    scatters = make_subplots(3, 2, start_cell="bottom-left")
    boxes = make_subplots(3, 2, start_cell="bottom-left")
    ridgelines = make_subplots(3, 2, start_cell="bottom-left", shared_xaxes=False)
    scatter_range, box_range, = [[], []], [[], []]

    for i, (change, comp) in enumerate(zip(_CHANGES, _COMPS)):
        ratio = f"{change}{comp}_ratio".lower()

        plot_data = (mismatch_table.mismatch_tools.calculate_asymmetry(
            bias_type=bias_type,
            change_numerator=change,
            change_denominator=comp,
            log_transform_ratio=log_transform_ratio,
            normalization_counts=normalization_counts,
            normalization_factor=normalization_factor
            ))

        for j, proj in enumerate(_proj_numbers.keys()):

            # Pre-subset and treat data.
            proj_data = plot_data.loc[idx[proj, :],]
            if filter_outliers:
                proj_data = filter_by_fences(proj_data, [change, comp], iqr_factor)
            if minimum_raw_count is not None:
                proj_data = proj_data.loc[(minimums[change[0]] <= proj_data[change])
                                          & (minimums[comp[0]] <= proj_data[comp])]

            # Make scatter.
            scatters.add_trace(go.Scatter(x=proj_data[comp], y=proj_data[change],
                                          legendgroup=proj, showlegend=False,
                                          name=proj[5:], mode="markers"),
                               row=i // 2 + 1, col=i % 2 + 1)
            scatter_range[0] += [proj_data[change].min(), proj_data[comp].min()]
            scatter_range[1] += [proj_data[change].max(), proj_data[comp].max()]

            # Make box.
            proj_data = proj_data.loc[~ pd.isna(proj_data[ratio])]
            proj_data = proj_data.loc[proj_data[ratio] != inf]
            boxes.add_trace(go.Box(y=proj_data[ratio], name=proj[5:],
                                   legendgroup=proj, showlegend=False),
                            row=i // 2 + 1, col=i % 2 + 1)
            box_range[0].append(proj_data[ratio].min())
            box_range[1].append(proj_data[ratio].max())

            # Make ridgeline.
            ridgelines.add_trace(go.Violin(x=proj_data[ratio], name=proj[5:], legendgroup=proj,
                                           showlegend=False, zorder=j),
                                 row=i // 2 + 1, col=i % 2 + 1)

            # Fix colors for consistency across traces.
            for plot in (scatters, boxes, ridgelines):
                plot.update_traces(
                    marker=dict(color=plot.layout["template"]["layout"]["colorway"][j]),
                    selector=dict(name=proj[5:])
                    )

        # Add y=x line to scatter plot.
        scatters.add_trace(go.Scatter(x=[0, plot_data[comp].max()],
                                      y=[0, plot_data[comp].max()],
                                      mode="lines", marker_color="rgba(0, 0, 0, .25)",
                                      line={"dash": "dash"},
                                      name="y=x", legendgroup="y=x", showlegend=False,
                                      zorder=len(_proj_numbers)+1),
                           row=i // 2 + 1, col=i % 2 + 1)

    # Adjust axes and axis titles.
    box_range = _make_axis_range(*box_range)
    scatter_range = _make_axis_range(*scatter_range)
    for i, (change, comp) in enumerate(zip(_CHANGES, _COMPS)):
        change, comp = make_change_label(change), make_change_label(comp)
        box_y_label = f"{change} / {comp}"
        if log_transform_ratio:
            # box_y_label = f"log<sub>2</sub>({box_y_label})"
            box_y_label = r"$\log_2\frac{" + change + r"}{" + comp + r"}$"
        boxes.layout[f"yaxis{i+1}"]["title"]["text"] = box_y_label
        ridgelines.layout[f"xaxis{i+1}"]["title"]["text"] = box_y_label
        if shared_yscale in ["box", "both", True]:
            boxes.layout[f"yaxis{i + 1}"]["range"] = box_range
            ridgelines.layout[f"xaxis{i + 1}"]["range"] = box_range
            # plot_box_ratio_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = [0, 1]
        if shared_yscale in ["scatter", "both", True]:
            scatters.layout[f"yaxis{i + 1}"]["range"] = scatter_range
        if normalization_counts is not None and use_scaled_labels:
            change, comp = (make_scaled_label(x, normalization_factor) for x in (change, comp))
        scatters.layout[f"xaxis{i+1}"]["title"]["text"] = comp
        scatters.layout[f"yaxis{i+1}"]["title"]["text"] = change
    for i in range(1,7):
        ridgelines.update_layout(**{f"xaxis{i}_tickmode": "array",
                                    f"xaxis{i}_tickvals":list(range(0, 13))})
    # if shared_yscale in ["box", "both", True]:
    #     plot_box_ratio_by_nt_change.update_yaxes(nticks=20)

    # Show one trace legend in each plot.
    scatters.update_traces(showlegend=True, row=1, col=1)
    boxes.update_traces(showlegend=True, row=1, col=1)
    ridgelines.update_traces(showlegend=True, row=1, col=1)

    # Change plot layouts.
    boxes.update_layout(legend=dict(
        orientation="h", yanchor="bottom", title_text="TCGA Project:", x=0.5, xanchor="center",
        ))
    boxes.update_xaxes(showticklabels=False)
    ridgelines.update_traces(orientation="h", side="positive", width=3, points=False)
    ridgelines.update_yaxes(showticklabels=False)
    ridgelines.update_legends(traceorder="reversed")
    ridgelines.data = ridgelines.data[::-1]

    # Add titles.
    boxes.update_layout(
        title=f"TCGA {seq_type.upper()}<br>"
              f"Asymmetry by {bias_type} strand<br>"
              f"Mismatch vs. complement {'log<sub>2</sub>' if log_transform_ratio else ''}ratio"
        )
    ridgelines.update_layout(title=boxes.layout.title)
    scatters.update_layout(title=f"TCGA {seq_type.upper()}<br>"
                                 f"Asymmetry by {bias_type} strand<br>"
                                 f"Mismatch vs. complement counts",
                           legend_title="TCGA Project")

    return boxes, ridgelines, scatters


# def plot_horizon(mismatch_table, filter_outliers=False, iqr_factor=1.5,
#                  normalize_by_nt_content=False, normalization_factor=1,
#                  normalization_data=None,
#                  use_scaled_labels=True, log_transform_ratio=False, shared_yscale=False,
#                  minimum_raw_count=None):
#     _proj_numbers = {y: x for x, y in enumerate(mismatch_table.index.levels[0], start=-4)}
#     normalization_counts = normalization_data.sum(axis=0) if normalize_by_nt_content else None
#     minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
#     plot_scatter_by_nt_change = make_subplots(1, 6, start_cell="bottom-left")
#     scatter_range = [0, 0]
#
#     comps = [reverse_complement(change) for change in _CHANGES]
#     for i, (change, comp) in enumerate(zip(_CHANGES, comps)):
#
#         change_comparison_data = (mismatch_table.mismatch_tools.calculate_reference_asymmetry(
#             change_numerator=change,
#             change_denominator=comp,
#             log_transform_ratio=log_transform_ratio,
#             normalization_counts=normalization_counts,
#             normalization_factor=normalization_factor
#             ))
#
#         for j, proj in enumerate(_proj_numbers.keys()):
#
#             # Pre-subset and treat data.
#             proj_data = change_comparison_data.loc[idx[proj, :],]
#             if filter_outliers:
#                 proj_data = filter_by_fences(proj_data, [change, comp], iqr_factor)
#             if minimum_raw_count is not None:
#                 proj_data = proj_data.loc[(minimums[change[0]] <= proj_data[change])
#                                           & (minimums[comp[0]] <= proj_data[comp])]
#
#             # Make scatter plot of change for project.
#             plot_scatter_by_nt_change.add_trace(go.Scatter(x=proj_data[comp], y=proj_data[change],
#                                                            legendgroup=proj, showlegend=False,
#                                                            name=proj, mode="markers"),
#                                                 row=1, col=i+1)
#             scatter_range = [min(scatter_range[0], proj_data[change].min(), proj_data[comp].min()),
#                              max(scatter_range[1], proj_data[change].max(), proj_data[comp].max())]
#
#             # Fix colors for consistency across traces for both scatter and box.
#             plot_scatter_by_nt_change.update_traces(
#                 marker=dict(color=plot_scatter_by_nt_change.layout["template"]["layout"]["colorway"][j]),
#                 selector=dict(name=proj)
#                 )
#
#         # Add y=x line to scatter plot.
#         plot_scatter_by_nt_change.add_trace(go.Scatter(x=[0, change_comparison_data[comp].max()],
#                                                        y=[0, change_comparison_data[comp].max()],
#                                                        mode="lines", marker_color="rgba(0, 0, 0, .25)",
#                                                        name="y=x", legendgroup="y=x", showlegend=False),
#                                             row=1, col=i + 1)
#
#     # Adjust axes and add titles.
#     for i, (change, comp) in enumerate(zip(_CHANGES, comps)):
#         change, comp = make_change_label(change), make_change_label(comp)
#         if shared_yscale in ["scatter", "both", True]:
#             plot_scatter_by_nt_change.layout[f"yaxis{i + 1}"]["range"] = scatter_range
#         if normalization_counts is not None:
#             if use_scaled_labels:
#                 change, comp = (make_scaled_label(x, normalization_factor) for x in (change, comp))
#         plot_scatter_by_nt_change.layout[f"xaxis{i+1}"]["title"]["text"] = comp
#         plot_scatter_by_nt_change.layout[f"yaxis{i+1}"]["title"]["text"] = change
#
#     # Show one trace legend in each plot.
#     plot_scatter_by_nt_change.update_traces(showlegend=True, row=1, col=1)
#
#     # Add titles.
#     plot_scatter_by_nt_change.update_layout(title="Mismatch vs. complement counts per 1k nts per TCGA project")
#
#     # Remove box plot tick labels.
#     # plot_box_ratio_by_nt_change.update_xaxes(showticklabels=False)
#
#     return plot_scatter_by_nt_change


def plot_single_call_asymmetry_box(call_table, bias_type,
                                   numerator, denominator,
                                   filter_outliers=False, iqr_factor=1.5,
                                   normalize_by_nt_content=False, normalization_factor=1,
                                   normalization_data=None,
                                   log_transform_ratio=False, shared_yscale=False,
                                   minimum_raw_count=None,
                                   filters=None, seq_type="WES"):
    if filters is None:
        filters = ["TOTAL", "PASS"]
    _proj_numbers = {y: x for x, y in enumerate(call_table.index.levels[0], start=-4)}
    normalization_counts = None
    if normalize_by_nt_content:
        normalization_counts = normalization_data.capkit.calculate_counts(bias_type)
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    boxes = go.Figure()
    boxes.add_hline(y=1, line=dict(color="darkred", width=1.5, dash="dash"),
                    name="ratio = 1", showlegend=True)
    # box_range = [[], []]

    ratio = f"{numerator}{denominator}_ratio".lower()
    change_comparison_data = call_table.call_tools.calculate_asymmetry(
        bias_type=bias_type,
        change_numerator=numerator,
        change_denominator=denominator,
        add_fraction_column=False,
        log_transform_ratio=log_transform_ratio,
        normalization_counts=normalization_counts,
        normalization_factor=normalization_factor
        )
    box_max = 0
    for j, proj in enumerate(_proj_numbers.keys()):
        for k, filtered in enumerate(filters):

            # Pre-subset and filter data.
            proj_data = change_comparison_data.loc[idx[proj, :, :],][filtered]
            if filter_outliers:
                proj_data = filter_by_fences(proj_data, [numerator, denominator], iqr_factor)
            if minimum_raw_count is not None:
                proj_data = proj_data.loc[(minimums[numerator[0]] <= proj_data[numerator])
                                          & (minimums[denominator[0]] <= proj_data[denominator])]
            # Make box plot of change for project.
            name = proj[5:] if not k else proj[5:] + f" {filtered.title()}"
            proj_data = proj_data.loc[~ pd.isna(proj_data[ratio])]
            proj_data = proj_data.loc[proj_data[ratio].abs() != inf]
            box_max = max(box_max, proj_data[ratio].max())
            boxes.add_trace(go.Box(y=proj_data[ratio], name=name,
                                   legendgroup=filtered, showlegend=not k))
            # box_range[0].append(proj_data[ratio].min())
            # box_range[1].append(proj_data[ratio].max())

            boxes.for_each_trace(lambda trace: trace.update(
                marker_color=boxes.layout["template"]["layout"]["colorway"][j]
                ) if proj[5:] in trace.name else ())

    # Adjust axes and add titles.
    change, comp = make_change_label(numerator), make_change_label(denominator)
    box_y_label = f"{change} / {comp}"
    if log_transform_ratio:
        box_y_label = f"log<sub>2</sub>({box_y_label})"

    # Show one trace legend in each plot.
    # boxes.update_traces(showlegend=True, row=1, col=1)

    # Move boxplot legend.
    # boxes.update_layout(legend=dict(
    #     orientation="h", yanchor="bottom", title_text="TCGA Project:", x=0.5, xanchor="center",
    #     ))
    boxes.update_layout(legend_title_text="TCGA Project:")
    # Add titles.
    boxes.update_layout(
        title=f"TCGA {seq_type.upper()}<br>"
              f"Asymmetry by {bias_type} strand<br>"
              f"Call vs. complement {'log<sub>2</sub>' if log_transform_ratio else ''}ratio",
        yaxis_title=box_y_label,
        )

    boxes.update_layout(yaxis_tickmode="array",
                        yaxis_tickvals=list(range(0, round(box_max)+1)))

    # Remove box plot tick labels.
    boxes.update_xaxes(showticklabels=False)

    return boxes

def plot_call_asymmetry(call_table, bias_type,
                        filter_outliers=False, iqr_factor=1.5,
                        normalize_by_nt_content=False, normalization_factor=1,
                        normalization_data=None,
                        log_transform_ratio=False, shared_yscale=False,
                        minimum_raw_count=None,
                        filters=None, seq_type="WES"):
    if filters is None:
        filters = ["TOTAL", "PASS"]
    _proj_numbers = {y: x for x, y in enumerate(call_table.index.levels[0], start=-4)}
    normalization_counts = None
    if normalize_by_nt_content:
        normalization_counts = normalization_data.capkit.calculate_counts(bias_type)
    minimums = _normalize_minimum(minimum_raw_count, normalize_by_nt_content, normalization_factor)
    boxes = make_subplots(3, 2, start_cell="bottom-left")
    ridgelines = make_subplots(3, 2, start_cell="bottom-left", shared_xaxes=False)
    scatters = make_subplots(3, 4, start_cell="bottom-left")
    box_range = [[], []]
    scatter_range = [[], []]

    comps = [make_complement(change) for change in _CHANGES]
    for i, (change, comp) in enumerate(zip(_CHANGES, comps)):
        ratio = f"{change}{comp}_ratio".lower()

        change_comparison_data = (
            call_table.call_tools.calculate_asymmetry(bias_type=bias_type,
                                                      change_numerator=change,
                                                      change_denominator=comp,
                                                      add_fraction_column=False,
                                                      log_transform_ratio=log_transform_ratio,
                                                      normalization_counts=normalization_counts,
                                                      normalization_factor=normalization_factor)
            )

        for j, proj in enumerate(_proj_numbers.keys()):
            for k, filtered in enumerate(filters, start=1):
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
                scatters.add_trace(go.Scatter(x=proj_data[comp], y=proj_data[change],
                                              legendgroup=proj, showlegend=False,
                                              name=proj, mode="markers"),
                                   row=i // 2 + 1, col=i % 2 * 2 + k)
                scatter_range[0] += [proj_data[change].min(), proj_data[comp].min()]
                scatter_range[1] += [proj_data[change].max(), proj_data[comp].max()]

                # Make box plot of change for project.
                proj_data = proj_data.loc[~ pd.isna(proj_data[ratio])]
                proj_data = proj_data.loc[proj_data[ratio].abs() != inf]
                boxes.add_trace(go.Box(y=proj_data[ratio], name=f"{proj[5:]}, {filtered.title()} calls",
                                       legendgroup=filtered, showlegend=False),
                                row=i // 2 + 1, col=i % 2 + 1)
                box_range[0].append(proj_data[ratio].min())
                box_range[1].append(proj_data[ratio].max())

                # Make ridgeline plot.
                ridgelines.add_trace(go.Violin(x=proj_data[ratio], name=f"{proj[5:]}, {filtered.title()} calls",
                                               legendgroup=filtered, showlegend=False,
                                               zorder=j),
                                     row=i // 2 + 1, col=i % 2 + 1)

                # Fix colors for consistency across traces for both scatter and box.
                scatters.update_traces(
                    marker=dict(color=scatters.layout["template"]["layout"]["colorway"][j]),
                    selector=dict(legendgroup=proj)
                    )

                for plot in (boxes, ridgelines):
                    plot.for_each_trace(lambda trace: trace.update(
                        marker_color=plot.layout["template"]["layout"]["colorway"][j]
                        ) if proj in trace.name else ())

        # Add y=x line to scatter plot.
        for k, filtered in enumerate(["TOTAL", "PASS"], start=1):
            xy_line = go.Scatter(
                x=[0, change_comparison_data[(filtered, change)].max()],
                y=[0, change_comparison_data[(filtered, comp)].max()],
                mode="lines", marker_color="rgba(0, 0, 0, .25)",
                name="y = x", legendgroup="y = x", showlegend=False
                )
            scatters.add_trace(xy_line, row=i // 2 + 1, col=i % 2 * 2 + k)

    # Adjust axes and add titles.
    box_range = _make_axis_range(*box_range)
    scatter_range = _make_axis_range(*scatter_range)
    for i, (change, comp) in enumerate(zip(_CHANGES, comps), start=1):
        change, comp = make_change_label(change), make_change_label(comp)
        box_y_label = f"{change} / {comp}"
        if log_transform_ratio:
            box_y_label = f"log<sub>2</sub>({box_y_label})"
        boxes.layout[f"yaxis{i}"]["title"]["text"] = box_y_label
        ridgelines.layout[f"xaxis{i}"]["title"]["text"] = box_y_label
        if shared_yscale in ["box", "both", True]:
            boxes.layout[f"yaxis{i}"]["range"] = box_range
            ridgelines.layout[f"xaxis{i}"]["range"] = box_range
    if shared_yscale in ["scatter", "both", True]:
        for k in range(1, 13):
            scatters.layout[f"yaxis{k}"]["range"] = scatter_range

    for i, (change, comp) in enumerate(zip(_CHANGES, comps)):
        if normalization_counts is not None:
            change, comp = (make_scaled_label(x, normalization_factor) for x in (change, comp))
        for k in range(1, 3):
            scatters.layout[f"xaxis{i*2+k}"]["title"]["text"] = comp
            scatters.layout[f"yaxis{i*2+k}"]["title"]["text"] = change

    # Show one trace legend in each plot.
    boxes.update_traces(showlegend=True, row=1, col=1)
    scatters.update_traces(showlegend=True, row=1, col=1)
    ridgelines.update_traces(showlegend=True, row=1, col=1)

    # Move boxplot legend.
    boxes.update_layout(legend=dict(
        orientation="h", yanchor="bottom", title_text="TCGA Project:", x=0.5, xanchor="center",
        ))

    # Add titles.
    boxes.update_layout(
        title=f"TCGA {seq_type.upper()}<br>"
              f"Asymmetry by {bias_type} strand<br>"
              f"Call vs. complement {'log<sub>2</sub>' if log_transform_ratio else ''}ratio",
        )
    scatters.update_layout(title=f"TCGA {seq_type}<br>"
                                 f"Asymmetry by {bias_type} strand<br>"
                                 f"Call vs. complement counts")

    # Remove box plot tick labels.
    boxes.update_xaxes(showticklabels=False)

    ridgelines.update_traces(orientation="h", side="positive", width=3, points=False)
    ridgelines.update_yaxes(showticklabels=False)
    ridgelines.update_legends(traceorder="reversed")
    ridgelines.update_layout(title=boxes.layout.title, legend_title="TCGA Project:")
    ridgelines.data = ridgelines.data[::-1]

    return boxes, ridgelines, scatters


def plot_oxog_calls_vs_mismatches_per_project(mismatch_table, call_table,
                                              minimum_mismatch_count=None):
    minimums = _normalize_minimum(minimum_mismatch_count)
    _proj_numbers = {y: x for x, y in enumerate(mismatch_table.index.levels[0], start=-4)}
    plots = dict()
    mismatches = mismatch_table.mismatch_tools.calculate_reference_asymmetry("GT", "CA").droplevel("file_id")

    calls = filter_raw_call_data_dups(call_table)
    calls = calls.call_tools.calculate_reference_asymmetry("GT", "CA")

    for filter_status in ["TOTAL", "PASS"]:
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


def plot_heatmap_of_unique_site_counts(call_table):
    ref_table = log2(call_table.call_tools.count_unique_calls("reference", filters=["TOTAL", "PASS"]))
    trans_table = log2(call_table.call_tools.count_unique_calls("transcription", filters=["TOTAL", "PASS"]))
    plot = make_subplots(2, 2, vertical_spacing=0.05,
                         subplot_titles=["By reference strand<br>Total calls",
                                         "By transcription strand<br>Total calls",
                                         "Pass calls", "Pass calls"])
    # x_labels = ["C>G/G>C", "C>T/G>A", "G>T/C>A", "T>A/A>T", "T>C/A>G", "T>G/A>C"]
    x_labels = [r"$\frac{C➔G}{G➔C}$", r"$\frac{C➔T}{G➔A}$", r"$\frac{G➔T}{C➔A}$",
                r"$\frac{T➔A}{A➔T}$", r"$\frac{T➔C}{A➔G}$", r"$\frac{T➔G}{A➔C}$"]
    x_labels = ["C➔G<br>G➔C", "C➔T<br>G➔A", "G➔T<br>C➔A",
                "T➔A<br>A➔T", "T➔C<br>A➔G", "T➔G<br>A➔C"]
    y_labels = [x[5:] for x in ref_table.index.values]
    plot.add_trace(go.Heatmap(y=y_labels, z=ref_table.TOTAL.ratio,
                              text=ref_table.TOTAL.ratio, texttemplate="%{text:.2f}",
                              colorscale="RdBu_r", zmin=-3, zmax=3,
                              colorbar_title="log2(ratio)"), row=1, col=1)
    plot.add_trace(go.Heatmap(z=trans_table.TOTAL.ratio,
                              text=trans_table.TOTAL.ratio, texttemplate="%{text:.2f}",
                              colorscale="RdBu_r", zmin=-3, zmax=3,
                              colorbar_title="log2(ratio)"), row=1, col=2)
    plot.add_trace(go.Heatmap(x=x_labels, y=y_labels, z=ref_table.PASS.ratio,
                              text=ref_table.PASS.ratio, texttemplate="%{text:.2f}",
                              colorscale="RdBu_r", zmin=-3, zmax=3,
                              colorbar_title="log2(ratio)"), row=2, col=1)
    plot.add_trace(go.Heatmap(x=x_labels, z=trans_table.PASS.ratio,
                              text=trans_table.PASS.ratio, texttemplate="%{text:.2f}",
                              colorscale="RdBu_r", zmin=-3, zmax=3,
                              colorbar_title="log2(ratio)"), row=2, col=2)
    plot.update_xaxes(showticklabels=False, row=1)
    plot.update_yaxes(showticklabels=False, col=2)
    plot.update_layout(title="Unique variant sites, change vs. complement")
    plot.update_yaxes(autorange="reversed")
    return plot


class PlotDirectory:
    def __init__(self, mismatch_table, variant_tables=None):
        self.mismatches = PlotsByBiasType(mismatch_table, "mismatches")
        if variant_tables:
            self.variants = PlotsByBiasType(variant_tables)
            self.correlation = CorrPlotsByFilterStatus(mismatch_table, variant_tables)

    def show_all(self):
        for plots in [self.mismatches, self.variants, self.correlation]:
            plots.show_all()

    def show_corr_plots(self):
        self.correlation.show_all()


class PlotsByBiasType:
    def __init__(self, data_table, data_type, inverse=False):
        params = dict(data_table=data_table, data_type=data_type, inverse=inverse)
        self.reference_bias = PlotsByPlotType(bias_type="reference", **params)
        self.transcription_bias = PlotsByPlotType(bias_type="transcription", **params)

    def show_all(self):
        for plots in [self.reference_bias, self.transcription_bias]:
            plots.show_all()


class PlotsByPlotType:
    def __init__(self, data_table, data_type, bias_type, inverse=False, filtered=False):
        if data_type not in {"variants", "mismatches"}:
            raise ValueError("'data_type' must be 'variants'"
                             f" or 'mismatches', not '{data_type}'.")
        self.box = data_table.__getattr__(attr_name)._plot_bias_box(**params)
        self.scatter = data_table.__getattr__(attr_name)._plot_bias_scatter(**params)
        self.hexbin = data_table.__getattr__(attr_name)._plot_bias_hexbin(**params)

    def show_all(self):
        self.box.show()
        self.scatter.show()
        if hasattr(self, "hexbin"):
            self.hexbin.show()


class CorrPlotsByFilterStatus:
    def __init__(self, mismatch_table, variant_tables):
        self.unfiltered = CorrPlotsByBiasType(mismatch_table, variant_tables.unfiltered, filtered=False)
        self.filtered = CorrPlotsByBiasType(mismatch_table, variant_tables.filtered, filtered=True)

    def show_all(self):
        for plots in [self.unfiltered, self.filtered]:
            plots.show_all()


class CorrPlotsByBiasType:
    def __init__(self, mismatch_table, variant_table, filtered):
        self.reference_bias = CorrPlots(mismatch_table, variant_table, "reference", filtered)
        self.transcription_bias = CorrPlots(mismatch_table, variant_table, "transcription", filtered)

    def show_all(self):
        for plots in [self.reference_bias, self.transcription_bias]:
            plots.show_all()


class CorrPlots:
    def __init__(self, mismatch_table, variant_table, bias_type, filtered):
        self.scatter = self.make_corr_scatter(mismatch_table, variant_table, bias_type, filtered)

    @staticmethod
    def make_corr_scatter(mismatch_table, variant_table, bias_type, filtered=False):
        mismatch_data = mismatch_table.ukb_mismatches.calculate_bias(bias_type, True, False)[["ratio"]].sort_index()
        variant_data = variant_table.ukb_variants.calculate_bias(bias_type, True).sort_index().Hets
        data = mismatch_data.join(variant_data, rsuffix="_var", lsuffix="_mis")
        corr_plot = go.Figure()
        for proj in data.index.unique("project_id"):
            sub_data = data.loc[proj]
            corr_plot.add_trace(go.Scatter(x=sub_data.ratio_mis, y=sub_data.ratio_var,
                                           mode="markers", name=proj))
        title = f"{bias_type.title()} strand 8-oxo-G mismatch vs. heterozygous variant ratio"
        if filtered:
            title += ", with common SNPs removed"
        corr_plot.update_layout(title=title,
                                xaxis_title="mismatch ratio",
                                yaxis_title="variant_ratio")
        return corr_plot

    def show_all(self):
        self.scatter.show()
