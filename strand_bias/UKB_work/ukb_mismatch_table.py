
import pandas as pd
from numpy import log2
from numpy.random import normal
from pandas import IndexSlice as idx
from pandas.api.extensions import register_dataframe_accessor
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from hexbin import plot_hexbin

from strand_bias.TCGA_analysis.capture_kit_counts import CaptureKit
from strand_bias.TCGA_analysis.utilities import rectify_ratios

pio.renderers.default = "browser"


def filter_by_ancestry(data_table, ancestry_file):
    ancestry = pd.read_csv(ancestry_file, sep="\t", index_col=[0, 1])
    ancestry = ancestry.loc[~ pd.isna(ancestry.British)]
    ancestry = ancestry.loc[ancestry.British & (ancestry.ancestry_codes == "1001")]
    ancestry.index = ancestry.index.remove_unused_levels()
    british = data_table.loc[list(data_table.reset_index("sample_id").sample_id.isin(ancestry.index.levels[1]))]
    british.sort_index(inplace=True)
    return british


@register_dataframe_accessor("ukb_mismatches")
class UkbMismatchTable:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(file_path):
        asym = pd.read_csv(file_path, sep="\t", index_col=[0, 1], header=[0, 1])
        return asym

    def reference_bias_counts(self, add_ratio=True, add_fraction=False,
                              normalization_data=None, normalization_factor=1,
                              ratio_by_template=None):
        ref_bias = self._obj.groupby("mutation", axis=1).agg(sum)
        if normalization_data is not None:
            norm_counts = normalization_data.capkit.calculate_reference_counts()
            ref_bias["GT"] = (ref_bias.GT / norm_counts.G) * normalization_factor
            ref_bias["CA"] = (ref_bias.CA / norm_counts.C) * normalization_factor
        if add_ratio:
            ref_bias["ratio"] = ref_bias.GT / ref_bias.CA
        if add_fraction:
            ref_bias["fraction"] = ref_bias.GT / ref_bias.sum(axis=1)
        return ref_bias

    def transcription_bias_counts(self, add_ratio=True, add_fraction=False,
                                  normalization_data=None, normalization_factor=1,
                                  ratio_by_template=True):
        gt = self._obj.forward.GT + self._obj.reverse.CA
        gt.name = "GT"
        ca = self._obj.forward.CA + self._obj.reverse.GT
        ca.name = "CA"
        trans_bias = pd.concat([gt, ca], axis=1)
        if normalization_data is not None:
            norm_counts = normalization_data.capkit.calculate_transcription_counts()
            trans_bias["GT"] = (trans_bias.GT / norm_counts.G) * normalization_factor
            trans_bias["CA"] = (trans_bias.CA / norm_counts.C) * normalization_factor
        if ratio_by_template:
            trans_bias["GT"], trans_bias["CA"] = trans_bias.CA, trans_bias.GT
        if add_ratio:
            trans_bias["ratio"] = trans_bias.GT / trans_bias.CA
            # if ratio_by_template:
            #     trans_bias["ratio"] = trans_bias.CA / trans_bias.GT
            # else:
            #     trans_bias["ratio"] = trans_bias.GT / trans_bias.CA
        if add_fraction:
            trans_bias["fraction"] = trans_bias.GT / trans_bias.sum(axis=1)
        return trans_bias

    def calculate_bias(self, bias_type, add_ratio=True, add_fraction=True,
                       normalization_data=None, normalization_factor=1,
                       ratio_by_template=True):
        if bias_type == "transcription":
            bias_data = self._obj.ukb_mismatches.transcription_bias_counts(
                add_ratio=add_ratio,
                add_fraction=add_fraction,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor,
                ratio_by_template=ratio_by_template
                )
        elif bias_type == "reference":
            bias_data = self._obj.ukb_mismatches.reference_bias_counts(
                add_ratio=add_ratio,
                add_fraction=add_fraction,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor)
        else:
            raise ValueError(f"Invalid bias_type '{bias_type}'. bias_type "
                             f"must be one of 'transcription' or 'reference'.")
        return bias_data

    def sort_by_batch_ratio_median(self, bias_type, normalization_data=None,
                                   normalization_factor=1, ratio_by_template=True):
        bias_data = self._obj.ukb_mismatches.calculate_bias(
            bias_type=bias_type,
            add_ratio=True,
            add_fraction=False,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        medians = {batch: bias_data.loc[batch].ratio.median()
                   for batch in bias_data.index.unique("project_id")}
        bias_data.sort_values(by="ratio", inplace=True)
        bias_data.sort_index(level=0, key=lambda col: [medians[batch] for batch in col],
                             inplace=True)
        return bias_data

    def sort_by_case_ratio(self, bias_type, normalization_data=None, normalization_factor=1,
                           ratio_by_template=True):
        bias_data = self._obj.ukb_mismatches.calculate_bias(
            bias_type=bias_type,
            add_ratio=True,
            add_fraction=False,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        bias_data.sort_values(by="ratio", inplace=True)
        return bias_data

    @staticmethod
    def _make_scatter_titles(x, y, bias_type):
        title = f"{bias_type.title()} strand 8-oxo-G mismatch bias"
        xaxis_title = ">".join(x) + " mismatches"
        yaxis_title = ">".join(y) + " mismatches"
        titles = dict(title=title, xaxis_title=xaxis_title, yaxis_title=yaxis_title)
        return titles

    def _plot_bias_scatter(self, bias_type, normalization_data=None, normalization_factor=1,
                           inverse=False, ratio_by_template=True, **kwargs):
        """Bias type must be one of 'transcription' or 'reference'."""
        plot = go.Figure()
        bias_data = self._obj.ukb_mismatches.calculate_bias(
            bias_type=bias_type,
            add_ratio=False,
            add_fraction=False,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        xcol, ycol = ("CA", "GT") if not inverse else ("GT", "CA")
        for batch in bias_data.index.unique("project_id"):
            batch_subset = bias_data.loc[batch]
            plot.add_trace(go.Scatter(x=batch_subset[xcol], y=batch_subset[ycol],
                                      name=batch, mode="markers"))
        max_point = max(bias_data.CA.max(), bias_data.GT.max())
        plot.add_trace(go.Scatter(x=[0, max_point], y=[0, max_point], name="y = x",
                                  mode="lines", line=dict(color="lightgray")))
        titles = self._make_scatter_titles(xcol, ycol, bias_type)
        plot.update_layout(**titles)
        return plot


    def _plot_bias_box(self, bias_type, sort_plot=True, normalization_data=None,
                       normalization_factor=1, inverse=False,
                       ratio_by_template=True, **kwargs):
        plot = go.Figure()
        if sort_plot:
            bias_data = self._obj.ukb_mismatches.sort_by_batch_ratio_median(
                bias_type=bias_type,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor,
                ratio_by_template=ratio_by_template
                )
        else:
            bias_data = self._obj.ukb_mismatches.calculate_bias(
                bias_type=bias_type,
                add_ratio=True,
                add_fraction=False,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor,
                ratio_by_template=ratio_by_template
                )
        yaxis_title = "G>T / C>A"
        if inverse:
            yaxis_title = "C>A / G>T"
            bias_data["ratio"] = 1/bias_data.ratio
        title = f"{bias_type.title()} strand 8-oxo-G mismatch bias by flowcell"
        for batch in bias_data.index.unique("project_id"):
            batch_subset = bias_data.loc[batch]
            plot.add_trace(go.Box(y=batch_subset.ratio, name=batch,
                                  marker=dict(color="blue", opacity=0.2),
                                  line=dict(color="rgba(0, 0, 255, 0.75)"),
                                  showlegend=False))

        median_x, median_y = zip(*[(batch, bias_data.loc[batch].ratio.median())
                                   for batch in bias_data.index.unique("project_id")])
        plot.add_trace(go.Scatter(x=median_x, y=median_y, marker=dict(color="cyan"),
                                  mode="markers", name="Medians"))

        plot.update_traces(marker=dict(size=3))

        plot.layout.xaxis2 = go.layout.XAxis(overlaying="x", range=[0, 2],
                                             showticklabels=False)
        plot.add_scatter(x=[0, 2], y=[1, 1], mode="lines", xaxis="x2", showlegend=False,
                         line=dict(dash="dash", color="firebrick", width=2))
        plot.update_layout(title=title,
                           xaxis_title=f"UKB Flowcell (n={len(bias_data.index.unique('project_id'))})",
                           yaxis_title=yaxis_title)
        plot.update_xaxes(tickfont=dict(color="rgba(0,0,0,0)", size=1))
        return plot

    def _plot_bias_hexbin(self, bias_type, inverse=False, normalization_data=None,
                          normalization_factor=1, ratio_by_template=True,
                          return_graph_objects=False, **kwargs):
        data = self._obj.ukb_mismatches.calculate_bias(
            bias_type=bias_type,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        xcol, ycol = ("CA", "GT") if not inverse else ("GT", "CA")
        fig = plot_hexbin(data[xcol], data[ycol], return_graph_objects=return_graph_objects)
        if return_graph_objects:
            return fig
        titles = self._make_scatter_titles(xcol, ycol, bias_type)
        titles["title"] += ", hexbinned"
        fig.update_layout(**titles)
        return fig

    def plot_reference_bias_scatter(self, normalization_data=None, normalization_factor=1):
        return self._plot_bias_scatter("reference", normalization_data=normalization_data,
                                       normalization_factor=normalization_factor)

    def plot_transcription_bias_scatter(self,normalization_data=None, normalization_factor=1, ratio_by_template=True):
        return self._plot_bias_scatter("transcription", normalization_data=normalization_data,
                                       normalization_factor=normalization_factor,
                                       ratio_by_template=ratio_by_template)

    def plot_reference_bias_box(self, inverse=False, normalization_data=None, normalization_factor=1):
        return self._plot_bias_box("reference", inverse=inverse, normalization_data=normalization_data,
                                   normalization_factor=normalization_factor)

    def plot_transcription_bias_box(self, inverse=False, normalization_data=None,
                                    normalization_factor=1, ratio_by_template=True):
        return self._plot_bias_box("transcription", inverse=inverse, normalization_data=normalization_data,
                                   normalization_factor=normalization_factor,
                                   ratio_by_template=ratio_by_template)

    def calculate_aggregate_bias(self, bias_type, normalization_data=None,
                                 normalization_factor=1, ratio_by_template=True):
        table = self.calculate_bias(bias_type, False, False)
        table = table.sum()
        if normalization_data is not None:
            norm_counts = normalization_data.capkit.calculate_counts(bias_type)
            table["GT"] = (table.GT / norm_counts.G) * normalization_factor
            table["CA"] = (table.CA / norm_counts.C) * normalization_factor
        if ratio_by_template:
            table["ratio"] = table.CA / table.GT
        else:
            table["ratio"] = table.GT / table.CA
        return table


    def plot_combined_bias_figure(self, normalization_data=None,
                                  normalization_factor=1,
                                  log_transform_boxplots=False,
                                  ratio_by_template=True,
                                  manual_yx_line=None,
                                  highlight_top10=False,
                                  add_subplot_titles=True):
        subplot_titles = None
        if add_subplot_titles:
            subplot_titles = [
                "Reference strand asymmetry:<br>Mismatch counts per 1k nts",
                "Transcription strand asymmetry:<br>Mismatch counts per 1k nts",
                "Reference strand asymmetry:<br>Mismatch ratio by flowcell",
                "Transcription strand asymmetry:<br>Mismatch ratio by flowcell",
                ]
        else:
            subplot_titles = [f"<b>{label+')':<80}</b>" for label in "abcd"]
        plot = make_subplots(2, 2, vertical_spacing=0.1, subplot_titles=subplot_titles)
        ref_data = self._obj.ukb_mismatches.sort_by_batch_ratio_median(
            bias_type="reference",
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        trans_data = self._obj.ukb_mismatches.sort_by_batch_ratio_median(
            bias_type="transcription",
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        for i, data in ((1, ref_data), (2, trans_data)):
            # Prepare data for boxplots.
            if log_transform_boxplots:
                data["ratio"] = log2(data.ratio)
            batches = data.index.unique("project_id")
            top10_batches = set()
            if highlight_top10 and i == 2:
                top10_batches = set(batches[-(len(batches) // 10):])
                # projects = projects[:cutoff]

                # proj_ids = data.index.unique("project_id")
                # top10_proj = proj_ids[-(len(proj_ids) // 10):]
                # non_top10_proj = list(set(proj_ids)-set(top10_proj))
                # data = data.loc[non_top10_proj]

            medians = []
            # for batch in data.index.unique("project_id"):
            # Plot boxplot per batch.
            for batch in batches:
                batch_data = data.loc[batch]
                medians.append((batch, batch_data.ratio.median()))

                # Set marker color for top batches.
                if highlight_top10 and batch in top10_batches:
                    color = "rgba(255, 0, 0, 0.75)"
                else:
                    color = "rgba(0, 0, 255, 0.75)"

                # Plot upper scatterplot.
                plot.add_trace(go.Scatter(x=batch_data["CA"], y=batch_data["GT"],
                                          mode="markers", marker_color=color,
                                          showlegend=False),
                               row=1, col=i)

                plot.add_trace(go.Box(y=batch_data.ratio, name=batch,
                                      marker=dict(color="grey", opacity=0.4, size=4),
                                      line=dict(color=color), showlegend=False),
                               row=2, col=i)

            # Add y=x line for scatterplots.
            if manual_yx_line:
                minmax = manual_yx_line
            else:
                minmax = (data[["CA", "GT"]].min().min(), data[["CA", "GT"]].max().max())
            plot.add_trace(go.Scatter(x=minmax, y=minmax, mode="lines",
                                      line=dict(color="rgb(211, 211, 211, 0.25)", dash="dash"),
                                      showlegend=bool(i%2), name="y=x"),
                           row=1, col=i)

            median_x, median_y = zip(*medians)

            # Plot median values onto boxplots.
            plot.add_trace(go.Scatter(x=median_x, y=median_y,
                                      marker=dict(color="cyan", size=4),
                                      mode="markers", name="Medians",
                                      showlegend=False),
                           row=2, col=i)

        plot.update_xaxes(showticklabels=False, row=2)
        # scatter_x_title = r"$\huge{\frac{C➔A}{1k\;C\:nts}}$"
        scatter_x_title = "C➔A per 1000 C"
        # scatter_y_title = r"$\huge{\frac{G>T}{1k\;G\:nts}}$"
        scatter_y_title = "G➔T per 1000 G"
        box_x_title = "Flowcell"
        # box_y_title = r"$\huge{\log_{2}\frac{G>T}{C>A}}$"
        box_y_title = "log2(G➔T / C➔A)"
        plot.update_traces(showlegend=True, row=2, col=2, selector=dict(type="scatter"))
        plot.update_layout(
            title="UKB: G➔T vs. C➔A asymmetry",
            font_size=28,
            xaxis1_title=scatter_x_title, yaxis1_title=scatter_y_title,
            xaxis2_title=scatter_x_title, yaxis2_title=scatter_y_title,
            xaxis3_title=box_x_title, yaxis3_title=box_y_title,
            xaxis4_title=box_x_title, yaxis4_title=box_y_title,
            xaxis2_matches="x",
            yaxis2_matches="y",
            yaxis4_matches="y3",
            legend=dict(xanchor="right", orientation="v", x=1.15, yanchor="middle", y=0.5),
            # legend=dict(yanchor="bottom", y=-0.05, orientation="h", x=0.47),
            height=2400, width=2400, margin=dict(t=180, l=100)
            )
        plot.update_annotations(font_size=28)

        # plot.layout.xaxis5 = go.layout.XAxis(overlaying="x", range=[0, 2],
        #                                      showticklabels=False)
        # plot.add_scatter(x=[0, 2], y=[1, 1], mode="lines", xaxis="x5", showlegend=False,
        #                  line=dict(dash="dash", color="firebrick", width=2))
        # plot.update_layout(title=title,
        #                    xaxis_title=f"UKB Flowcell (n={len(bias_data.index.unique('project_id'))})",
        #                    yaxis_title=yaxis_title)
        # plot.update_xaxes(tickfont=dict(color="rgba(0,0,0,0)", size=1))
        return plot


if __name__ == '__main__':
    xgen = CaptureKit.read_capture_kit_nucleotide_summary(
        "/home/tyler/Documents/Resource_Data/capture_kits/IDT_xGen_Exome_Hyb_Panel/nt_counts.tsv")
    mismatches = UkbMismatchTable.read_csv("~/StrandBias/UKB_analysis/asymmetry.tsv")
