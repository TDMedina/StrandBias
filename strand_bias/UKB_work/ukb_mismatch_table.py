
import pandas as pd
from pandas import IndexSlice as idx
from pandas.api.extensions import register_dataframe_accessor
import plotly.graph_objects as go

from hexbin import plot_hexbin

from strand_bias.TCGA_analysis.capture_kit_counts import CaptureKit

xgen = CaptureKit.read_capture_kit_nucleotide_summary("/home/tyler/Documents/Resource_Data/capture_kits/IDT_xGen_Exome_Hyb_Panel/nt_counts.tsv")


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

    def reference_bias_counts(self, add_ratio=True, add_fraction=False, normalization_data=None):
        ref_bias = self._obj.groupby("mutation", axis=1).agg(sum)
        if normalization_data is not None:
            norm_counts = normalization_data.capkit.calculate_reference_counts()
            ref_bias["GT"] = ref_bias.GT * norm_counts.G / norm_counts.sum()
            ref_bias["CA"] = ref_bias.CA * norm_counts.C / norm_counts.sum()
        if add_ratio:
            ref_bias["ratio"] = ref_bias.GT / ref_bias.CA
        if add_fraction:
            ref_bias["fraction"] = ref_bias.GT / ref_bias.sum(axis=1)
        return ref_bias

    def transcription_bias_counts(self, add_ratio=True, add_fraction=False, normalization_data=None):
        gt = self._obj.forward.GT + self._obj.reverse.CA
        gt.name = "GT"
        ca = self._obj.forward.CA + self._obj.reverse.GT
        ca.name = "CA"
        trans_bias = pd.concat([gt, ca], axis=1)
        if normalization_data is not None:
            norm_counts = normalization_data.capkit.calculate_transcription_counts()
            trans_bias["GT"] = trans_bias.GT / (norm_counts.G / norm_counts.sum())
            trans_bias["CA"] = trans_bias.CA / (norm_counts.C / norm_counts.sum())
        if add_ratio:
            trans_bias["ratio"] = trans_bias.GT / trans_bias.CA
        if add_fraction:
            trans_bias["fraction"] = trans_bias.GT / trans_bias.sum(axis=1)
        return trans_bias

    def calculate_bias(self, bias_type, add_ratio=True, add_fraction=True, normalization_data=None):
        if bias_type == "transcription":
            bias_data = self._obj.ukb_mismatches.transcription_bias_counts(add_ratio, add_fraction, normalization_data)
        elif bias_type == "reference":
            bias_data = self._obj.ukb_mismatches.reference_bias_counts(add_ratio, add_fraction, normalization_data)
        else:
            raise ValueError(f"Invalid bias_type '{bias_type}'. bias_type "
                             f"must be one of 'transcription' or 'reference'.")
        return bias_data

    def sort_by_batch_ratio_median(self, bias_type, normalization_data=None):
        bias_data = self._obj.ukb_mismatches.calculate_bias(bias_type, True, False, normalization_data)
        medians = {batch: bias_data.loc[batch].ratio.median()
                   for batch in bias_data.index.unique("project_id")}
        bias_data.sort_values(by="ratio", inplace=True)
        bias_data.sort_index(level=0, key=lambda col: [medians[batch] for batch in col],
                             inplace=True)
        return bias_data

    def sort_by_case_ratio(self, bias_type):
        bias_data = self._obj.ukb_mismatches.calculate_bias(bias_type, True, False)
    def sort_by_case_ratio(self, bias_type, normalization_data=None):
        bias_data = self._obj.ukb_mismatches.calculate_bias(bias_type, True, False, normalization_data)
        bias_data.sort_values(by="ratio", inplace=True)
        return bias_data

    @staticmethod
    def _make_scatter_titles(x, y, bias_type):
        title = f"{bias_type.title()} strand 8-oxo-G mismatch bias"
        xaxis_title = ">".join(x) + " mismatches"
        yaxis_title = ">".join(y) + " mismatches"
        titles = dict(title=title, xaxis_title=xaxis_title, yaxis_title=yaxis_title)
        return titles

    def _plot_bias_scatter(self, bias_type, normalization_data=None, inverse=False, **kwargs):
        """Bias type must be one of 'transcription' or 'reference'."""
        plot = go.Figure()
        bias_data = self._obj.ukb_mismatches.calculate_bias(bias_type, False, False, normalization_data)
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

    def _plot_bias_box(self, bias_type, sort_plot=True, normalization_data=None, inverse=False, **kwargs):
        plot = go.Figure()
        if sort_plot:
            bias_data = self._obj.ukb_mismatches.sort_by_batch_ratio_median(bias_type, normalization_data=normalization_data)
        else:
            bias_data = self._obj.ukb_mismatches.calculate_bias(bias_type, True, False, normalization_data)
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

    def _plot_bias_hexbin(self, bias_type, inverse=False, normalization_data=None, **kwargs):
        data = self._obj.ukb_mismatches.calculate_bias(bias_type, normalization_data=normalization_data)
        xcol, ycol = ("CA", "GT") if not inverse else ("GT", "CA")
        fig = plot_hexbin(data[xcol], data[ycol])
        titles = self._make_scatter_titles(xcol, ycol, bias_type)
        titles["title"] += ", hexbinned"
        fig.update_layout(**titles)
        return fig

    def plot_reference_bias_scatter(self, normalization_data=None):
        return self._plot_bias_scatter("reference", normalization_data=normalization_data)

    def plot_transcription_bias_scatter(self,normalization_data=None):
        return self._plot_bias_scatter("transcription", normalization_data=normalization_data)

    def plot_reference_bias_box(self, inverse=False, normalization_data=None):
        return self._plot_bias_box("reference", inverse=inverse, normalization_data=normalization_data)

    def plot_transcription_bias_box(self, inverse=False, normalization_data=None):
        return self._plot_bias_box("transcription", inverse=inverse, normalization_data=normalization_data)

    def calculate_aggregate_bias(self, bias_type, normalization_data=None):
        table = self.calculate_bias(bias_type, False, False)
        table = table.sum()
        if normalization_data is not None:
            norm_counts = normalization_data.capkit.calculate_counts(bias_type)
            table["GT"] = table.GT / (norm_counts.G / norm_counts.sum())
            table["CA"] = table.CA / (norm_counts.C / norm_counts.sum())
        table["ratio"] = table.GT / table.CA
        return table
