
import pandas as pd
from pandas.api.extensions import register_dataframe_accessor
from pandas import IndexSlice as idx
import plotly.graph_objects as go
import plotly.io as pio

from hexbin import plot_hexbin
from strand_bias.TCGA_analysis.capture_kit_counts import CaptureKit

pio.renderers.default = "browser"


# @register_dataframe_accessor("ukb_variants")
# class UkbVariantTable:
#     def __init__(self, pandas_obj):
#         self._obj = pandas_obj
#         self._has_mut_col = "mutation" in pandas_obj.index.names
#
#     @staticmethod
#     def read_csv(file_path):
#         table = pd.read_csv(file_path, sep="\t", header=0, index_col=[0, 1, 2, 3])
#         if all(pd.isna(table.index.unique("mutation"))):
#             table = table.droplevel("mutation")
#         # table.rename(index=lambda x: int(x.split("_")[0]), level=0, inplace=True)
#         return table
#
#     def unstack(self):
#         if self._has_mut_col:
#             table = self._obj.unstack([-2, -1])
#             table.columns.names = ["stat", "mutation", "transcribed"]
#         else:
#             table = self._obj.unstack()
#             table.columns.names = ["stat", "transcribed"]
#         return table


@register_dataframe_accessor("ukb_variants")
class UkbVariantsUnstacked:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        self._has_mut_header = "mutation" in pandas_obj.columns.names

    # def _unstack(self):
    #     if "mutation" in self._obj.index.names:
    #         table = self._obj.unstack([-2, -1])
    #         table.columns.names = ["stat", "mutation", "transcribed"]
    #     else:
    #         table = self._obj.unstack()
    #         table.columns.names = ["stat", "transcribed"]
    #     return table

    @staticmethod
    def read_csv(file_path):
        table = pd.read_csv(file_path, sep="\t", header=0, index_col=[0, 1, 2, 3])
        if all(pd.isna(table.index.unique("mutation"))):
            table = table.droplevel("mutation")
            table = table.unstack()
            table.columns.names = ["stat", "transcribed"]
        else:
            table = table.unstack([-2, -1])
            table.columns.names = ["stat", "mutation", "transcribed"]
        return table

    # @staticmethod
    # def read_csv(file_path):
    #     table = pd.read_csv(file_path, )
    #
    #     table = table.ukb_variants.unstack()
    #     return table

    def calculate_ref_bias(self, drop_counts=True,
                           normalization_data=None,
                           normalization_factor=1,
                           ratio_by_template=None):
        if not self._has_mut_header:
            return
        table = self._obj.groupby(["stat", "mutation"], axis=1).agg(sum)
        if normalization_data is not None:
            normalization_data = normalization_data.capkit.calculate_reference_counts()
        for stat in table.columns.levels[0]:
            if normalization_data is not None:
                table[(stat, "GT")] = (table[(stat, "GT")] / normalization_data.G) * normalization_factor
                table[(stat, "CA")] = (table[(stat, "CA")] / normalization_data.C) * normalization_factor
            table[(stat, "ratio")] = table[stat].GT / table[stat].CA
        if drop_counts:
            table = table.loc[:, idx[:, "ratio"]]
        return table

    def calculate_transcription_bias(self, drop_counts=True,
                                     normalization_data=None,
                                     normalization_factor=1,
                                     ratio_by_template=True):
        if not self._has_mut_header:
            return
        table = self._obj.copy()
        if normalization_data is not None:
            normalization_data = normalization_data.capkit.calculate_transcription_counts()
        for stat in table.columns.levels[0]:
            stat_table = table[stat]
            table[(stat, "GT", "combined")] = stat_table.GT.forward + stat_table.CA.reverse
            table[(stat, "CA", "combined")] = stat_table.CA.forward + stat_table.GT.reverse
            if normalization_data is not None:
                table[(stat, "GT", "combined")] = (table[(stat, "GT", "combined")]
                                                   / normalization_data.G) * normalization_factor
                table[(stat, "CA", "combined")] = (table[(stat, "CA", "combined")]
                                                   / normalization_data.C) * normalization_factor
            if ratio_by_template:
                table[(stat, "GT", "combined")], table[(stat, "CA", "combined")] = (table[(stat, "CA", "combined")],
                                                                                    table[(stat, "GT", "combined")])
            table[(stat, "ratio", "combined")] = (table[(stat, "GT", "combined")]
                                                  / table[(stat, "CA", "combined")])

        if drop_counts:
            table = table.loc[:, idx[:, "ratio", "combined"]]
        else:
            table = table.loc[:, idx[:, :, "combined"]]
        table = table.droplevel("transcribed", axis=1)
        return table

    def calculate_bias(self, bias_type, drop_counts=True, normalization_data=None,
                       normalization_factor=1, ratio_by_template=True):
        if bias_type == "transcription":
            return self.calculate_transcription_bias(
                drop_counts=drop_counts,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor,
                ratio_by_template=ratio_by_template
                )
        elif bias_type == "reference":
            return self.calculate_ref_bias(
                drop_counts=drop_counts,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor)
        else:
            raise ValueError(f"'bias_type' must be 'transcription' or 'reference'.")

    def plot_ref_bias_histogram(self, normalization_data=None, normalization_factor=1):
        ref_bias = self.calculate_ref_bias(normalization_data=normalization_data,
                                           normalization_factor=normalization_factor)
        fig = go.Figure(go.Histogram(x=ref_bias.Hets.ratio))
        return fig

    def plot_transcription_bias_histogram(self, normalization_data=None, normalization_factor=1,
                                          ratio_by_template=True):
        transcription_bias = self.calculate_transcription_bias(
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        fig = go.Figure(go.Histogram(x=transcription_bias.Hets.ratio))
        return fig

    def sort_by_batch_ratio_median(self, bias_type, drop_counts=True, normalization_data=None,
                                   normalization_factor=1, ratio_by_template=True):
        bias_data = self._obj.ukb_variants.calculate_bias(
            bias_type=bias_type,
            drop_counts=drop_counts,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        medians = {batch: bias_data.loc[batch].Hets.ratio.median()
                   for batch in bias_data.index.unique("project_id")}
        bias_data.sort_values(by=("Hets", "ratio"), inplace=True)
        bias_data.sort_index(level=0, key=lambda col: [medians[batch] for batch in col],
                             inplace=True)
        return bias_data

    def _plot_bias_box(self, bias_type, sort_plot=True, inverse=False, filtered=False,
                       normalization_data=None, normalization_factor=1, ratio_by_template=True, **kwargs):
        plot = go.Figure()
        if sort_plot:
            bias_data = self._obj.ukb_variants.sort_by_batch_ratio_median(
                bias_type=bias_type,
                drop_counts=True,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor,
                ratio_by_template=ratio_by_template
                )
        else:
            bias_data = self._obj.ukb_variants.calculate_bias(
                bias_type=bias_type,
                drop_counts=True,
                normalization_data=normalization_data,
                normalization_factor=normalization_factor,
                ratio_by_template=ratio_by_template)

        if inverse:
            bias_data = 1/bias_data
        title = f"{bias_type.title()} strand 8-oxo-G heterozygous variant bias by flowcell"
        if filtered:
            title += ", with common SNPs removed"
        yaxis_title = "G>T / C>A"
        if inverse:
            yaxis_title = "C>A / G>T"
        if bias_type == "transcription":
            yaxis_title = "coding strand " + yaxis_title
        for batch in bias_data.index.unique("project_id"):
            batch_subset = bias_data.loc[batch]
            plot.add_trace(go.Box(y=batch_subset.Hets.ratio, name=batch,
                                  marker=dict(color="blue", opacity=0.2),
                                  line=dict(color="rgba(0, 0, 255, 0.75)"),
                                  showlegend=False))

        median_x, median_y = zip(*[(batch, bias_data.loc[batch].Hets.ratio.median())
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

    @staticmethod
    def _make_scatter_titles(x, y, bias_type, filtered):
        title = f"{bias_type.title()} strand 8-oxo-G heterozygous variant bias"
        if filtered:
            title += ", with common SNPs removed"
        xaxis_title = ">".join(x) + " variants"
        yaxis_title = ">".join(y) + " variants"
        titles = dict(title=title, xaxis_title=xaxis_title, yaxis_title=yaxis_title)
        return titles

    def _plot_bias_scatter(self, bias_type, inverse=False, filtered=False,
                           normalization_data=None, normalization_factor=1, ratio_by_template=True, **kwargs):
        fig = go.Figure()
        data = self.calculate_bias(
            bias_type=bias_type,
            drop_counts=False,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        xcol, ycol = ("CA", "GT") if not inverse else ("GT", "CA")
        for proj in self._obj.index.unique("project_id"):
            fig.add_trace(go.Scatter(x=data.loc[idx[proj], idx["Hets", xcol]],
                                     y=data.loc[idx[proj], idx["Hets", ycol]],
                                     mode="markers"))
        maximum = data.loc[:, idx["Hets", ["GT", "CA"]]].max().max()
        fig.add_trace(go.Scatter(x=[0, maximum], y=[0, maximum], mode="lines",
                                 line=dict(color="lightgray"), name="y = x"))
        fig.update_layout(**self._make_scatter_titles(xcol, ycol, bias_type, filtered))
        return fig

    def _plot_bias_hexbin(self, bias_type, inverse=False, filtered=False,
                          normalization_data=None, normalization_factor=1,
                          ratio_by_template=True, showlegend=True,
                          return_graph_objects=False, **kwargs):
        data = self.calculate_bias(
            bias_type=bias_type,
            drop_counts=False,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )
        xcol, ycol = ("CA", "GT") if not inverse else ("GT", "CA")
        fig = plot_hexbin(data.Hets[xcol], data.Hets[ycol],
                          showlegend=showlegend,
                          return_graph_objects=return_graph_objects)
        if return_graph_objects:
            gos = fig
        maximum = data.loc[:, idx["Hets", ["GT", "CA"]]].max().max()
        minimum = data.loc[:, idx["Hets", ["GT", "CA"]]].min().min()
        yx_line = go.Scatter(x=[minimum, maximum], y=[minimum, maximum], mode="lines",
                             line=dict(color="lightgray", dash="dash"), name="y = x",
                             showlegend=False)
        if return_graph_objects:
            gos.append(yx_line)
            return gos
        fig.add_trace(yx_line)
        titles = self._make_scatter_titles(xcol, ycol, bias_type, filtered)
        titles["title"] += ", hexbinned"
        fig.update_layout(**titles)
        return fig

    def plot_transcription_bias_scatter(self, normalization_data=None, normalization_factor=1, ratio_by_template=True):
        return self._plot_bias_scatter(
            bias_type="transcription",
            normalization_data=normalization_data,
            normalization_factor=normalization_factor,
            ratio_by_template=ratio_by_template
            )

    def plot_reference_bias_scatter(self, normalization_data=None, normalization_factor=1):
        return self._plot_bias_scatter(
            bias_type="reference",
            normalization_data=normalization_data,
            normalization_factor=normalization_factor
            )

    def calculate_aggregate_bias(self, bias_type, normalization_data=None,
                                 normalization_factor=1, ratio_by_template=True):
        table = self.calculate_bias(
            bias_type=bias_type,
            drop_counts=False,
            normalization_data=normalization_data,
            normalization_factor=normalization_factor
            )
        table = table.loc[:, idx[:, ["GT", "CA"]]].sum()
        for stat in table.index.unique("stat"):
            if ratio_by_template:
                table[stat, "ratio"] = table[stat, "CA"] / table[stat, "GT"]
            else:
                table[stat, "ratio"] = table[stat, "GT"] / table[stat, "CA"]
        table = table.sort_index()
        return table

if __name__ == '__main__':
    xgen = CaptureKit.read_capture_kit_nucleotide_summary(
        "/home/tyler/Documents/Resource_Data/capture_kits/IDT_xGen_Exome_Hyb_Panel/nt_counts.tsv")
    variants = UkbVariantsUnstacked.read_csv("~/StrandBias/UKB_analysis/variant_reanalysis/oxog_vars.snp_filtered.stats.project_ids.tsv")
