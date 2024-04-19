
from numpy import log2
import pandas as pd
from pandas import IndexSlice as idx
from pandas.api.extensions import register_dataframe_accessor
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"


# @register_dataframe_accessor("pileup_tools")
# class PileupTable:
#     def __init__(self, pandas_obj):
#         self._obj = pandas_obj
#
#     @staticmethod
#     def read_csv(file_path):
#         table = pd.read_csv(file_path, sep="\t", index_col=[0, 1, 2], header=[0, 1])
#         return table
#
#     def make_asymmetry_summary_table(self, by="coding_strand",
#                                      as_proportion=False, as_ratio=False):
#         summary_table = (self._obj
#                          .groupby(["reference", by]).agg(sum)
#                          .groupby(axis=1, level="alt").agg(sum))
#         if as_proportion:
#             summary_table = summary_table.groupby("reference").agg(self._make_summary_proportion)
#         elif as_ratio:
#             summary_table = summary_table.groupby("reference").agg(self._make_summary_ratio)
#         summary_table = summary_table.loc[["C", "G"], ["A", "T"]]
#         return summary_table
#
#     def calculate_orientation_bias_by_coding_strand(self, drop_non_ox=True):
#         results = self._obj.groupby("alt", axis=1).agg(sum)
#         results = results.loc[idx[:, :, "F1R2"]] / results.loc[idx[:, :, "F2R1"]]
#         if drop_non_ox:
#             results = results.loc[idx[["C", "G"], :], ["A", "T"]]
#         return results
#
#     def calculate_strand_bias_by_coding_strand(self, drop_non_ox=True):
#         results = self._obj.groupby(["reference", "coding_strand"]).agg(sum)
#         results = results.groupby("alt", axis=1).agg(self._agg_div_alignment)
#         if drop_non_ox:
#             results = results.loc[idx[["C", "G"], :], ["A", "T"]]
#         return results
#
#     @staticmethod
#     def _agg_div_alignment(df):
#         df = df.droplevel("alt", axis=1)
#         return df.forward / df.reverse
#
#     @staticmethod
#     def _make_summary_proportion(pair):
#         if not len(pair) > 1 or min(pair) == 0:
#             return pd.NA
#         norm = min(pair)
#         pair = ["1" if val == norm else f"{val/norm:.2f}" for val in pair]
#         pair = f"{pair[0]}:{pair[1]}"
#         return pair
#
#     @staticmethod
#     def _make_summary_ratio(pair):
#         if not len(pair) > 1 or min(pair) == 0:
#             return pd.NA
#         pair = pair[0] / pair[1]
#         return pair


@register_dataframe_accessor("mismatch_tools")
class ConcatenatedPileupTable:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(file_path):
        table = pd.read_csv(file_path, sep="\t", index_col=list(range(6)), header=[0, 1])
        return table

    def simplify_for_reference_asymmetry(self, change_numerator, change_denominator,
                                         add_ratio_column=True, log_transform_ratio=False,
                                         add_fraction_column=True,
                                         normalization_counts=None, normalization_factor=1):
        table = (self._obj
                 .groupby(["project_id", "case_id", "file_id", "reference"])
                 .agg(sum)
                 .groupby("alt", axis=1)
                 .agg(sum))
        table = table.unstack(level=-1)[[tuple(reversed(change_numerator)),
                                         tuple(reversed(change_denominator))]]
        table.columns = [change_numerator, change_denominator]

        if normalization_counts is not None:
            for col in table.columns:
                table[col] = table[col] / normalization_counts[col[0]] * normalization_factor

        if add_ratio_column:
            ratio_label = f"{change_numerator}{change_denominator}_ratio".lower()
            ratios = table[change_numerator] / table[change_denominator]
            if log_transform_ratio:
                ratios = log2(ratios)
            table[ratio_label] = ratios
        if add_fraction_column:
            ratio_label = f"{change_numerator}{change_denominator}_fraction".lower()
            ratios = (table[change_numerator]
                      / (table[change_denominator] + table[change_numerator]))
            table[ratio_label] = ratios
        return table

    def subset_oxo_nucleotides(self):
        table = self._obj.loc[idx[:, ["C", "G"], :, :], ["A", "T"]]
        return table

    def collapse_coding_regions(self):
        index_levels = list(self._obj.index.names)
        index_levels.remove("coding_strand")
        table = self._obj.groupby(index_levels).agg(sum)
        return table

    def collapse_alignment_direction_columns(self):
        table = self._obj.groupby("alt", axis=1).agg(sum)
        return table

    def calculate_ref_strand_ox_bias(self):
        results = (self._obj.loc[idx[:, ["C", "G"], :, :], ["A", "T"]]
                   .groupby("alt", axis=1).agg(sum)
                   .groupby(["file_id", "reference"]).agg(sum))
        results = results.loc[idx[:, "G"], "T"].droplevel("reference") / results.loc[idx[:, "C"], "A"].droplevel("reference")
        return results

    def plot_ref_strand_ox_scatter(self):
        results = (self._obj.loc[idx[:, ["C", "G"], :, :], ["A", "T"]]
                   .groupby("alt", axis=1).agg(sum)
                   .groupby(["file_id", "reference"]).agg(sum))
        x = results.loc[idx[:, "G"], "T"]
        y = results.loc[idx[:, "C"], "A"]
        plot = go.Figure(go.Scatter(x=x, y=y, mode="markers", name="Samples"))
        plot.add_trace(go.Scatter(x=[0, max(x)], y=[0, max(x)], mode="lines", name="y = x"))
        plot.update_xaxes(title="G->T mismatch count")
        plot.update_yaxes(title="C->A mismatch count")
        return plot

    def plot_ox_nt_ori_boxplots(self):
        plot = go.Figure()
        for base, mut in (("C", "A"), ("G", "T")):
            for ori in ("F1R2", "F2R1"):
                data = self._obj.loc[idx[:, base, ori], mut]
                plot.add_trace(go.Box(y=data, name=f"{base}➔{mut}, {ori}"))
        return plot
