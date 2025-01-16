
from itertools import product
from statistics import mode

from matplotlib import colormaps
import pandas as pd
from pandas.api.extensions import register_dataframe_accessor
from pandas import CategoricalDtype
from pandas import IndexSlice as idx
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
from scipy import stats

from plotting import make_hexbin_data

pio.renderers.default = "browser"


# %% Hexbinning.

def plot_hexagons(hexbin_data, colormap="viridis", dims=(1200, 1200)):
    width, height = dims

    dx = mode(np.diff(sorted(set(hexbin_data.x))))
    dy = mode(np.diff(sorted(set(hexbin_data.y))))

    fig = go.Figure()
    viridis = colormaps[colormap]
    max_z = max([np.log2(z) for z in hexbin_data.counts])
    for x, y, z in zip(hexbin_data.x, hexbin_data.y, hexbin_data.counts):
        color = "rgba" + str(viridis(np.log2(z) / max_z))
        hex_points = [(x, y+2*dy/3), (x+dx, y+dy/3), (x+dx, y-dy/3), (x, y-2*dy/3),
                      (x-dx, y-dy/3), (x-dx, y+dy/3), (x, y+2*dy/3)]
        hex_points = list(zip(*hex_points))
        fig.add_trace(go.Scatter(x=hex_points[0], y=hex_points[1], mode="none",
                                 fill="toself", fillcolor=color, showlegend=False))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", showlegend=False,
                             marker=dict(colorbar=dict(title="Sample Count"),
                                         colorscale=colormap,
                                         showscale=True,
                                         cmin=1, cmax=max(hexbin_data.counts))))
    return fig


def plot_hexbin(x, y, colormap="viridis", dims=(1200, 1200)):
    hexbin_data = make_hexbin_data(x, y)
    fig = plot_hexagons(hexbin_data, colormap, dims)
    return fig


# %% The rest.
def filter_by_ancestry(data_table, ancestry_file):
    ancestry = pd.read_csv(ancestry_file, sep="\t", index_col=[0, 1])
    ancestry = ancestry.loc[~ pd.isna(ancestry.British)]
    ancestry = ancestry.loc[ancestry.British & (ancestry.ancestry_codes == "1001")]
    ancestry.index = ancestry.index.remove_unused_levels()
    british = data_table.loc[list(data_table.reset_index("case_id").case_id.isin(ancestry.index.levels[1]))]
    british.sort_index(inplace=True)
    return british


@register_dataframe_accessor("ukb")
class UkbTable:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(file_path, ancestry_file=None):
        asym = pd.read_csv(file_path, sep="\t", index_col=[0, 1], usecols=list(range(6)))
        asym.index.names = ["case_id", "project_id"]
        asym = asym.reset_index().set_index(["project_id", "case_id"])
        asym.columns = pd.MultiIndex.from_tuples([("forward", "GT"), ("reverse", "GT"),
                                                  ("forward", "CA"), ("reverse", "CA")])
        asym.columns.names = ["coding_region", "mutation"]
        asym = asym[sorted(asym.columns)]
        asym = asym.sort_index()
        if ancestry_file:
            asym = asym.ukb.filter_by_ancestry(ancestry_file)
        return asym

    def filter_by_ancestry(self, ancestry_file):
        british = filter_by_ancestry(self._obj, ancestry_file)
        return british

    def reference_bias_counts(self, add_ratio=True, add_fraction=True):
        ref_bias = self._obj.groupby("mutation", axis=1).agg(sum)
        if add_ratio:
            ref_bias["ratio"] = ref_bias.GT / ref_bias.CA
        if add_fraction:
            ref_bias["fraction"] = ref_bias.GT / ref_bias.sum(axis=1)
        return ref_bias

    def transcription_bias_counts(self, add_ratio=True, add_fraction=True):
        gt = self._obj.forward.GT + self._obj.reverse.CA
        gt.name = "GT"
        ca = self._obj.forward.CA + self._obj.reverse.GT
        ca.name = "CA"
        trans_bias = pd.concat([gt, ca], axis=1)
        if add_ratio:
            trans_bias["ratio"] = trans_bias.GT / trans_bias.CA
        if add_fraction:
            trans_bias["fraction"] = trans_bias.GT / trans_bias.sum(axis=1)
        return trans_bias

    def calculate_bias(self, bias_type, add_ratio=True, add_fraction=True):
        if bias_type == "transcription":
            bias_data = self._obj.ukb.transcription_bias_counts(add_ratio, add_fraction)
        elif bias_type == "reference":
            bias_data = self._obj.ukb.reference_bias_counts(add_ratio, add_fraction)
        else:
            raise ValueError(f"Invalid bias_type '{bias_type}'. bias_type "
                             f"must be one of 'transcription' or 'reference'.")
        return bias_data

    def sort_by_batch_ratio_median(self, bias_type):
        bias_data = self._obj.ukb.calculate_bias(bias_type, True, False)
        medians = {batch: bias_data.loc[batch].ratio.median()
                   for batch in bias_data.index.unique("project_id")}
        bias_data.sort_values(by="ratio", inplace=True)
        bias_data.sort_index(level=0, key=lambda col: [medians[batch] for batch in col],
                             inplace=True)
        return bias_data

    def sort_by_case_ratio(self, bias_type):
        bias_data = self._obj.ukb.calculate_bias(bias_type, True, False)
        bias_data.sort_values(by="ratio", inplace=True)
        return bias_data

    def _plot_bias_scatter(self, bias_type):
        """Bias type must be one of 'transcription' or 'reference'."""
        plot = go.Figure()
        bias_data = self._obj.ukb.calculate_bias(bias_type, False, False)
        title = f"{bias_type.title()} strand 8-oxo-G mismatch bias"
        for batch in bias_data.index.unique("project_id"):
            batch_subset = bias_data.loc[batch]
            plot.add_trace(go.Scatter(x=batch_subset.CA, y=batch_subset.GT,
                                      name=batch, mode="markers"))
        max_point = max(bias_data.CA.max(), bias_data.GT.max())
        plot.add_trace(go.Scatter(x=[0, max_point], y=[0, max_point], name="y = x",
                                  mode="lines", line=dict(color="lightgray")))
        plot.update_layout(title=title, xaxis_title="C>A mismatches",
                           yaxis_title="G>T mismatches")
        return plot

    def _plot_bias_box(self, bias_type, sort_plot=True, inverse=False):
        plot = go.Figure()
        if sort_plot:
            bias_data = self._obj.ukb.sort_by_batch_ratio_median(bias_type)
        else:
            bias_data = self._obj.ukb.calculate_bias(bias_type, True, False)
        if inverse:
            bias_data["ratio"] = 1/bias_data.ratio
        title = f"{bias_type.title()} strand 8-oxo-G mismatch bias"
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
        plot.update_layout(title=title)
        return plot

    def plot_reference_bias_scatter(self):
        return self._plot_bias_scatter("reference")

    def plot_transcription_bias_scatter(self):
        return self._plot_bias_scatter("transcription")

    def plot_reference_bias_box(self, inverse=False):
        return self._plot_bias_box("reference", inverse=inverse)

    def plot_transcription_bias_box(self, inverse=False):
        return self._plot_bias_box("transcription", inverse=inverse)


@register_dataframe_accessor("ukb_var")
class UkbVariants:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(variant_count_file):
        varcounts = pd.read_csv(variant_count_file, sep="\t",
                                usecols=[0, 1, 2, 3, 4, 5, 6, 7, 11])
        # cat_cols = ["coding_strand", "ref", "alt"]
        # cat_vals = [list("ACGT"), list("ACGT")]
        # for col, vals in zip(cat_cols, cat_vals):
        #     varcounts[col] = varcounts[col].astype(CategoricalDtype(vals))
        varcounts.set_index(["project_id", "case_id", "coding_strand", "ref", "alt"],
                            inplace=True)
        # varcounts = varcounts.groupby(["project_id", "case_id", "ref", "alt"]).agg(sum)
        varcounts.sort_index(inplace=True)
        return varcounts

    def filter_by_ancestry(self, ancestry_file):
        british = filter_by_ancestry(self._obj, ancestry_file)
        return british

    def subset_oxog(self):
        subset = self._obj.loc[idx[:, :, :, ["G"], ["T"]]]
        subset = pd.concat([subset, self._obj.loc[idx[:, :, :, ["C"], ["A"]]]])
        subset.sort_index(inplace=True)
        return subset

    def calculate_bias(self, bias_type, add_ratio=True, add_fraction=True):
        if bias_type == "transcription":
            bias_data = self._obj.ukb_var.transcription_bias_counts()
        elif bias_type == "reference":
            bias_data = self._obj.ukb_var.reference_bias_counts()
        else:
            raise ValueError(f"Invalid bias_type '{bias_type}'. bias_type "
                             f"must be one of 'transcription' or 'reference'.")
        return bias_data

    def transcription_bias_counts(self, add_ratio=True, add_fraction=True):
        gt = (self._obj.loc[idx[:, :, "forward", "G", "T"]]
              + self._obj.loc[idx[:, :, "reverse", "C", "A"]])
        gt.columns = pd.MultiIndex.from_tuples(product(["GT"], gt.columns))
        ca = (self._obj.loc[idx[:, :, "forward", "C", "A"]]
              + self._obj.loc[idx[:, :, "reverse", "G", "T"]])
        ca.columns = pd.MultiIndex.from_tuples(product(["CA"], ca.columns))
        trans_bias = pd.concat([gt, ca], axis=1)
        # trans_bias = trans_bias.reset_index().set_index(["project_id", "case_id"])
        trans_bias.sort_index(inplace=True)
        # if add_ratio:
        #     trans_bias["ratio"] = trans_bias.GT / trans_bias.CA
        # if add_fraction:
        #     trans_bias["fraction"] = trans_bias.GT / trans_bias.sum(axis=1)
        return trans_bias

    def reference_bias_counts(self, change_numerator, change_denominator):
        ref_bias = self._obj.groupby(["project_id", "case_id", "ref", "alt"]).agg(sum)
        gt = ref_bias.loc[idx[:, :, "G", "T"]]
        gt.columns = pd.MultiIndex.from_tuples(product(["GT"], gt.columns))
        ca = ref_bias.loc[idx[:, :, "C", "A"]]
        ca.columns = pd.MultiIndex.from_tuples(product(["CA"], ca.columns))
        ref_bias = pd.concat([gt, ca], axis=1)
        ref_bias.sort_index(inplace=True)
        # if add_ratio:
        #     ref_bias["ratio"] = ref_bias.GT / ref_bias.CA
        # if add_fraction:
        #     ref_bias["fraction"] = ref_bias.GT / ref_bias.sum(axis=1)
        return ref_bias

    def test_binomial(self, bias_type):
        bias_data = self._obj.ukb_var.calculate_bias(bias_type)
        binoms = dict()
        for col in bias_data.columns.unique(1):
            binoms[col] = stats.binomtest(bias_data.GT[col].sum(),
                                          (bias_data.GT[col] + bias_data.CA[col]).sum())
        return binoms

    def _plot_bias_scatter(self, bias_type, variant_type):
        """Bias type must be one of 'transcription' or 'reference'."""
        plot = go.Figure()
        bias_data = (self._obj.ukb_var.calculate_bias(bias_type)
                     .loc[:, idx[:, variant_type]])
        # bias_data.columns = bias_data.columns.levels[0]
        bias_data = bias_data.droplevel(1, axis=1)
        title = f"{bias_type.title()} strand 8-oxo-G {variant_type} variant call bias"
        for batch in bias_data.index.unique("project_id"):
            batch_subset = bias_data.loc[batch]
            plot.add_trace(go.Scatter(x=batch_subset.CA,
                                      y=batch_subset.GT,
                                      name=batch, mode="markers"))
        max_point = max([bias_data.CA.max(), bias_data.GT.max()])
        min_point = min([bias_data.CA.min(), bias_data.GT.min()])
        plot.add_trace(go.Scatter(x=[min_point, max_point],
                                  y=[min_point, max_point], name="y = x",
                                  mode="lines", line=dict(color="lightgray")))
        plot.update_layout(title=title, xaxis_title="C>A calls",
                           yaxis_title="G>T calls")
        return plot

    def _plot_bias_hexbin(self, bias_type, variant_type):
        bias_data = (self._obj.ukb_var.calculate_bias(bias_type)
                     .loc[:, idx[:, variant_type]])
        bias_data = bias_data.droplevel(1, axis=1)
        hex = plot_hexbin(bias_data.CA, bias_data.GT)
        max_point = max([bias_data.CA.max(), bias_data.GT.max()])
        min_point = min([bias_data.CA.min(), bias_data.GT.min()])
        hex.add_trace(go.Scatter(x=[min_point, max_point],
                                 y=[min_point, max_point], name="y = x",
                                 mode="lines", line=dict(color="lightgray")))
        return hex
    # def plot_reference_bias_scatter(self):


# %% Compare mismatches and calls.

def plot_mismatch_vs_call_hexbin(mismatch_data, call_data, bias_type, variant_type):
    mismatch_bias = mismatch_data.ukb.calculate_bias(bias_type, add_fraction=False)
    call_bias = call_data.ukb_var.calculate_bias(bias_type)
    call_bias = call_bias.GT / call_bias.CA
    hex = plot_hexbin(mismatch_bias.ratio, call_bias[variant_type])
    return hex


def test_mismatch_vs_call_correlation(mismatch_data, call_data, bias_type, variant_type,
                                      method="spearman"):
    mismatch_bias = mismatch_data.ukb.calculate_bias(bias_type, add_fraction=False)
    call_bias = call_data.ukb_var.calculate_bias(bias_type)
    call_bias = call_bias.GT / call_bias.CA
    if method == "spearman":
        stattest = stats.spearmanr
    elif method == "pearson":
        stattest = stats.pearsonr
    else:
        raise ValueError(f"Method must be 'pearson' or 'spearman', got '{method}'.")
    stat_test = stattest(mismatch_bias.ratio, call_bias[variant_type])
    return stat_test
