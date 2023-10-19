
from statistics import mode

from matplotlib.pyplot import hexbin
from matplotlib import colormaps
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def make_hexbin_data(x, y):
    hexbin_data = hexbin(x, y)
    hexbin_data = np.array(np.hstack([hexbin_data.get_offsets(),
                                      hexbin_data.get_array().reshape([-1, 1])]))
    hexbin_data = pd.DataFrame(hexbin_data, columns=["x", "y", "counts"])
    hexbin_data = hexbin_data.loc[hexbin_data.counts > 0]
    hexbin_data["log2count"] = [np.log2(x) for x in hexbin_data.counts]
    return hexbin_data


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
                                 fill="toself", fillcolor=color, showlegend=False,
                                 text=f"C➔A Ratio: {round(x, 3)}<br>"
                                      f"SNVs: {round(y, 3)}<br>"
                                      f"Count: {int(z)}"))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", showlegend=False,
                             marker=dict(colorbar=dict(title="Sample Count"),
                                         colorscale=colormap,
                                         showscale=True,
                                         cmin=1, cmax=max(hexbin_data.counts))))
    fig.update_layout(xaxis_title="C➔A coding/C➔A template",
                      yaxis_title="Singleton SNVs",
                      width=width, height=height)
    return fig


def plot_hexbin(x, y, colormap="viridis", dims=(1200, 1200)):
    hexbin_data = make_hexbin_data(x, y)
    fig = plot_hexagons(hexbin_data, colormap, dims)
    return fig


def plot_asymmetry_data(asym_data, x, y):
    fig = make_subplots(2, 2, shared_xaxes=True, shared_yaxes=True)
    asym_scatter = go.Scatter(
        x=asym_data[x],
        y=asym_data[y],
        mode="markers",
        marker=dict(opacity=0.5),
        customdata=asym_data.index,
        hovertemplate="ID: %{customdata}<br>SNPs: %{y}<br>Ratio: %{x}<extra></extra>")

    y_hist = go.Histogram(y=asym_data[y])
    y_box = go.Box(y=asym_data[y])
    x_hist = go.Histogram(x=asym_data[x])
    x_box = go.Box(x=asym_data[x])

    fig.add_trace(asym_scatter, row=2, col=1)
    fig.add_trace(y_box, row=2, col=1)
    fig.add_trace(x_box, row=2, col=1)
    fig.add_trace(y_hist, row=2, col=2)
    fig.add_trace(x_hist, row=1, col=1)
    fig.update_layout(xaxis_title="C➔A coding/C➔A template",
                      yaxis_title="Singleton heterozygous SNPs")
    return fig


def plot_per_variant_ratios(variant_counts: pd.DataFrame):
    grouped = variant_counts.groupby(["coding_strand", "variant"]).agg(sum)
    ratios = ((grouped.loc["forward"] + grouped.loc["reverse_complement"])
              / (grouped.loc["reverse"] + grouped.loc["forward_complement"]))
    figs = []
    for feature in ["nHets", "nNonRefHom", "nNonRefAlleles"]:
        ratio_fig = go.Figure([go.Bar(x=ratios.index, y=ratios[feature],
                                      name="(F+RC)/(R+FC) ratio"),
                               go.Bar(x=ratios.index, y=ratios[feature] - 1,
                                      name="(F+RC)/(R+FC) ratio - 1"),
                               go.Bar(x=ratios.index, y=abs(ratios[feature] - 1),
                                      name="abs((F+RC)/(R+FC) ratio - 1)")])
        ratio_fig.update_layout(title=f"Ratio of {feature} genotypes in coding vs. template strands")
        figs.append(ratio_fig)
    return figs


def plot_by_batch_with_trendline(asym_data: pd.DataFrame, x, y):
    fig = px.scatter(asym_data,
                     x=x,
                     y=y,
                     color="batch",
                     trendline="rolling",
                     trendline_scope="overall",
                     trendline_options=dict(window=1000),
                     trendline_color_override="red"
                     )
    return fig
