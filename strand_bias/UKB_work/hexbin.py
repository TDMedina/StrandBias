
from statistics import mode

from matplotlib.pyplot import hexbin
from matplotlib import colormaps
import numpy as np
import pandas as pd
import plotly.graph_objects as go


# %% Hexbinning.

def make_hexbin_data(x, y):
    hexbin_data = hexbin(x, y)
    hexbin_data = np.array(np.hstack([hexbin_data.get_offsets(),
                                      hexbin_data.get_array().reshape([-1, 1])]))
    hexbin_data = pd.DataFrame(hexbin_data, columns=["x", "y", "counts"])
    hexbin_data = hexbin_data.loc[hexbin_data.counts > 0]
    hexbin_data["log2count"] = [np.log2(x) for x in hexbin_data.counts]
    return hexbin_data


def make_hexagons(hexbin_data, colormap="viridis", dims=(1200, 1200),
                  showlegend=True, return_graph_objects=False):
    width, height = dims
    gos = []
    dx = mode(np.diff(sorted(set(hexbin_data.x))))
    dy = mode(np.diff(sorted(set(hexbin_data.y))))

    fig = go.Figure()
    viridis = colormaps[colormap]
    max_z = np.log2(hexbin_data.counts).max()
    for x, y, z in zip(hexbin_data.x, hexbin_data.y, hexbin_data.counts):
        color = ",".join(str(float(val)) for val in viridis(np.log2(z)/max_z))
        color = f"rgba({color})"
        hex_points = [(x, y+2*dy/3), (x+dx, y+dy/3), (x+dx, y-dy/3), (x, y-2*dy/3),
                      (x-dx, y-dy/3), (x-dx, y+dy/3), (x, y+2*dy/3)]
        hex_points = list(zip(*hex_points))
        hexagons = go.Scatter(x=hex_points[0], y=hex_points[1], mode="none",
                              fill="toself", fillcolor=color, showlegend=False)
        if return_graph_objects:
            gos.append(hexagons)
            continue
        fig.add_trace(hexagons)
    if not showlegend:
        if return_graph_objects:
            return gos
        return fig
    legend_thing = go.Scatter(x=[None], y=[None], mode="markers", showlegend=False,
                             marker=dict(colorbar=dict(title="Sample Count"),
                                         colorscale=colormap,
                                         showscale=True,
                                         cmin=1, cmax=max(hexbin_data.counts)))
    if return_graph_objects:
        gos.append(legend_thing)
        return gos
    fig.add_trace(legend_thing)
    return fig


def plot_hexbin(x, y, colormap="viridis", dims=(1200, 1200),
                return_graph_objects=False, showlegend=True):
    hexbin_data = make_hexbin_data(x, y)
    fig = make_hexagons(hexbin_data, colormap, dims, showlegend=showlegend,
                        return_graph_objects=return_graph_objects)
    return fig
