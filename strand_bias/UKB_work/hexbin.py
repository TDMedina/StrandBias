
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
