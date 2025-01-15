
import pandas as pd
from pandas import IndexSlice as idx
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

import strand_bias.pileup_parser.aggregate_asym_tables as asym

pio.renderers.default = "browser"

vcrome_kit = asym.read_capture_kit_nucleotide_summary("/home/tyler/Documents/Projects/StrandBias/VCRome.hg38.nt_counts.tsv")

data = asym.read_concatenated_table("/home/tyler/Documents/Projects/StrandBias/concatenated_asym_tables.tsv")
data.reset_index(inplace=True)
data["file_id"] = [x.replace(".asym_table", "") for x in data.file_id]
data.set_index(["file_id", "reference", "coding_strand", "orientation"], inplace=True)


total_base_counts = pd.read_csv("/home/tyler/Documents/Projects/StrandBias/current_total_counts.tsv",
                                sep="\t", index_col=[0, 1, 2])
total_base_counts = total_base_counts.loc[idx[:, :, ["F1R2", "F2R1"]],]
total_base_counts.reset_index(inplace=True)
total_base_counts["file_id"] = [x.replace("_wxs_gdc_realn", "") for x in total_base_counts.file_id]
total_base_counts.set_index(["file_id", "coding_region", "orientation"], inplace=True)

# 1. Mismatch nucleotide enrichment
simple_aggregation = (data.groupby(["file_id", "reference"]).agg(sum)
                          .groupby("alt", axis=1).agg(sum))[list("ACGT")]

simple_aggregation_norm = simple_aggregation.reset_index().set_index("file_id")
total_base_counts_agg = total_base_counts.groupby("file_id").agg(sum)
for base in "ACGT":
    simple_aggregation_norm[base] = [1e6 * count / float(total_base_counts_agg.loc[file_id])
                                     for file_id, count in simple_aggregation_norm[base].items()]
simple_aggregation_norm = simple_aggregation_norm.reset_index().set_index(["file_id", "reference"])


def plot_4_by_4(table):
    plot = make_subplots(4, 4, shared_xaxes=True, shared_yaxes=True,
                         x_title="Reference", y_title="Mismatch",
                         start_cell="bottom-left")
    bases = list(enumerate("ACGT", start=1))
    for i, x in bases:
        for j, y in bases:
            if i == j:
                plot.add_trace(go.Scatter({}), j, i)
                continue
            data = table.loc[idx[:, x], y].sort_values()
            line_color = "black"
            if {x, y} == {"C", "T"} or {x, y} == {"G", "A"}:
                line_color = "red"
            scatter = go.Scatter(x=list(range(data.shape[0])),
                                 y=data,
                                 line_color=line_color)
            plot.add_trace(scatter, j, i)

    for i, x in bases:
        plot.update_xaxes(title_text=x, row=1, col=i)
        plot.update_yaxes(title_text=x, row=i, col=1)
    return plot


