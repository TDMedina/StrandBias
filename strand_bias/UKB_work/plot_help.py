
from collections import namedtuple

import plotly.graph_objects as go

VarCounts = namedtuple("VarCounts", ["filtered", "unfiltered"])


class PlotDirectory:
    def __init__(self, mismatch_table, variant_tables: VarCounts, normalization_data=None,
                 mismatch_normalization_factor=1, variant_normalization_factor=1,
                 ratio_by_template=True):
        self.mismatches = PlotsByBiasType(mismatch_table, "mismatches",
                                          normalization_data=normalization_data,
                                          normalization_factor=mismatch_normalization_factor,
                                          ratio_by_template=ratio_by_template)
        self.variants = PlotsByFilterStatus(variant_tables,
                                            normalization_data=normalization_data,
                                            normalization_factor=variant_normalization_factor,
                                            ratio_by_template=ratio_by_template)
        self.correlation = CorrPlotsByFilterStatus(mismatch_table, variant_tables)

    def show_all(self):
        for plots in [self.mismatches, self.variants, self.correlation]:
            plots.show_all()

    def show_corr_plots(self):
        self.correlation.show_all()


class PlotsByFilterStatus:
    def __init__(self, variant_table, normalization_data=None, normalization_factor=1, ratio_by_template=True):
        self.unfiltered = PlotsByBiasType(variant_table.unfiltered, "variants", filtered=False,
                                          normalization_data=normalization_data,
                                          normalization_factor=normalization_factor,
                                          ratio_by_template=ratio_by_template)
        self.filtered = PlotsByBiasType(variant_table.filtered, "variants", filtered=True,
                                        normalization_data=normalization_data,
                                        normalization_factor=normalization_factor,
                                        ratio_by_template=ratio_by_template)

    def show_all(self):
        for plots in [self.filtered, self.unfiltered]:
            plots.show_all()


class PlotsByBiasType:
    def __init__(self, data_table, data_type, inverse=False, filtered=False, normalization_data=None,
                 normalization_factor=1, ratio_by_template=True):
        params = dict(data_table=data_table, data_type=data_type, inverse=inverse, filtered=filtered,
                      normalization_data=normalization_data, normalization_factor=normalization_factor,
                      ratio_by_template=ratio_by_template)
        self.reference_bias = PlotsByPlotType(bias_type="reference", **params)
        self.transcription_bias = PlotsByPlotType(bias_type="transcription", **params)

    def show_all(self):
        for plots in [self.reference_bias, self.transcription_bias]:
            plots.show_all()


class PlotsByPlotType:
    def __init__(self, data_table, data_type, bias_type, inverse=False, filtered=False,
                 normalization_data=None, normalization_factor=1, ratio_by_template=True):
        if data_type not in {"variants", "mismatches"}:
            raise ValueError("'data_type' must be 'variants'"
                             f" or 'mismatches', not '{data_type}'.")
        params = dict(bias_type=bias_type, inverse=inverse, filtered=filtered,
                      normalization_data=normalization_data, normalization_factor=normalization_factor,
                      ratio_by_template=ratio_by_template)
        attr_name = f"ukb_{data_type}"
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
