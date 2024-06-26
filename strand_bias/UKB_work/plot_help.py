
from collections import namedtuple

VarCounts = namedtuple("VarCounts", ["filtered", "unfiltered"])


class PlotDirectory:
    def __init__(self, mismatch_table, variant_tables: VarCounts):
        self.mismatches = PlotsByBiasType(mismatch_table, "mismatches")
        self.variants = PlotsByFilterStatus(variant_tables)

    def show_all(self):
        for plots in [self.mismatches, self.variants]:
            plots.show_all()


class PlotsByFilterStatus:
    def __init__(self, variant_table):
        self.filtered = PlotsByBiasType(variant_table.filtered, "variants", filtered=True)
        self.unfiltered = PlotsByBiasType(variant_table.unfiltered, "variants", filtered=False)

    def show_all(self):
        for plots in [self.filtered, self.unfiltered]:
            plots.show_all()


class PlotsByBiasType:
    def __init__(self, data_table, data_type, inverse=False, filtered=False):
        params = dict(data_table=data_table, data_type=data_type, inverse=inverse, filtered=filtered)
        self.reference_bias = PlotsByPlotType(bias_type="reference", **params)
        self.transcription_bias = PlotsByPlotType(bias_type="transcription", **params)

    def show_all(self):
        for plots in [self.reference_bias, self.transcription_bias]:
            plots.show_all()


class PlotsByPlotType:
    def __init__(self, data_table, data_type, bias_type, inverse=False, filtered=False):
        if data_type not in {"variants", "mismatches"}:
            raise ValueError("'data_type' must be 'variants'"
                             f" or 'mismatches', not '{data_type}'.")
        params = dict(bias_type=bias_type, inverse=inverse, filtered=filtered)
        attr_name = f"ukb_{data_type}"
        self.box = data_table.__getattr__(attr_name)._plot_bias_box(**params)
        self.scatter = data_table.__getattr__(attr_name)._plot_bias_scatter(**params)
        self.hexbin = data_table.__getattr__(attr_name)._plot_bias_hexbin(**params)

    def show_all(self):
        self.box.show()
        self.scatter.show()
        if hasattr(self, "hexbin"):
            self.hexbin.show()
