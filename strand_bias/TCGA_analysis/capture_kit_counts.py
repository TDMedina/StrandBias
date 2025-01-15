
import pandas as pd
from pandas.api.extensions import register_dataframe_accessor

from strand_bias.TCGA_analysis.utilities import _COMP_DICT


@register_dataframe_accessor("capkit")
class CaptureKit:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_capture_kit_nucleotide_summary(file_path):
        table = pd.read_csv(file_path, sep="\t", index_col=0)
        return table

    def calculate_reference_counts(self):
        counts = self._obj.sum(axis=0)
        return counts

    def calculate_transcription_counts(self):
        counts = {x: self._obj.loc["forward_coding", x] + self._obj.loc["reverse_coding", y]
                  for x, y in _COMP_DICT.items()}
        counts["total"] = sum(counts.values())
        counts = pd.Series(counts)
        return counts

    def calculate_counts(self, count_type):
        if count_type == "reference":
            return self.calculate_reference_counts()
        elif count_type == "transcription":
            return self.calculate_transcription_counts()
        else:
            raise ValueError(f"'count_type' must be one of 'reference' or"
                             f" 'transcription', not '{count_type}'")
