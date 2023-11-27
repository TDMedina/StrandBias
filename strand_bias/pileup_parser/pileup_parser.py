
import argparse
from collections import Counter
import re

import numpy as np
import pandas as pd
from pandas import CategoricalDtype, MultiIndex
from pandas import IndexSlice as idx
from pandas.api.extensions import register_dataframe_accessor

from pileup_regex import PILEUP_POS_RE, FORWARD_BASES, REVERSE_BASES, MATCH_SET

_BASE_ORDER = ["A", "a", "C", "c", "G", "g", "T", "t", ".", ","]

_ORIENTATIONS = {"F1R2", "F2R1", "F1", "R2", "F2", "R1"}
_INCLUSION_MAP = {"all": _ORIENTATIONS,
                  "combined": {"F1R2", "F2R1"},
                  "non_combined": _ORIENTATIONS - {"F1R2", "F2R1"}}


class PileupBase:
    """Object representing a single pileup character.

    Not to be confused with a PileupPosition, which is a collection of all PileupBase's
    at a single locus in a pileup file.
    """
    def __init__(self, start_code="", base_code="", indel_code="", end_code="",
                 base_qual="", read_pos="", pileup_string=None):
        if pileup_string is not None:
            base = re.findall(PILEUP_POS_RE, pileup_string)[0]
            start_code, base_code, indel_code, end_code = base
        self.is_start, self.mapq = self._parse_start_code(start_code)
        self.base_code = base_code
        self.indel_code = indel_code
        self.is_end = True if end_code == "$" else False

        self.base_qual = base_qual
        self.read_pos = read_pos

    def __repr__(self):
        string = (f"PileupBase("
                  f"start_code='{self.start_code}', "
                  f"base_code='{self.base_code}', "
                  f"indel_code='{self.indel_code}', "
                  f"end_code='{'$' if self.is_end else ''}', "
                  f"base_qual='{self.base_qual}', "
                  f"read_pos={self.read_pos}"
                  f")")
        return string

    def __str__(self):
        string = f"PileupBase('{self.build_pileup_string()}')"
        return string

    @staticmethod
    def _parse_start_code(start_code):
        if start_code:
            return True, start_code[1]
        return False, None

    @property
    def start_code(self):
        string = "^" + self.mapq if self.is_start else ""
        return string

    def build_pileup_string(self):
        string = self.start_code + self.base_code + self.indel_code
        if self.is_end:
            string += "$"
        return string

    @property
    def is_forward(self):
        return self.base_code in FORWARD_BASES

    @property
    def is_reverse(self):
        return self.base_code in REVERSE_BASES

    @property
    def is_match(self):
        return self.base_code in MATCH_SET


class PileupPosition:
    """Object representing a single locus in a pileup file.

    Intended to contain all PileupBase's at a single locus.
    """
    def __init__(self, contig: str, position: int, reference: str, depth: int,
                 pileup_string: str, qual_string: str, alignment_pos_string: str = None,
                 skip_match_bases=False):
        self.contig = contig
        self.position = int(position)
        self.reference = reference
        self.depth = int(depth)
        self.match_filtered = skip_match_bases

        self.pileup_bases, self.pileup_counts = self._parse_pileup_bases(
            pileup_string=pileup_string,
            alignment_pos_string=alignment_pos_string,
            qual_string=qual_string,
            skip_match_bases=skip_match_bases
            )

    def __repr__(self):
        string = (f"PileupPosition("
                  f"contig='{self.contig}', "
                  f"position={self.position}, "
                  f"reference='{self.reference}', "
                  f"depth={self.depth}, "
                  f"pileup_string='{self.build_pileup_string()}', "
                  f"qual_string='{self.build_qual_string()}', "
                  f"alignment_pos_string='{self.build_alignment_pos_string()}'"
                  f")")
        return string

    def __str__(self):
        string = ", ".join(f"'{x}': {y}" for x, y in self.pileup_counts.items())
        string = (f"PileupPosition("
                  f"{self.contig}:{self.position}:{string}, "
                  f"depth={self.depth}"
                  f")")
        return string

    @staticmethod
    def _parse_alignment_pos_string(alignment_pos_string):
        alignment_pos = [int(x) for x in alignment_pos_string.split(",")]
        return alignment_pos

    def _parse_pileup_bases(self, pileup_string, alignment_pos_string, qual_string,
                            skip_match_bases=False):
        positions = self._parse_alignment_pos_string(alignment_pos_string)
        other_data = zip(qual_string, positions)
        bases = re.findall(PILEUP_POS_RE, pileup_string)
        bases = [base + other for base, other in zip(bases, other_data)]
        bases = [PileupBase(*base) for base in bases]
        if skip_match_bases:
            bases = [base for base in bases if not base.is_match]
        counts = Counter(base.base_code for base in bases)
        return bases, counts

    @property
    def is_only_matches(self):
        matches = (self.pileup_counts.keys() == MATCH_SET
                   or self.pileup_counts == Counter())
        return matches

    @property
    def locus(self):
        return self.contig, self.position

    def build_pileup_string(self):
        string = "".join(base.build_pileup_string() for base in self.pileup_bases)
        return string

    def build_qual_string(self):
        string = "".join(base.base_qual for base in self.pileup_bases)
        return string

    def build_alignment_pos_string(self):
        string = ",".join(str(base.read_pos) for base in self.pileup_bases)
        return string

    def count_strands(self):
        strand = ["reverse", "forward"]
        count = Counter(strand[base.is_forward] for base in self.pileup_bases)
        return count

    def count_mismatches(self):
        mismatches = sum(self.pileup_counts[x] for x in _BASE_ORDER[:-2])
        return mismatches


class Pileup:
    """Object containing the contents of one whole pileup file.

    Intended to contain all PileupPosition's from the file.
    """
    def __init__(self, pileup_positions: list[PileupPosition] = None, pileup_file=None,
                 skip_match_positions=False, skip_match_bases=False):
        self.pileup_positions = pileup_positions
        if pileup_file is not None:
            self.pileup_positions = self._read_pileup_file(pileup_file,
                                                           skip_match_positions,
                                                           skip_match_bases)
        self.size = len(self.pileup_positions)
        self.pileup_counts = None
        self.pileup_counts_against_ref = None

    @staticmethod
    def _read_pileup_file(pileup_file, skip_match_positions=False, skip_match_bases=False):
        """Read and parse a pileup file.

        This is a static method that returns a list of PileupPosition's. Using
        Pileup(pileup_file='filepath') is more useful, as it uses this method to
        create a Pileup object, which provides additional tools.
        """
        pileups = []
        with open(pileup_file) as infile:
            for line in infile:
                line = line.strip().split("\t")
                pileup_pos = PileupPosition(*line, skip_match_bases=skip_match_bases)
                if skip_match_positions and pileup_pos.is_only_matches:
                    continue
                pileups.append(pileup_pos)
        return pileups

    def count_all_pileup_bases(self):
        """Tally the total number of PileupBase's contained.

        This method counts how many of each PileupBase symbol (ACGT.acgt,) are
        contained within all contained PileupPosition's. The result is stored in the
        pileup_counts attribute, so unless recounting after filtering, access the
        attribute instead after running this method once.
        """
        # if self.pileup_counts is not None:
        #     return self.pileup_counts
        counter = Counter()
        for pileup_pos in self.pileup_positions:
            counter += pileup_pos.pileup_counts
        self.pileup_counts = counter
        return counter

    def count_all_pileup_bases_against_reference(self):
        if self.pileup_counts_against_ref is not None:
            return self.pileup_counts_against_ref
        ref_counter = {"A": Counter(), "C": Counter(), "G": Counter(), "T": Counter()}
        for pileup_pos in self.pileup_positions:
            ref_counter[pileup_pos.reference] += pileup_pos.pileup_counts
        self.pileup_counts_against_ref = ref_counter
        return ref_counter

    def tabulate_pileup_balance_against_ref(self):
        if self.pileup_counts_against_ref is None:
            self.count_all_pileup_bases_against_reference()
        counts = self.pileup_counts_against_ref
        table = {"A": dict(), "C": dict(), "G": dict(), "T": dict()}
        for base1 in table:
            for base2 in table:
                if base1 == base2:
                    table[base1][base2] = None
                    continue
                try:
                    table[base1][base2] = (counts[base1][base2]
                                           / counts[base1][base2.lower()])
                except ZeroDivisionError:
                    table[base1][base2] = np.inf
        table = pd.DataFrame.from_dict(table, orient="index")
        return table

    def tabulate_pileup_counts_against_ref(self):
        if self.pileup_counts_against_ref is None:
            self.count_all_pileup_bases_against_reference()
        table = pd.DataFrame.from_dict(self.pileup_counts_against_ref, orient="index")
        sorted_cols = sorted(table.columns, key=lambda x: _BASE_ORDER.index(x))
        table = table.loc[:, sorted_cols]
        return table

    def filter_pileup_positions(self, filter_dictionary):
        pileup_positions = [pos for pos in self.pileup_positions
                            if pos.locus[1] not in filter_dictionary[pos.locus[0]]]
        pileup = Pileup(pileup_positions)
        return pileup

    def filter_single_mismatches(self):
        singles = [pileup_pos for pileup_pos in self.pileup_positions
                   if pileup_pos.count_mismatches() == 1]
        pileup = Pileup(singles)
        return pileup


class CodingRegionPileups:
    def __init__(self, region_name,
                 F1: Pileup = None, R2: Pileup = None, F1R2: Pileup = None,
                 F2: Pileup = None, R1: Pileup = None, F2R1: Pileup = None):
        self.region_name = region_name
        self.F1, self.R2, self.F1R2 = F1, R2, F1R2
        self.F2, self.R1, self.F2R1 = F2, R1, F2R1

    @property
    def orientations(self) -> [Pileup]:
        orientations = [value for ori in _ORIENTATIONS
                        if (value := self.__getattribute__(ori)) is not None]
        return orientations


class SamplePileups:
    def __init__(self,
                 forward_coding: CodingRegionPileups,
                 reverse_coding: CodingRegionPileups):
        self.forward_coding = forward_coding
        self.reverse_coding = reverse_coding

    @staticmethod
    def read_pileups(file_prefix, include="all"):
        pileups = dict()
        orientations = _INCLUSION_MAP[include]
        for strand in {"forward_coding", "reverse_coding"}:
            pileups[strand] = dict()
            for orientation in orientations:
                file = f"{file_prefix}.filtered.{strand}.{orientation}.no_match_positions.pileup"
                pileups[strand][orientation] = Pileup(pileup_file=file)
                pileups[strand][orientation].count_all_pileup_bases_against_reference()
            pileups[strand] = CodingRegionPileups(strand, **pileups[strand])
        pileups = SamplePileups(**pileups)
        return pileups

    @property
    def regions(self):
        return [self.forward_coding, self.reverse_coding]

    def filter_single_mismatches(self):
        for region in self.regions:
            for pileup in region.orientations:
                pileup.filter_single_mismatches()

    def filter_positions(self, filter_dict):
        for region in self.regions:
            for pileup in region.orientations:
                pileup.filter_pileup_positions(filter_dict)

    def tabulate(self):
        tables = []
        col_dex = MultiIndex.from_product([["A", "C", "G", "T", "match"],
                                           ["forward", "reverse"]],
                                          names=["alt", "alignment"])
        for coding_region in [self.forward_coding, self.reverse_coding]:
            for orientation in _ORIENTATIONS:
                if coding_region.__getattribute__(orientation) is None:
                    continue
                pileup = coding_region.__getattribute__(orientation)
                table = pileup.tabulate_pileup_counts_against_ref()
                table = table.reset_index()
                table = table.rename(columns={"index": "reference"})
                table["coding_strand"] = coding_region.region_name.split("_")[0]
                table["orientation"] = orientation

                cat_cols = ["coding_strand", "orientation", "reference"]
                dtypes = [["forward", "reverse"], _ORIENTATIONS, list("ACGT")]
                for col, dtype in zip(cat_cols, dtypes):
                    table[col] = table[col].astype(CategoricalDtype(dtype))
                table = table.set_index(["reference", "coding_strand", "orientation"])

                tables.append(table)

        table = pd.concat(tables)
        table = table.sort_index()
        table = table.loc[:, sorted(table.columns, key=lambda x: _BASE_ORDER.index(x))]
        table.columns = col_dex
        return table


@register_dataframe_accessor("pileup_tools")
class PileupTable:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(file_path):
        table = pd.read_csv(file_path, sep="\t", index_col=[0, 1, 2], header=[0, 1])
        return table

    def make_asymmetry_summary_table(self, by="coding_strand",
                                     as_proportion=False, as_ratio=False):
        summary_table = (self._obj
                         .groupby(["reference", by]).agg(sum)
                         .groupby(axis=1, level="alt").agg(sum))
        if as_proportion:
            summary_table = summary_table.groupby("reference").agg(self._make_summary_proportion)
        elif as_ratio:
            summary_table = summary_table.groupby("reference").agg(self._make_summary_ratio)
        summary_table = summary_table.loc[["C", "G"], ["A", "T"]]
        return summary_table

    def calculate_orientation_bias_by_coding_strand(self):
        results = self._obj.groupby("alt", axis=1).agg(sum)
        results = results.loc[idx[:, :, "F1R2"]] / results.loc[idx[:, :, "F2R1"]]
        return results

    def calculate_strand_bias_by_coding_strand(self):
        results = self._obj.groupby(["reference", "coding_strand"]).agg(sum)
        results = results.groupby("alt", axis=1).agg(self._agg_div_alignment)
        return results

    @staticmethod
    def _agg_div_alignment(df):
        df = df.droplevel("alt", axis=1)
        return df.forward / df.reverse

    @staticmethod
    def _make_summary_proportion(pair):
        if not len(pair) > 1 or min(pair) == 0:
            return pd.NA
        norm = min(pair)
        pair = ["1" if val == norm else f"{val/norm:.2f}" for val in pair]
        pair = f"{pair[0]}:{pair[1]}"
        return pair

    @staticmethod
    def _make_summary_ratio(pair):
        if not len(pair) > 1 or min(pair) == 0:
            return pd.NA
        pair = pair[0] / pair[1]
        return pair

# def read_pileups_split_by_orientation(file_prefix_id, include="all"):
#     pileups = dict()
#     orientations = {"F1R2", "F2R1", "F1", "R2", "F2", "R1"}
#     if include == "all":
#         pass
#     elif include == "combined":
#         orientations = {"F1R2", "F2R1"}
#     elif include == "non_combined":
#         orientations -= {"F1R2", "F2R1"}
#     for strand in {"forward", "reverse"}:
#         pileups[strand] = dict()
#         for orientation in orientations:
#             file = f"{file_prefix_id}.filtered.{strand}_coding.{orientation}.no_match_positions.pileup"
#             pileups[strand][orientation] = Pileup(pileup_file=file)
#             pileups[strand][orientation].count_all_pileup_bases_against_reference()
#     return pileups


# def tabulate_pileups_split_by_orientation(pileups_dict):
#     tables = []
#     orientations = ["F1", "R2", "F1R2", "F2", "R1", "F2R1"]
#     col_dex = MultiIndex.from_product([["A", "C", "G", "T", "match"],
#                                       ["forward", "reverse"]],
#                                       names=["alt", "alignment"])
#
#     for strand, strand_dict in pileups_dict.items():
#         for orientation, pileup in strand_dict.items():
#             table = pileup.tabulate_pileup_counts_against_ref()
#
#             table = table.reset_index()
#             table = table.rename(columns={"index": "reference"})
#             table["coding_strand"] = strand
#             table["orientation"] = orientation
#
#             cat_cols = ["coding_strand", "orientation", "reference"]
#             dtypes = [["forward", "reverse"], orientations, list("ACGT")]
#             for col, dtype in zip(cat_cols, dtypes):
#                 table[col] = table[col].astype(CategoricalDtype(dtype))
#             table = table.set_index(["reference", "coding_strand", "orientation"])
#
#             tables.append(table)
#     tables = pd.concat(tables)
#     tables = tables.sort_index()
#     tables = tables.loc[:, sorted(tables.columns, key=lambda x: _BASE_ORDER.index(x))]
#     # tables = tables.loc[:, ["A", "a", "C", "c", "G", "g", "T", "t", ".", ","]]
#     tables.columns = col_dex
#     return tables


def read_filter_bed(filter_file):
    """Read in a BED file and create large sets of loci.

    Warning: Very memory intensive.
    """
    filter_dict = dict()
    with open(filter_file) as infile:
        for line in infile:
            contig, start, stop = line.strip().split("\t")
            if contig not in filter_dict:
                filter_dict[contig] = set()
            filter_range = range(int(start), int(stop))
            filter_dict[contig].update(set(filter_range))
    return filter_dict


# def filter_single_mismatches_in_pileup_dict(pileup_dict):
#     pileups = {
#         strand: {ori: pileup.filter_single_mismatches()
#                  for ori, pileup in ori_dict.items()}
#         for strand, ori_dict in pileup_dict.items()
#         }
#     return pileups


# def filter_pileup_positions_in_pileup_dict(pileup_dict, filter_dict):
#     pileups = {
#         strand: {ori: pileup.filter_pileup_positions(filter_dict)
#                  for ori, pileup in ori_dict.items()}
#         for strand, ori_dict in pileup_dict.items()
#         }
#     return pileups


# def make_asymmetry_summary_table(table, by="coding_strand", as_proportion=False, as_ratio=False):
#     summary_table = (
#         table
#         .groupby(["reference", by]).agg(sum)
#         .groupby(axis=1, level="alt").agg(sum)
#         )
#     if as_proportion:
#         summary_table = summary_table.groupby("reference").agg(_make_summary_proportion)
#     elif as_ratio:
#         summary_table = summary_table.groupby("reference").agg(_make_summary_ratio)
#     summary_table = summary_table.loc[["C", "G"], ["A", "T"]]
#     return summary_table
#
#
# def calculate_orientation_bias_by_coding_strand(table):
#     results = table.groupby("alt", axis=1).agg(sum)
#     results = results.loc[idx[:, :, "F1R2"]] / results.loc[idx[:, :, "F2R1"]]
#     return results
#
#
# def calculate_strand_bias_by_coding_strand(table):
#     results = table.groupby(["reference", "coding_strand"]).agg(sum)
#     results = results.groupby("alt", axis=1).agg(_agg_div_alignment)
#     return results
#
#
# def _agg_div_alignment(df):
#     df = df.droplevel("alt", axis=1)
#     return df.forward / df.reverse
#
#
# def _make_summary_proportion(pair):
#     if not len(pair) > 1 or min(pair) == 0:
#         return pd.NA
#     norm = min(pair)
#     pair = ["1" if val == norm else f"{val/norm:.2f}" for val in pair]
#     pair = f"{pair[0]}:{pair[1]}"
#     return pair
#
#
# def _make_summary_ratio(pair):
#     if not len(pair) > 1 or min(pair) == 0:
#         return pd.NA
#     pair = pair[0] / pair[1]
#     return pair


# def read_asymmetry_table(table_path):
#     table = pd.read_csv(table_path, sep="\t", index_col=[0, 1, 2], header=[0, 1])
#     return table


def main(file_prefix, include="all",
         single_mismatches_only=False, filter_bed=None,
         export_path=None, summary_path=None):
    sample = SamplePileups.read_pileups(file_prefix, include)
    if single_mismatches_only:
        sample.filter_single_mismatches()
    if filter_bed is not None:
        filter_dict = read_filter_bed(filter_bed)
        sample.filter_positions(filter_dict)

    sample_table = sample.tabulate()

    if export_path is not None:
        sample_table.to_csv(export_path, sep="\t", index=True)
    if summary_path is not None:
        summary_table = sample_table.pileup_tools.make_asymmetry_summary_table()
        summary_table.to_csv(summary_path, sep="\t", index=True)
    return sample


# def main(file_prefix, filter_bed=None, export_path=None, summary_path=None, single_mismatches_only=False,
#          include="all"):
#     pileups = read_pileups_split_by_orientation(file_prefix, include)
#     if filter_bed is not None:
#         filter_dict = read_filter_bed(filter_bed)
#         pileups = filter_pileup_positions_in_pileup_dict(pileups, filter_dict)
#     if single_mismatches_only:
#         pileups = filter_single_mismatches_in_pileup_dict(pileups)
#     table = tabulate_pileups_split_by_orientation(pileups)
#     if export_path is not None:
#         table.to_csv(export_path, sep="\t", index=True)
#     if summary_path is not None:
#         summary = make_asymmetry_summary_table(table)
#         summary.to_csv(summary_path, sep="\t", index=True)
#     return table


def _setup_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file-prefix", required=True)
    parser.add_argument("-b", "--filter-bed")
    parser.add_argument("-o", "--export-path")
    parser.add_argument("--summary-path")
    parser.add_argument("-s", "--singles-only", dest="single_mismatches_only",
                        action="store_true")
    parser.add_argument("-m", "--include-mode", dest="include",
                        choices=["all", "combined", "non_combined"])
    return parser


if __name__ == '__main__':
    import sys
    pileup_sample = main(**vars(_setup_argparse().parse_args(sys.argv[1:])))
