
import argparse
from collections import Counter
import re

import numpy as np
import pandas as pd
from pandas import CategoricalDtype

from pileup_regex import PILEUP_POS_RE, FORWARD_BASES, REVERSE_BASES, MATCH_SET


class Pileup:
    def __init__(self, pileup_positions=None, pileup_file=None,
                 skip_match_positions=False, skip_match_bases=False):
        self.pileup_positions = pileup_positions
        if pileup_file is not None:
            self.pileup_positions = self.read_pileup_file(pileup_file,
                                                          skip_match_positions,
                                                          skip_match_bases)
        self.size = len(self.pileup_positions)
        self.pileup_counts = None
        self.pileup_counts_against_ref = None

    @staticmethod
    def read_pileup_file(pileup_file, skip_match_positions=False, skip_match_bases=False):
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
        if self.pileup_counts is not None:
            return self.pileup_counts
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
        table = table.loc[:, ["A", "a", "C", "c", "G", "g", "T", "t", ".", ","]]
        return table

    def filter_pileup_positions(self, filter_dictionary):
        pileup_positions = [pos for pos in self.pileup_positions
                            if pos.locus[1] not in filter_dictionary[pos.locus[0]]]
        pileup = Pileup(pileup_positions)
        return pileup


class PileupPosition:
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


class PileupBase:
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


def read_pileups_split_by_orientation(file_prefix_id):
    pileups = dict()
    for strand in {"forward", "reverse"}:
        pileups[strand] = dict()
        for orientation in {"F1R2", "F2R1", "F1", "R2", "F2", "R1"}:
            file = f"{file_prefix_id}.filtered.{strand}_coding.{orientation}.no_match_positions.pileup"
            pileups[strand][orientation] = Pileup(pileup_file=file)
            pileups[strand][orientation].count_all_pileup_bases_against_reference()
    return pileups


def tabulate_pileups_split_by_orientation(pileups_dict):
    tables = []
    orientations = ["F1", "R2", "F1R2", "F2", "R1", "F2R1"]
    for strand in {"forward", "reverse"}:
        for orientation in orientations:
            table = pileups_dict[strand][orientation].tabulate_pileup_counts_against_ref()

            table = table.reset_index()
            table = table.rename(columns={"index": "reference"})
            table["coding_strand"] = [strand]*4
            table["orientation"] = [orientation]*4

            cat_cols = ["coding_strand", "orientation", "reference"]
            dtypes = [["forward", "reverse"], orientations, list("ACGT")]
            for col, dtype in zip(cat_cols, dtypes):
                table[col] = table[col].astype(CategoricalDtype(dtype))
            table = table.set_index(["reference", "coding_strand", "orientation"])

            tables.append(table)
    tables = pd.concat(tables)
    tables = tables.sort_index()
    tables = tables.loc[:, ["A", "a", "C", "c", "G", "g", "T", "t", ".", ","]]
    return tables


# def calculate_gtca_asymmetry(pileup_table):


def read_filter_bed(filter_file):
    '''Read in a BED file and create large sets of loci.

    Warning: Very memory intensive.
    '''
    filter_dict = dict()
    with open(filter_file) as infile:
        for line in infile:
            contig, start, stop = line.strip().split("\t")
            if contig not in filter_dict:
                filter_dict[contig] = set()
            filter_range = range(int(start), int(stop))
            filter_dict[contig].update(set(filter_range))
    return filter_dict


def main(file_prefix, filter_bed=None, export_path=None):
    pileups = read_pileups_split_by_orientation(file_prefix)
    if filter_bed is not None:
        filter_dict = read_filter_bed(filter_bed)
        pileups = {
            strand: {ori: pileup.filter_pileup_positions(filter_dict)
                     for ori, pileup in ori_dict.items()}
            for strand, ori_dict in pileups.items()
            }
    table = tabulate_pileups_split_by_orientation(pileups)
    if export_path is not None:
        table.to_csv(export_path, sep="\t", index=True)
    return table


def _setup_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file-prefix", required=True)
    parser.add_argument("-b", "--filter-bed")
    parser.add_argument("-o", "--export-path")
    return parser


# if __name__ == '__main__':
#     import sys
#     table = main(**vars(_setup_argparse().parse_args(sys.argv[1:])))
#     print()
