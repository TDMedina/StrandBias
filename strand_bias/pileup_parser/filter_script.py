
import sys

import pileup_parser

with open(sys.argv[1]) as infile:
    prefix_list = infile.readlines()
prefix_list = [prefix.strip() for prefix in prefix_list]

filter_dict = pileup_parser.read_filter_bed(sys.argv[2])

for prefix in prefix_list:
    pileups = pileup_parser.read_pileups_split_by_orientation(prefix, include="combined")
    pileups = pileup_parser.filter_single_mismatches_in_pileup_dict(pileups)
    pileups = pileup_parser.filter_pileup_positions_in_pileup_dict(pileups, filter_dict)
    table = pileup_parser.tabulate_pileups_split_by_orientation(pileups)
    table.to_csv(f"{prefix}.asym_table.filtered.tsv", sep="\t", index=True)
