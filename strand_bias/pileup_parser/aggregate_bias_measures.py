
import argparse
from pathlib import Path

from pileup_parser import PileupTable


def main(path_file):
    with open(path_file) as infile:
        paths = infile.readlines()
    paths = [path.rstrip() for path in paths]

    for path in paths:
        path = Path(path)
        dest_dir = path.parent
        table = PileupTable.read_csv(path)

        bias = table.pileup_tools.calculate_orientation_bias_by_coding_strand()
        bias.to_csv(dest_dir/path.stem+".ori_bias.tsv", sep="\t")

        bias = table.pileup_tools.calculate_strand_bias_by_coding_strand()
        bias.to_csv(dest_dir/path.stem+".strand_bias.tsv", sep="\t")


def _setup_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--path-file")
    return parser


if __name__ == '__main__':
    args = _setup_argparser().parse_args()
    main(args.path_file)
