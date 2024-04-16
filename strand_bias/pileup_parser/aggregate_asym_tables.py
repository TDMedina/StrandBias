
import argparse
from pathlib import Path

import pandas as pd


def read_capture_kit_nucleotide_summary(file_path):
    table = pd.read_csv(file_path, sep="\t", index_col=0)
    return table


def read_individual_table(file_path, add_file_id_col=False):
    table = pd.read_csv(file_path, sep="\t", header=[0, 1], index_col=[0, 1, 2])
    if add_file_id_col:
        file_id = Path(file_path).name.replace("_wxs_gdc_realn.asym_summary_table.tsv", "")
        table["file_id"] = file_id
        table = table.reset_index().set_index(
            ["file_id", "reference", "coding_strand", "orientation"]
            )
    return table


def read_concatenated_table(file_path):
    table = pd.read_csv(file_path, sep="\t", header=[0, 1], index_col=list(range(4)))
    return table


def main(path_file, output=None, existing_concatenated_table=None) -> pd.DataFrame:
    with open(path_file) as infile:
        paths = infile.readlines()
    paths = [path.rstrip() for path in paths]

    if existing_concatenated_table:
        table = read_concatenated_table(existing_concatenated_table)
    else:
        path_0 = paths.pop()
        table = read_individual_table(path_0, add_file_id_col=True)

    for path in paths:
        table2 = read_individual_table(path, add_file_id_col=True)
        table = pd.concat([table, table2])

    if output:
        table.to_csv(output, index=True, sep="\t")
    return table


class CustomHelp(argparse.HelpFormatter):
    """Custom help formatter_class that only displays metavar once."""

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar
        parts = []
        if action.nargs == 0:
            parts.extend(action.option_strings)
        else:
            default = self._get_default_metavar_for_optional(action)
            args_string = self._format_args(action, default)
            for option_string in action.option_strings:
                parts.append(f"{option_string}")
            parts[-1] += f" {args_string} "
        return ", ".join(parts)

    def _format_action(self, action):
        parts = super()._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts


def _setup_argparser():
    argparser = argparse.ArgumentParser()
    argparser.formatter_class = CustomHelp
    argparser.add_argument("-l", "--path-file",
                           help="Path to file containing a list of paths.")
    argparser.add_argument("-o", "--output",
                           help="Path to output TSV file.")
    return argparser


if __name__ == "__main__":
    import sys
    parser = _setup_argparser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    main(**vars(args))
