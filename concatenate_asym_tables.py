
import argparse
import pandas as pd


def read_pileup_table(file, wgs=False):
    index_col = [0, 1]
    if not wgs:
        index_col.append(2)
    table = pd.read_csv(file, sep="\t", index_col=index_col, header=[0, 1])
    return table


def read_all_tables(file_list, wgs=False, output=None):
    with open(file_list) as infile:
        files = infile.readlines()
    files = [file.rstrip() for file in files]
    table = read_pileup_table(files[0], wgs)
    for file in files[1:]:
        table = pd.concat([table, read_pileup_table(file, wgs)])
    if output:
        table.to_csv(output, sep="\t", index=True)
    return table


def _setup_argparser():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-f", "--file-list", required=True)
    argparser.add_argument("-o", "--output")
    argparser.add_argument("-w", "--wgs", action="store_true", default=False)
    return argparser


if __name__ == '__main__':
    import sys
    argparser = _setup_argparser()
    if len(sys.argv) == 1:
        argparser.print_help()
        sys.exit()

    args = argparser.parse_args()
    result_table = read_all_tables(**vars(args))
