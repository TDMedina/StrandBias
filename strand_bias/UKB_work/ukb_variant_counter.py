
import argparse
from pathlib import Path

import pandas as pd
from pandas import IndexSlice as idx
from pandas.api.extensions import register_dataframe_accessor
from pysam import VariantFile, VariantRecord


def _make_var_id_col(pandas_obj):
    new_col = [f"{x.chromosome}_{x.pos}_{x.mutation}"
               for x in pandas_obj.index.to_frame().itertuples(index=False)]
    return new_col

def _is_oxog_mut(record):
    return ((record.ref == "G" and record.alts[0] == "T")
            or (record.ref == "C" and record.alts[0] == "A"))


def tally_variants(vcf_path, transcribed=None):
    hets = {(0, 1), (1, 0), (1, None), (None, 1)}
    results = []
    vcf = VariantFile(vcf_path)
    for record in vcf.fetch():
        if not (_is_simple_snv(record) and _is_oxog_mut(record)):
            continue
        for sample, genotype in record.samples.items():
            alleles = genotype.allele_indices
            if 1 not in alleles:
                continue
            if len(alleles) == 1:
                genotypes = [0, 0, 1, 1]
            elif alleles == (1, 1):
                genotypes = [0, 1, 0, 1]
            elif alleles in hets:
                genotypes = [1, 0, 0, 1]
            entry = [sample.split("_")[0],
                     record.chrom, record.pos,
                     record.ref + record.alts[0],
                     transcribed] + genotypes
            results.append(entry)
    if not results:
        print(f"No variants found for '{vcf_path}'.")
        exit()
    results = pd.DataFrame(results)
    results.columns = ["sample_id", "chromosome", "pos", "mutation", "transcribed", "het", "hom", "hap", "count"]
    results = results.set_index(["sample_id", "chromosome", "pos", "mutation", "transcribed"]).sort_index()
    return results


def _is_simple_snv(variant_record: VariantRecord):
    if any([
        variant_record.rlen != 1,
        len(variant_record.alts) != 1,
         ]):
        return False
    return True


# def tally_all_vcfs(vcf_paths, export_path=None):
#     table = tally_variants(vcf_paths[0], transcribed=vcf_paths[0].split(".")[2])
#     for vcf_path in vcf_paths[1:]:
#         table = pd.concat(tally_variants(vcf_path, transcribed=vcf_path[0].split(".")[2]))
#     if export_path is not None:


# def read_vcf_paths(vcf_path_file):
#     with open(vcf_path_file) as infile:
#         vcf_paths = vcf_path_file.readlines()
#     vcf_paths = [x.rstrip() for x in vcf_paths]
#     return vcf_paths


def _setup_argparser():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", "--input-vcf", required=True, dest="vcf_path")
    argparser.add_argument("-o", "--output", required=True, dest="export_path",
                           help="Output TSV file path.")
    argparser.add_argument("-t", "--transcribed")
    return argparser


@register_dataframe_accessor("ukb_tally")
class UKBVariantTally:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(file_path):
        table = pd.read_csv(file_path, sep="\t", index_col=list(range(5)))
        return table

    def make_var_id_col(self):
        var_id_col = _make_var_id_col(self._obj)
        return var_id_col

    def identify_singleton_loci(self, consider=None):
        if consider is None:
            consider = ["count"]
        singletons = self._obj[consider].groupby(["chromosome", "pos", "mutation"]).agg(sum).sum(axis=1)
        singletons = singletons.loc[singletons == 1].to_frame()
        singletons.columns = ["count"]
        return singletons

    def filter_for_singletons(self, consider=None):
        singleton_loci = self.identify_singleton_loci(consider)
        singleton_loci["var_id"] = _make_var_id_col(singleton_loci)
        table = self._obj.copy()
        table["var_id"] = _make_var_id_col(table)
        table = table.loc[table.var_id.isin(singleton_loci.var_id)]
        return table


def main(vcf_path, export_path, transcribed=None):
    vcf_path = Path(vcf_path)
    table = tally_variants(vcf_path, transcribed)
    table.to_csv(export_path, sep="\t", index=True)
    return table


if __name__ == '__main__':
    args = _setup_argparser().parse_args()
    counts = main(**vars(args))
