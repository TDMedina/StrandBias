
import argparse
from collections import Counter, defaultdict

import mygene
import pandas as pd
from pandas import IndexSlice as idx
from pysam import VariantFile, VariantRecord

from gtf_parser import read_gtf_file


_MYGENESERVICE = mygene.MyGeneInfo()


def _overlap(range1, range2):
    if range1.start <= range2.stop-1 and range2.start <= range1.stop-1:
        return True
    return False


def make_count_table(project_id=None, case_id=None, file_id=None):
    nts = [] if project_id is None else list("ACGT")
    project_id = [] if project_id is None else [project_id]
    case_id = [] if case_id is None else [case_id]
    file_id = [] if file_id is None else [file_id]

    counts = pd.DataFrame(0,
                          index=pd.MultiIndex.from_product(
                              [project_id, case_id, file_id, nts],
                              names=["project_id", "case_id", "file_id", "reference"]
                              ),
                          columns=pd.MultiIndex.from_product(
                              [["forward_coding", "reverse_coding"], list("ACGT"), ["PASS", "FAIL"]],
                              names=["coding_region", "alt", "filtered"]))
    return counts


def _is_simple_snv(variant_record: VariantRecord):
    if any([variant_record.rlen != 1,
            len(variant_record.alts) != 1,
            variant_record.alts[0] not in variant_record.samples[0].alleles]):
        return False
    return True


def count_region_snv_types(vcf_path, coding_region, project_id, case_id, file_id):
    counts = make_count_table(project_id, case_id, file_id)
    vcf = VariantFile(vcf_path)
    for record in vcf.fetch():
        if not _is_simple_snv(record):
            continue
        alt = record.alts[0]
        filtered = "PASS" if set(record.filter) == {"PASS"} else "FAIL"
        counts.loc[idx[project_id, case_id, file_id, record.ref],
                   idx[coding_region, alt, filtered]] += 1
    return counts


def count_snv_types(vcf_forward, vcf_reverse, project_id, case_id, file_id):
    counts = (count_region_snv_types(vcf_forward, "forward_coding", project_id, case_id, file_id)
              + count_region_snv_types(vcf_reverse, "reverse_coding", project_id, case_id, file_id))
    return counts


def _make_file_paths(sample):
    file_base = sample.file_name.replace(".vcf.gz", "")
    file_base = "/".join([".", "samples", sample.case_id, sample.file_id, file_base])
    vcf_forward = file_base + ".tumor_snvs.forward.vcf.gz"
    vcf_reverse = file_base + ".tumor_snvs.reverse.vcf.gz"
    return vcf_forward, vcf_reverse


def count_all_sample_snv_types(sample_table: pd.DataFrame):
    counts = make_count_table()
    total = sample_table.shape[0]
    for i, sample in sample_table.iterrows():
        print(f"Counting files: {i}/{total}.\r", end="")
        vcf_forward, vcf_reverse = _make_file_paths(sample)
        sample_counts = count_snv_types(vcf_forward, vcf_reverse, sample.project_id, sample.case_id, sample.file_id)
        counts = pd.concat([counts, sample_counts])
    return counts


def tally_variant_information(vcf_path, sample_info, coding_region):
    vcf = VariantFile(vcf_path)
    counts = {(sample_info.project_id, sample_info.case_id, sample_info.file_id,
               record.contig, record.pos, coding_region,
               record.ref, record.alts[0]):
              Counter(record.filter)
              for record in vcf.fetch() if _is_simple_snv(record)}
    for counter in counts.values():
        counter["TOTAL"] = 1
        if counter["PASS"] == 0:
            counter["FAIL"] = 1
    return counts


def fetch_gene_symbol(contig, position):
    query = f"{contig}:{position}-{position+1}&species:human"
    fields = ["symbol", "ensembl.gene"]
    results = _MYGENESERVICE.query(query, fields=fields, as_dataframe=True)
    results.loc[pd.isna(results["ensembl.gene"]), "ensembl.gene"] = ""
    ensg_results = results.loc[results["ensembl.gene"].str.startswith("ENSG")]
    if ensg_results.empty:
        symbols = results.symbol[0]
    else:
        symbols = ", ".join([str(symbol) for symbol in ensg_results.symbol])
    return symbols


def count_all_variants_in_aggro(sample_table: pd.DataFrame, gtf_path=None):
    all_counts = dict()
    # all_counts = defaultdict(Counter)
    total = sample_table.shape[0]
    for i, sample in enumerate(sample_table.itertuples(), start=1):
        print(f"Counting files: {i}/{total}.\r", end="")
        file_paths = _make_file_paths(sample)
        for path, region in zip(file_paths, ["forward", "reverse"]):
            variant_counts = tally_variant_information(path, sample, region)
            all_counts |= variant_counts
            # for id_info, count in variant_counts.items():
            #     all_counts[id_info] += count
            #     all_counts[id_info]["total"] += 1
    print("\nMaking count table...")
    all_counts = pd.DataFrame.from_dict(all_counts, orient="index")
    all_counts.index.names = ["project_id", "case_id", "file_id",
                              "contig", "pos",
                              "gene_orientation", "ref", "alt"]
    if gtf_path is not None:
        print("Assigning gene symbols...")
        symbols = _lookup_gene_positions(all_counts, gtf_path)
        all_counts["genes"] = symbols
    # symbols = [fetch_gene_symbol(row.contig, row.pos)
    #            for row in all_counts.reset_index().itertuples()]
    # all_counts["gene"] = symbols
    # print("Done.")
    return all_counts


def _lookup_gene_positions(count_table: pd.DataFrame, gtf_path):
    gene_info = read_gtf_file(gtf_path)
    symbols = [", ".join(_lookup_gene_position(row.contig, row.pos, gene_info))
               for row in count_table.reset_index().itertuples()]
    return symbols


def _lookup_gene_position(contig, position, gene_table):
    contig_results = gene_table.loc[contig]
    contig_results = contig_results.loc[(contig_results.start <= position)
                                        & (position <= contig_results.stop)]
    symbols = set(contig_results.gene_name)
    # symbols = [gene.gene_name for gene in bin_results.itertuples()
    #            if _overlap(lookup_range, range(gene.start, gene.stop+1))]
    return symbols


def read_all_allele_frequencies(sample_table: pd.DataFrame):
    all_freqs = dict()
    total = sample_table.shape[0]
    for i, sample in enumerate(sample_table.itertuples(), start=1):
        print(f"Counting files: {i}/{total}.\r", end="")
        file_paths = _make_file_paths(sample)
        for path, region in zip(file_paths, ["forward", "reverse"]):
            freqs = read_allele_frequency(path, sample, region)
            all_freqs |= freqs
    print("\nMaking frequency table...")
    all_freqs = pd.DataFrame.from_dict(all_freqs, orient="index")
    all_freqs.index = pd.MultiIndex.from_tuples(all_freqs.index)
    all_freqs.index.names = ["project_id", "case_id", "file_id",
                             "contig", "pos",
                             "gene_orientation", "ref", "alt"]
    all_freqs.columns = ["frequency", "ref_depth", "alt_depth"]
    return all_freqs


def read_allele_frequency(vcf_path, sample_info, coding_region):
    vcf = VariantFile(vcf_path)
    freqs = {(sample_info.project_id, sample_info.case_id, sample_info.file_id,
              record.contig, record.pos, coding_region,
              record.ref, record.alts[0]):
             [(geno := record.samples[0])["AF"], geno["AD"][0], geno["AD"][1]]
             for record in vcf.fetch() if _is_simple_snv(record)}
    return freqs


def read_sample_table(sample_table_path):
    samples = pd.read_csv(sample_table_path, sep="\t")
    return samples


def _setup_argparser():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-s", "--sample-table", required=True)
    argparser.add_argument("-o", "--output", help="Output TSV file path.")
    argparser.add_argument("-a", "--in-aggro", help="Tally variants by position, aggregating across samples.",
                           action="store_true")
    argparser.add_argument("-f", "--read-frequencies",
                           help="Read variant frequencies by position for all samples.",
                           action="store_true")
    argparser.add_argument("-g", "--gtf-path", help="Path to GTF file of gene information.")
    return argparser


def main(sample_table, output=None, in_aggro=False, read_frequencies=False, gtf_path=None):
    sample_table = read_sample_table(sample_table)
    if in_aggro:
        table = count_all_variants_in_aggro(sample_table, gtf_path=gtf_path)
    elif read_frequencies:
        table = read_all_allele_frequencies(sample_table)
    else:
        table = count_all_sample_snv_types(sample_table)
    if output is not None:
        table.to_csv(output, sep="\t", index=True)
    return table


if __name__ == '__main__':
    args = _setup_argparser().parse_args()
    counts = main(**vars(args))
