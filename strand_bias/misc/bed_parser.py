
import argparse

import mygene

STRAND_DICT = {-1: "-", 1: "+"}


class GeneAnnotation:
    def __init__(self, annotation_type, gene_id, strand=None):
        self.annotation_type = annotation_type
        self.gene_id = gene_id
        self.strand = strand

    def __str__(self):
        string = f"{self.annotation_type}={self.gene_id}"
        return string

    def __repr__(self):
        string = (f"GeneAnnotation("
                  f"annotation_type='{self.annotation_type}', "
                  f"gene_id='{self.gene_id}'"
                  f")")
        return string

    # def lookup_strand(self):
    #     if self.annotation_type != "ensembl_gene_id":
    #         return None
    #     mg = mygene.MyGeneInfo()
    #     results = mg.query(self.gene_id, fields="genomic_pos",
    #                        species="human", scope="ensembl.gene")
    #     if not results["hits"]:
    #         return None
    #     gen_pos = results["hits"][0]["genomic_pos"]
    #
    #     if isinstance(gen_pos, dict):
    #         return STRAND_DICT[gen_pos["strand"]]
    #     elif isinstance(gen_pos, list):
    #         strands = set()
    #         for locus in gen_pos:
    #             if locus["ensemblgene"] != self.gene_id:
    #                 continue
    #             else:
    #                 strands.add(STRAND_DICT[locus["strand"]])
    #         strands = "/".join(strands)
    #         return strands

    def lookup_strand(self, strand_dict):
        if self.annotation_type != "ensembl_gene_id":
            return None
        if self.gene_id not in strand_dict:
            return
        strand = "/".join(strand_dict[self.gene_id])
        return strand


class BedEntry:
    def __init__(self, chromosome, start, stop, annotations=None, strand=None, gene_sep: str = ";"):
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.annotations = annotations
        self.strand = strand
        self.genes = self.unpack_annotations(gene_sep)
        self._packed_annotations = True
    
    def __str__(self):
        fields = [self.chromosome, self.start, self.stop]
        if self.annotations is not None:
            fields.append(self.annotations)
        if self.strand is not None:
            fields.append(self.strand)
        string = "\t".join(fields) + "\n"
        return string

    def unpack_annotations(self, gene_sep: str = ";"):
        if self.annotations == "." or self.annotations is None:
            genes = []
        else:
            genes = [x for x in self.annotations.split(gene_sep) if x]
            genes = [x.split("=") for x in genes]
            genes = [GeneAnnotation(x, y) for x, y in genes]
        return genes

    def set_gene_strands(self, strand_dict):
        for gene in self.genes:
            gene.strand = gene.lookup_strand(strand_dict)
    
    def compile_gene_strands(self, reduce=False):
        strands = [gene.strand for gene in self.genes]
        if reduce:
            strands = sorted({strand for strand in strands if strand})
        return strands

    def set_entry_strand(self, reduce=True):
        compiled = self.compile_gene_strands(reduce)
        compiled = ",".join(compiled)
        self.strand = compiled


def parse_bed_file(bed_file: str, gene_sep: str = ";"):
    with open(bed_file) as infile:
        bed = infile.readlines()
    bed = [line.strip().split("\t") for line in bed]
    bed = [BedEntry(*line, gene_sep=gene_sep) for line in bed]
    return bed


def build_strand_dict(bed_entries: [BedEntry]):
    gene_ids = {gene.gene_id for bed in bed_entries for gene in bed.genes
                if gene.annotation_type == "ensembl_gene_id"}
    mg = mygene.MyGeneInfo()
    results = mg.querymany(gene_ids, fields="genomic_pos",
                           species="human", scope="ensembl.gene")
    strand_dict = dict()

    for result in results:
        if "notfound" in result and result["notfound"]:  # Skip notfound results.
            continue

        strands = set()
        query = result["query"]

        if isinstance(result["genomic_pos"], dict):  # Normalize types by placing dicts inside lists.
            result["genomic_pos"] = [result["genomic_pos"]]

        for locus in result["genomic_pos"]:
            if locus["ensemblgene"] != query:  # Ignore loci that don't match the query.
                continue
            strands.add(STRAND_DICT[locus["strand"]])

        if query not in strand_dict:
            strand_dict[query] = set()
        strand_dict[query] |= strands

    return strand_dict


def write_bed(bed_entries: [BedEntry], output: str):
    with open(output, "w") as outfile:
        outfile.writelines([str(entry) for entry in bed_entries])


def _setup_argparser():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-b", "--bed", dest="bed_file", help="Input BED file.")
    argparser.add_argument("-r", "--reduce-strand", action="store_true", dest="reduce",
                           default=False,
                           help="Unique-ify strand annotations.")
    argparser.add_argument("-o", "--output", help="Output BED file.")
    argparser.add_argument("-d", "--delimiter", default=";",
                           help="BED file gene list (col 4) delimiter.")
    return argparser


def _main():
    argparser = _setup_argparser()
    args = argparser.parse_args()
    beds = parse_bed_file(args.bed_file, args.delimiter)
    strand_dict = build_strand_dict(beds)
    for entry in beds:
        entry.set_gene_strands(strand_dict)
        entry.set_entry_strand(args.reduce)
    write_bed(beds, args.output)


if __name__ == '__main__':
    _main()
