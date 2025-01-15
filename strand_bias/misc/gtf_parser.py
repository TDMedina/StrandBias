
import gzip
import pandas as pd


def read_gtf_file(gtf_path):
    with open(gtf_path) as infile:
        gene_data = infile.readlines()
    gene_data = [line.rstrip("\n").split("\t") for line in gene_data]
    for i, entry in enumerate(gene_data):
        end_dict = entry[-1].rstrip("; ").split(";")
        end_dict = [x.strip().split(" ") for x in end_dict]
        end_dict = dict([[x.strip('" ') for x in y] for y in end_dict])
        entry = {"contig": entry[0], "start": int(entry[3]),
                 "stop": int(entry[4]), "orientation": entry[6]} | end_dict
        gene_data[i] = entry
    gene_data = pd.DataFrame.from_dict(gene_data)
    # gene_data["bin"] = gene_data.start // 10000 * 10000
    gene_data.sort_values(by=["start", "stop"], inplace=True)
    gene_data.set_index("contig", inplace=True)
    gene_data.sort_index(inplace=True)
    return gene_data
