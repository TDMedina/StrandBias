# Strand Bias Asymmetry Tools

This repository contains tools related to the article "Consistent asymmetry in DNA damage artefacts across target regions in exome sequencing data", which details observed asymmetry in whole-exome sequencing alignment mismatches against the reference genome.

Scripts and tools related to the generation of pileup files and mismatch counts are found in `count_pipeline/` and `count_pipeline_wgs/`.

Tools related to parsing and analyzing these counts are found `strand_bias/pileup_parser/`.

For scripts to reproduce the analysis from the related article, see `TCGA_analysis/tcga_main_analysis.py` and `UKB_work/ukb_main_analysis.py`.
