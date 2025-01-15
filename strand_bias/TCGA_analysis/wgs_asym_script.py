
from strand_bias.pileup_parser.pileup_table import ConcatenatedPileupTable
from capture_kit_counts import CaptureKit
from strand_bias.TCGA_analysis.plotting import (
    plot_mismatch_asymmetry,
    # plot_mismatch_asym_ridgeline
    )

_dir = "/home/tyler/Documents/Projects/StrandBias/Analysis/WGS/"
wgs_exome = ConcatenatedPileupTable.read_csv(_dir + "all_samples.asym_table.tsv")
wgs_chr21 = ConcatenatedPileupTable.read_csv(_dir + "all_samples.asym_table.wgs.tsv")
vcrome = CaptureKit.read_capture_kit_nucleotide_summary("/home/tyler/Documents/Projects/StrandBias/VCRome.hg38.nt_counts.tsv")

wgs_ref_plots = plot_mismatch_asymmetry(wgs_exome, "reference", shared_yscale="box",
                                        normalize_by_nt_content=True,
                                        normalization_data=vcrome,
                                        normalization_factor=1000,
                                        seq_type="WGS")
wgs_trans_plots = plot_mismatch_asymmetry(wgs_exome, "transcription", shared_yscale="box",
                                          normalize_by_nt_content=True,
                                          normalization_data=vcrome,
                                          seq_type="WGS")

chr21_ref_plots = plot_mismatch_asymmetry(wgs_chr21, "reference", shared_yscale="box",
                                          seq_type="WGS")


# %% Paper figures.
sfig4 = wgs_ref_plots[0].update_layout(title=None, height=1000, width=1000)
sfig5 = wgs_trans_plots[0].update_layout(title=None, height=1000, width=1000)
