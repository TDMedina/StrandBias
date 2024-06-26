
from pathlib import Path

from ukb_mismatch_table import UkbMismatchTable
from ukb_variant_table import UkbVariantsUnstacked

from plot_help import VarCounts, PlotDirectory

# %% Read data.
data_dir = Path("~/StrandBias/UKB_analysis/")
var_dir = data_dir/"variant_reanalysis"

mismatches = UkbMismatchTable.read_csv(data_dir/"asymmetry.tsv")
oxog_variants = VarCounts(
    filtered=UkbVariantsUnstacked.read_csv(var_dir/"oxog_vars.snp_filtered.stats.project_ids.tsv"),
    unfiltered=UkbVariantsUnstacked.read_csv(var_dir/"oxog_vars.unfiltered.stats.project_ids.tsv")
    )


# %% Plot batch boxplots.

plots = PlotDirectory(mismatches, oxog_variants)


def test_correlation(mismatch_table, variant_table, method="pearson"):
    pass
