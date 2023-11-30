# -*- coding: utf-8 -*-
# %% [markdown]
#  # Small utility for converting the result files to excel
#
# The data converted here can then be copy paste data into corresponding .xlsx file from
# the original report to reproduce the result figures.
#
# NOTE: This script must run after the `aggregater.py` script has been run.
#
# Results are stored at results/excel_report_format
#
# Requires pandas
#

# %% [markdown]
# Import required libraries

# %%
from pathlib import Path
import pandas as pd

# %% [markdown]
# Set the work path to the directory where this script is located
# and if this is not available to the current working directory

# %%
try:
    work_path = Path(__file__).parent.absolute()  # when running as script
except NameError:
    work_path = Path.cwd()

# %% [markdown]
# Folder definitions

# %%
results_folder = work_path / "results"

report_folder = results_folder / "excel_report_format"
report_folder.mkdir(exist_ok=True)

# %% [markdown]
# Read result files and convert to table format

# %%
_tot_all = pd.read_csv(results_folder / "total_accounts.tsv", sep="\t")
ind_col = [c for c in _tot_all.columns if c != "value"]
all_totals = _tot_all.set_index(ind_col, append=False, drop=True).squeeze()

for emis in all_totals.index.get_level_values("emission").unique():
    tot = all_totals.loc[
        all_totals.index.get_level_values("emission") == emis, :
    ].unstack("year")
    if "Biomass" in emis:
        filename = report_folder / "bio_totals.tsv"
    else:
        filename = report_folder / f"{emis}_totals.tsv"
    tot.to_csv(filename, sep="\t")

all_sources_agg = pd.read_csv(results_folder / "footprint_sources_agg.tsv", sep="\t")
ghg_source_reg_agg = (
    all_sources_agg.query("emission=='GHG'").query("year==2020").drop("year", axis=1)
)
ghg_source_reg_agg.set_index(
    ["emission", "unit", "sector", "region"], inplace=True, append=True
)
ghg_source_reg_agg = ghg_source_reg_agg.droplevel(0, axis=0)
ghg_source_reg_agg_pivot = ghg_source_reg_agg.unstack("region")
ghg_source_reg_agg_pivot = ghg_source_reg_agg_pivot.droplevel(0, axis=1)
ghg_source_reg_agg_pivot.to_csv(report_folder / "ghg_source_reg_agg_2020.tsv", sep="\t")


all_sector_breakdown_agg = pd.read_csv(
    results_folder / "sector_accounts_agg.tsv", sep="\t"
)
for loop_component in all_sector_breakdown_agg.component.unique():
    ghg_sector = (
        all_sector_breakdown_agg.query("year==2020")
        .query("component==@loop_component")
        .drop("year", axis=1)
        .drop("component", axis=1)
        .drop("unit", axis=1)
    )
    ghg_sector.set_index(["emission", "sector"], inplace=True, append=True)
    ghg_sector = ghg_sector.droplevel(0, axis=0)
    ghg_sector_pivot = ghg_sector.unstack("emission")
    ghg_sector_pivot = ghg_sector_pivot.droplevel(0, axis=1)
    ghg_sector_pivot = ghg_sector_pivot[
        ["GHG", "CO2", "Biomass CO2", "N2O", "SF6_NF3", "CH4", "HFC", "PFC"]
    ]
    output_file_path = report_folder / (
        "sector_accounts_agg_2020_" + str(loop_component) + ".tsv"
    )
    ghg_sector_pivot.to_csv(output_file_path, sep="\t")


all_findem_breakdown_agg = pd.read_csv(
    results_folder / "footprints_final_demand_breakdown_agg.tsv", sep="\t"
)
ghg_findem = (
    all_findem_breakdown_agg.query("emission=='GHG'")
    .query("year==2020")
    .query("component=='total_footprint'")
    .drop("year", axis=1)
    .drop("component", axis=1)
)
ghg_findem.set_index(
    ["emission", "unit", "sector", "category"], inplace=True, append=True
)
ghg_findem = ghg_findem.droplevel(0, axis=0)
ghg_findem_pivot = ghg_findem.unstack("category")
ghg_findem_pivot = ghg_findem_pivot.droplevel(0, axis=1)
output_file_path = report_folder / (
    "sector_accounts_agg_GHG_final_demand_breakdown_2020.tsv"
)
ghg_findem_pivot.to_csv(output_file_path, sep="\t")
