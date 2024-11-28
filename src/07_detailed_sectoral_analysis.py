# -*- coding: utf-8 -*-
# %% [markdown]
#  # Detailed analysis of sectors
#
# Script that uses the deflators to provide a better understanding of 
# the changes in the consumption-based emissions accounts for each sector
# based on the changes in the emission intensities of the sectors in its supply chain.
# It uses deflators to produce constant price stressors that can be used alongside
# the current price stressors to interpret the numbers. 
#
# It also provides a function to analyse specific sectors. 
# Two examples are provided at the end of the script.

# %% [markdown]
# Copyright (C) 2024 XIO Sustainability Analytics, Inc
#
# Written by
#
# - Kajwan Rasul
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# %% [markdown]
# Import config file.
# %%
import config

# %% [markdown]
# Import own functions
# %%
from functions.general import (
    aggregate_chemical_sectors,
    diagonalise_series,
    analyse_sector,
    make_results
)

# %% [markdown]
# Import of required external packages. These are all available through pip/conda/mamba install.
# %%
import pandas as pd
import numpy as np
import pymrio

# %% [markdown]
# Python internal packages (this don't need to be installed, part of standard python)
# %%
from pathlib import Path
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


# %% [markdown]
# ## Settings
# %% [markdown]
# ### Locations / folder definitions

# %% [markdown]
# Here we specify the path to the processed norwegian IO data
# in current prices and the territorial emissions

# %%
norwegian_IO_interim_path: Path = (
    config.data_path
    / "02_interim" 
    / "norwegian IO"
)

norwegian_territorial_emissions_path: Path = (
    config.data_path
    / "02_interim"
    / "territorial emissions"
)
deflators_interim_path: Path = (
    config.data_path
    / "02_interim" 
    / "deflators"
)

# %% [markdown]
# Here we specify the path to the EXIOBASE data.
# %%
exiobase_interim_path: Path = (
    config.data_path
    / "02_interim"
    / "exiobase"
)

# %% [markdown]
# Specify the path the results are to be saved.
# %%
deflated_results_path: Path = (
    config.data_path
    / "03_results"
    / "Detailed analysis"
)
deflated_results_path.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Load data
# %%
sector_matching = pd.read_excel(
    (
        config.data_path
        / "00_auxiliary"
        / "concordances"
        / "exiobase_and_io.xlsx"
    ),
    index_col=[0, 1, 2, 3],
    header=[0, 1, 2, 3]
)
exiobase_indudstry_mapper = (
    dict(zip(
        sector_matching.index.get_level_values(2),
        sector_matching.index.get_level_values(3),
    ))
)
exiobase_norway_concordance = (
    sector_matching
    .droplevel([0, 1, 2], axis=0)
    .droplevel([0, 1, 3], axis=1)
)

# %%
sector_aggregator = (
    pd.read_csv(
        (
            config.data_path
            / "00_auxiliary" 
            / "classifications" 
            / "sector_agg_spec.tsv"
        ),
        sep="\t"
    )
    .set_index("src")
    .rename(lambda x: x.strip())
    .loc[:, "agg"]
    .to_dict()
)

# %%
exchange_rates = (
    pd.read_csv(
        (
            config.data_path
            / "00_auxiliary"
            / "other"
            / "exchange_rates.csv"
        ),
        sep="\t"
    )
    .set_index(["year"])
    .loc[:, "euro_to_nok"]
)

# %%
final_demand_emissions = (
    pd.read_excel(
        norwegian_territorial_emissions_path
        / "final demand emissions.xlsx",
        index_col=[0, 1]
    )
    .rename_axis(["year", "account"], axis=0)
)
industry_emissions = (
    pd.read_excel(
        norwegian_territorial_emissions_path
        / "industry emissions.xlsx",
        index_col=[0, 1]
    )
    .rename_axis(["year", "account"], axis=0)  # We should instead be consistent accross all datasets.
)
industry_emissions = aggregate_chemical_sectors(industry_emissions, axis=1, level=0)

# %%
import_deflators = (
    pd.read_excel(
        deflators_interim_path / "import.xlsx",
        index_col=[0, 1]
    ).iloc[:, 0]
)

gross_domestic_deflators = (
    pd.read_excel(
        deflators_interim_path / "domestic.xlsx",
        index_col=[0, 1]
    ).iloc[:, 0]
)

# %%
# Prepare placeholders for all results of interest
F_domestic_list = []
F_import_list = []
F_import_exiobase_list = []
x_domestic_list = []
x_domestic_constant_list = []
x_import_list = []
x_import_constant_list = []
x_import_exiobase_list = []
x_import_exiobase_constant_list = []
S_domestic_list = []
S_domestic_constant_list = []
S_import_list = []
S_import_constant_list = []

import_footprint_list = []
domestic_footprint_list = []

for year in config.years_io:
    # Load Norwegian domestic IO data
    Z_domestic = pd.read_excel(
        norwegian_IO_interim_path / str(year) / "Z_domestic.xlsx",
        index_col=[0]
    )

    Z_import = pd.read_excel(
        norwegian_IO_interim_path / str(year) / "Z_import.xlsx",
        index_col=[0]
    )

    Y_domestic_breakdown = pd.read_excel(
        norwegian_IO_interim_path / str(year) / "Y_domestic.xlsx",
        index_col=[0]
    )
    Y_domestic_breakdown[Y_domestic_breakdown < 0] = 0

    Y_domestic_export = (
        Y_domestic_breakdown
        .loc[:, "Exports fob (2)"]
    )

    Y_domestic = (
        Y_domestic_breakdown
        .drop("Exports fob (2)", axis=1)
        .sum(axis=1)
    )

    Y_import_breakdown = pd.read_excel(
        norwegian_IO_interim_path / str(year) / "Y_import.xlsx",
        index_col=[0]
    )
    Y_import_breakdown[Y_import_breakdown < 0] = 0

    Y_import_export = (
        Y_import_breakdown
        .loc[:, "Exports fob (2)"]
    )

    Y_import = (
        Y_import_breakdown
        .drop("Exports fob (2)", axis=1)
        .sum(axis=1)
    )

    # Calculating x after removing negatives
    x_domestic = (
        Z_domestic.sum(axis=1) 
        + Y_domestic_breakdown.sum(axis=1)
    )
    x_domestic_constant = (
        x_domestic
        .mul(
            gross_domestic_deflators.loc[year],
        )
    )
    x_import = (
        Z_import.sum(axis=1)
        + Y_import_breakdown.sum(axis=1)
    )
    x_import_constant = (
        x_import
        .mul(
            import_deflators.loc[year],
        )
    )

    # Calculate intermediate matrices
    A_domestic = pymrio.calc_A(Z_domestic, x_domestic)
    A_import = pymrio.calc_A(Z_import, x_domestic)
    L_domestic = pymrio.calc_L(A_domestic)

    F_domestic = industry_emissions.loc[year]
    S_domestic = pymrio.calc_S(F_domestic, x_domestic)
    S_domestic_constant = (
        pymrio.calc_S(F_domestic, x_domestic_constant)
    )

    SL_domestic = pymrio.calc_M(S=S_domestic, L=L_domestic) 

    # Load EXIOBASE intermediate results.
    exiobase_import_data = (
        pd.read_excel(
            exiobase_interim_path / "exio_no_bp_raw.xlsx",
            index_col=[0],
            sheet_name=f"imp_{year}",
            header=[0, 1]
        )
        .rename(exiobase_indudstry_mapper, axis=0)
        .rename_axis(["variable", "unit"], axis=1)
        .sort_index(axis=1, level="unit")
    )

    exiobase_import_data.loc[:, ("imports", "M.NOK")] = (
        exiobase_import_data
        .loc[:, ("imports", "M.EUR")]
        .mul(exchange_rates.loc[year])
    )

    # Make weighted concordance and use it to adjust the import data
    weighted_concordance = (
        exiobase_norway_concordance
        .mul(
            exiobase_import_data.loc[:, ("imports", "M.EUR")],
            axis=0
        )
    )
    weighted_concordance_normalised = (
        weighted_concordance
        .div(
            weighted_concordance.sum(axis=1),
            axis=0
        )
        .replace([np.inf, -np.inf, np.nan], 0)
    )

    import_data = (
        weighted_concordance_normalised.T
        .dot(exiobase_import_data)
    )
    import_data = (
        aggregate_chemical_sectors(import_data, axis=0, level=0)
    )

    # Create relevant data from EXIOBASE
    x_import_exiobase = (
        import_data.loc[:, "imports"]
    )

    x_import_exiobase_constant = (
        x_import_exiobase
        .mul(import_deflators.loc[year], axis=0)
    )

    F_import_exiobase = (
        import_data
        .drop("imports", axis=1, level=0)
    )

    S_import = (
        F_import_exiobase
        .div(
            x_import_exiobase.loc[:, "M.NOK"],
            axis=0
        )
        .replace([np.nan, np.inf, -np.inf], 0)
    )

    F_import = (
        S_import
        .mul(x_import, axis=0)
    )

    # Missing deflator for "Water transport services" sets that to infinity (zero).
    S_import_constant = (
        F_import_exiobase
        .div(
            x_import_exiobase_constant.loc[:, "M.NOK"],
            axis=0
        )
        .replace([np.nan, np.inf, -np.inf], 0)
    )

    # Store intermediate results
    F_domestic_list.append(F_domestic)
    F_import_list.append(F_import)
    F_import_exiobase_list.append(F_import_exiobase)
    x_domestic_list.append(x_domestic)
    x_domestic_constant_list.append(x_domestic_constant)
    x_import_list.append(x_import)
    x_import_constant_list.append(x_import_constant)
    x_import_exiobase_list.append(x_import_exiobase)
    x_import_exiobase_constant_list.append(x_import_exiobase_constant)
    S_domestic_list.append(S_domestic)
    S_domestic_constant_list.append(S_domestic_constant)
    S_import_list.append(S_import)
    S_import_constant_list.append(S_import_constant)

    # Calculate footprints
    domestic_to_domestic_footprint = (
        diagonalise_series(S_domestic.loc[config.substance, :])
        .dot(L_domestic)
        .dot(diagonalise_series(Y_domestic))
    )
    import_to_domestic_footprint = (
        diagonalise_series(S_import.loc[:, (config.substance, config.substance_unit)])
        .dot(A_import)
        .dot(L_domestic)
        .dot(diagonalise_series(Y_domestic))
    )
    import_direct_footprint = (
        diagonalise_series(S_import.loc[:, (config.substance, config.substance_unit)])
        .dot(diagonalise_series(Y_import))
    )
    import_footprint = (
        import_to_domestic_footprint
        + import_direct_footprint
    )


    import_footprint_list.append(
        import_footprint
        .rename_axis(["Source sector"], axis=0)
        .rename_axis(["Final sector"], axis=1)
    )
    domestic_footprint_list.append(
        domestic_to_domestic_footprint
        .rename_axis(["Source sector"], axis=0)
        .rename_axis(["Final sector"], axis=1)
    )


# %%
# Combines detailed data on sectors in the supply chain.
# %%
F_domestic = (
    pd.concat(
        F_domestic_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis("Source sector", axis=1)
    .unstack("Year")
    .loc[config.substance]
    .reorder_levels(["Year", "Source sector"])
    .to_frame("Source sector: F")
)

x_domestic = (
    pd.concat(
        x_domestic_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Year", "Source sector"])
    .to_frame("Source sector: x (current)")
)

x_domestic_constant = (
    pd.concat(
        x_domestic_constant_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Year", "Source sector"])
    .to_frame("Source sector: x (constant)")
)

F_import = (
    pd.concat(
        F_import_list,
        keys=config.years_io,
        names=["Year"]
    )
    .loc[:, (config.substance, config.substance_unit)]
    .rename_axis(["Year", "Source sector"])
    .to_frame("Source sector: F")
)

x_import = (
    pd.concat(
        x_import_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Year", "Source sector"])
    .to_frame("Source sector: x (current)")
)

x_import_constant = (
    pd.concat(
        x_import_constant_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Year", "Source sector"])
    .to_frame("Source sector: x (constant)")
)

x_import_exiobase = (
    pd.concat(
        x_import_exiobase_list, 
        keys=config.years_io,
        names=["Year"]
    )
    .loc[:, "M.NOK"]
    .rename_axis(["Year", "Source sector"])
    .to_frame("Source sector: x (current EXIOBASE)")
)

x_import_exiobase_constant = (
    pd.concat(
        x_import_exiobase_constant_list,
        keys=config.years_io,
        names=["Year"]
    )
    .loc[:, "M.NOK"]
    .rename_axis(["Year", "Source sector"])
    .to_frame("Source sector: x (constant EXIOBASE)")
)

# %%
# Combine the data
source_sector_domestic_data = (
    pd.concat(
        [
            F_domestic,
            x_domestic,
            x_domestic_constant
        ],
        axis=1
    )
)
source_sector_import_data = (
    pd.concat(
        [
            F_import,
            x_import,
            x_import_constant,
            x_import_exiobase,
            x_import_exiobase_constant
        ],
        axis=1
    )
)

source_sector_data = (
    pd.concat(
        [
            source_sector_domestic_data,
            source_sector_import_data
        ],
        keys=[
            "Domestic",
            "Import"
        ],
        names=["Variable"],
        axis=0
    )
    .rename_axis(["Indicator"], axis=1)
)

# %%
# Calculate relative change numbers

source_sector_data_year_pivot = (
    source_sector_data
    .stack("Indicator")
    .unstack("Year")
)

source_sector_change_data = (
    source_sector_data_year_pivot
    .subtract(
        source_sector_data_year_pivot.loc[:, config.base_year], axis=0
    )
    .div(
        source_sector_data_year_pivot.loc[:, config.base_year], axis=0
    )
    .replace([np.nan, np.inf, -np.inf], 0)
    .stack("Year")
    .unstack("Indicator")
    .rename(lambda x: f"Percantage change in source sector: {x.split(':')[1]}", axis=1)
    .mul(100)
    .round(2)
)

# %%
# Make stressors and calculate relative changes
S_domestic_ts = (
    pd.concat(
        S_domestic_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Source sector"], axis=1)
    .unstack("Year")
    .loc[config.substance, :]
    .reorder_levels(["Year", "Source sector"])
)
S_domestic_constant_ts = (
    pd.concat(
        S_domestic_constant_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Source sector"], axis=1)
    .unstack("Year")
    .loc[config.substance, :]
    .reorder_levels(["Year", "Source sector"])
)

S_import_ts = (
    pd.concat(
        S_import_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Year", "Source sector"], axis=0)
    .loc[:, (config.substance, config.substance_unit)]
)

S_import_constant_ts = (
    pd.concat(
        S_import_constant_list,
        keys=config.years_io,
        names=["Year"]
    )
    .rename_axis(["Year", "Source sector"], axis=0)
    .loc[:, (config.substance, config.substance_unit)]
)

S_current = (
    pd.concat(
        [
            S_domestic_ts,
            S_import_ts
        ],
        keys=["Domestic", "Import"],
        names=["Variable"],
        axis=0
    )
    
)
S_current_absolute_difference = (
    S_current
    .unstack("Year")
    .subtract(
        S_current.unstack("Year").loc[:, config.base_year],
        axis=0
    )
)
S_current_relative_change = (
    S_current_absolute_difference
    .div(
        S_current.unstack("Year").loc[:, config.base_year],
        axis=0
    )
    .replace([np.nan, np.inf, -np.inf], 0)
    .stack("Year")
    .mul(100)
    .round(2)
)

S_constant = (
    pd.concat(
        [
            S_domestic_constant_ts,
            S_import_constant_ts
        ],
        keys=["Domestic", "Import"],
        names=["Variable"],
        axis=0
    )
)
S_constant_absolute_difference = (
    S_constant
    .unstack("Year")
    .subtract(
        S_constant.unstack("Year").loc[:, config.base_year],
        axis=0
    )
)
S_constant_relative_change = (
    S_constant_absolute_difference
    .div(
        S_constant.unstack("Year").loc[:, config.base_year],
        axis=0
    )
    .replace([np.nan, np.inf, -np.inf], 0)
    .stack("Year")
    .mul(100)
    .round(2)
)

# %%
# Combine results
S_results = (
    pd.concat(
        [
            S_current_relative_change,
            S_constant_relative_change
        ],
        keys=[
            f"Percentage change in stressor (current prices) from {config.base_year}",
            f"Percentage change in stressor (constant prices) from {config.base_year}"
        ],
        names=["Indicator"],
        axis=1
    )
)

import_footprint_ts = (
    pd.concat(
        import_footprint_list,
        keys=config.years_io,
        names=["Year"]
    )
)

domestic_footprint_ts = (
        pd.concat(
        domestic_footprint_list,
        keys=config.years_io,
        names=["Year"]
    )
)
# %%
# Make result datasets
results = make_results(
    year=config.base_year,
    domestic_footprint_ts=domestic_footprint_ts,
    import_footprint_ts=import_footprint_ts,
    S_results=S_results,
    source_sector_data=source_sector_data,
    source_sector_change_data=source_sector_change_data,
    sector_aggregator=sector_aggregator,
    aggregate_sectors=False
)

results_agg = make_results(
    year=config.base_year,
    domestic_footprint_ts=domestic_footprint_ts,
    import_footprint_ts=import_footprint_ts,
    S_results=S_results,
    source_sector_data=source_sector_data,
    source_sector_change_data=source_sector_change_data,
    sector_aggregator=sector_aggregator,
    aggregate_sectors=True
)


# %%
# Save results
results.reset_index().to_excel(deflated_results_path / "sectoral_analysis.xlsx")
results_agg.reset_index().to_excel(deflated_results_path / "sectoral_analysis_agg.xlsx")

# %%
# Example 1: Construction sector analysis
result_construction = (
    analyse_sector(
        year=2021,
        final_sector=1,
        dataset=results_agg
    )
)

result_construction.sort_values(by="Footprint", ascending=False).round(2)

# %%
# Example 2: Business Services analysis
result_business =  (
    analyse_sector(
        year=2021,
        final_sector="Business Services",
        dataset=results_agg
    )
)

result_business.sort_values(by="Footprint", ascending=False).round(2)
