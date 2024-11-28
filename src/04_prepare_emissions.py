# -*- coding: utf-8 -*-
# %% [markdown]
#  # Norwegian territorial emissions.
#
# Script that prepares the Norwegian territorial emissions
# as reported by SSB for the SNAC coupling.
#
# %% [markdown]
# Copyright (C) 2024 XIO Sustainability Analytics, Inc
#
# Originally written by
#
# - Richard Wood
# - Konstantin Stadler
#
# Updated for CaFEAN 2 by
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
# Import of required external packages. These are all available through pip/conda/mamba install.

# %%
import pandas as pd
from functions.general import aggregate_chemical_sectors 

# %% [markdown]
# Python internal packages (this don't need to be installed, part of standard python)

# %%
from pathlib import Path

# %% [markdown]
# ## Settings

# %% [markdown]
# ### Parameters
# Run the model for the following years:
#
# TODO: Ensure that columns in household emissions are named correctly!
# %% [markdown]
# ### Path definitions
# %% [markdown]
# The official Norwegian emission data file.
# %%
norwegian_emission_file: Path = (
    config.data_path 
    / "01_raw"
    / "ssb"
    / "territorial emissions"
    / "AEA Questionnaire_2023_Norway2023_til footprint.xlsm"
)

# %% [markdown]
# The path for storing the processed emissions data.
# %%
emissions_interim_path: Path = (
    config.data_path
    / "02_interim"
    / "territorial emissions"
)
emissions_interim_path.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ### Specifying/reading data sources
# This section defines all data sources needed for the calculation.
# Smaller dataset are read already here, larger datasets are read
# at the places they are needed.


# %% [markdown]
# The sector list file contains the names used for
# parsing the sectors in the emission data. It is also
# used for then renaming the sectors to the names
# used in the official economic data.
# *NOTE:* The official economic data is provided in
# industry by industry (ixi) format, but the names
# are "product names". We stay consistent to that
# and also use product names throughout the script.

# %%
sector_list = pd.read_excel(
    config.data_path
    / "00_auxiliary"
    / "classifications" 
    / "sector_list.xlsx"
)

industry_code_to_product_name_mapper = (
    sector_list
    .set_index(["IndustryCode"])["Product"]
    .str.strip()
    .dropna()
    .to_dict()
)

industry_codes = (
    list(
        industry_code_to_product_name_mapper
        .keys()
    )
)
# Industry L68  [Real estate activities (excluding imputed rents)] 
# is not available in the emissions. 
# We therefore have to add it manually.
industry_codes.remove("L68")

# %% [markdown]
# The emission characterization file contains
#   1. The conversion factors from EXIOBASE emission units/names to Norwegian emission units/names
#   2. The characterization of GHG to GWP100
# *NOTE: *: All unit conversion for emission data happens in that file!
# These units must match the units provided by the EXIOBASE and Norwegian emission data.
# Double-check these when updating the data!

# %%
charactarisation_factors_full = pd.read_excel(
    (
        config.data_path
        / "00_auxiliary"
        / "other" 
        / "emis_conv_char.xlsx"
    ),
    sheet_name="nor_char",
    index_col=0
)
charactarisation_factors = (
    charactarisation_factors_full
    .reset_index()
    .loc[:, ["stressor", "factor"]]
)

# %% [markdown]
# Small helper function used to rename the index
# of the households emissions.

# %%
def rename_household_names(x):
    if "-" in x:
        return x.strip("-").strip()
    else:
        return "Households, totals"

# %% [markdown]
# # Script
# %% [markdown]
# ## Load an process data for each emission type

# %%
household_emissions = []
industry_emissions = []

for emission_type in config.emission_types:
    raw_emission_data = (
        pd.read_excel(
            norwegian_emission_file, 
            sheet_name=emission_type,
            skiprows=3,
            index_col=[0, 1, 2, 3, 4]
        )
        .loc[:, config.years_io]
        .droplevel(0, axis=0)
    )

    household_emissions_index = (
        raw_emission_data
        .index
        .get_level_values(0)
        .str.contains("Households, totals", na=False)
    )
    household_emission = (
        raw_emission_data
        .loc[
            household_emissions_index, :
        ]
        .droplevel([0, 1, 3])
        .rename_axis(["year"], axis=1)
        .rename_axis(["sector"], axis=0)
    )
    household_emissions.append(household_emission)

    # Industry L68  [Real estate activities (excluding imputed rents)] 
    # is not available in the emissions. 
    # We therefore have to add it manually.
    industry_emission = (
        raw_emission_data
        .loc[
            industry_codes,
            : 
        ]
        .droplevel([1, 2, 3])
    )
    industry_emission.loc["L68", :] = 0

    industry_emission = (
        industry_emission

        .rename_axis(["sector"], axis=0)
        .rename_axis(["year"], axis=1)
    )

    industry_emission = aggregate_chemical_sectors(
        industry_emission, axis=0, level=0
    )
    industry_emissions.append(industry_emission)


# %% [markdown]
# ## Concatenate data and calculate emissions in GHG equivalents.
# %%
household_emissions_df = (
    pd.concat(
        household_emissions,
        keys=config.emission_types,
        names=["stressor"]
    )
    .stack("year")
)


household_emissions_df = (
    household_emissions_df
    .to_frame("emissions")
    .reset_index()
    .merge(
        charactarisation_factors,
        on="stressor",
        how="left"
    )
    .set_index(household_emissions_df.index.names)
)

household_emissions_df["emissions (GHG eqv.)"] = (
    household_emissions_df["emissions"]
    .mul(household_emissions_df["factor"])
)
household_ghg_df = (
    pd.concat(
        {"GHG": (
            household_emissions_df["emissions (GHG eqv.)"]
            .groupby(["sector", "year"]).sum()
        )},
        names=["stressor"]
    )
)
household_emissions_df = (
    pd.concat(
        [
            household_ghg_df,
            household_emissions_df["emissions"]
        ], axis=0
    )
    .unstack("sector")
    .reorder_levels([1, 0])
    .sort_index()
    .rename(industry_code_to_product_name_mapper, axis=1)
    .rename(rename_household_names, axis=1)
)

industry_emissions_df = (
        pd.concat(
        industry_emissions,
        keys=config.emission_types,
        names=["stressor"]
    )
    .stack("year")
)
industry_emissions_df = (
    industry_emissions_df
    .to_frame("emissions")
    .reset_index()
    .merge(
        charactarisation_factors,
        on="stressor",
        how="left"
    )
    .set_index(industry_emissions_df.index.names)
)

industry_emissions_df["emissions (GHG eqv.)"] = (
    industry_emissions_df["emissions"]
    .mul(industry_emissions_df["factor"])
)
industry_ghg_df = (
    pd.concat(
        {"GHG": (
            industry_emissions_df["emissions (GHG eqv.)"]
            .groupby(["sector", "year"]).sum()
        )},
        names=["stressor"]
    )
)
industry_emissions_df = (
    pd.concat(
        [
            industry_ghg_df,
            industry_emissions_df["emissions"]
        ], axis=0
    )
    .unstack("sector")
    .reorder_levels([1, 0])
    .sort_index()
    .rename(industry_code_to_product_name_mapper, axis=1)
)

# %% [markdown]
# ## Sort out units
# Next we make a summary dataframe holding the units for each emission type.
# We also order that and use this order during saving the results.

# %%
emission_units = (
    charactarisation_factors_full
    .loc[:, "stressor_unit"]
)
emission_units.loc["GHG"] = (
    charactarisation_factors_full["impact_unit"]
    .unique()[0]
)

emission_units = emission_units.sort_values()

# %% [markdown]
# ## Save results
# %%
# ? Should probably be saved in a single Excel file.
household_emissions_df.to_excel(
    emissions_interim_path / "final demand emissions.xlsx"
)
industry_emissions_df.to_excel(
    emissions_interim_path / "industry emissions.xlsx"
)
emission_units.to_excel(
    emissions_interim_path / "units.xlsx"
)
