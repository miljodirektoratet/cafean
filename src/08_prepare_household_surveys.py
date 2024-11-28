# -*- coding: utf-8 -*-
# %% [markdown]
#  # Household survey footprint preparation script
#
# This script disaggregates the household footprint calculated by the SNAC
# model using SSB's household consumption surveys at the highest level of resolution.
#
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
# Import of required external packages. These are all available through pip/conda/mamba install.

# %%
import numpy as np
import pandas as pd
from functions.general import load_survey
from functions.gras import gras

# %% [markdown]
# Python internal packages (this don't need to be installed, part of standard python)

# %%
from pathlib import Path

# %% [markdown]
# The household consumption surveys and related statistics.
# %%
surveys_path: Path = ( 
    config.data_path 
    / "01_raw"
    / "ssb"
    / "household consumption surveys"
)

household_statistics_paths: Path = (
    config.data_path
    / "01_raw"
    / "ssb"
    / "household statistics"
)

# %% [markdown]
# The processed Norwegian IO data.
# %%
norwegian_IO_path: Path = (
    config.data_path
    / "02_interim"
    / "norwegian IO"
    / str(config.survey_year)
)

# %% [markdown]
# The path to the household footprints (CaFEAN 1 results)
# %%
cafean_result_path: Path = (
    config.data_path
    / "03_results"
    / "CaFEAN"
)

# %% [markdown]
# The path to store the disaggregated footprint results
# %%
disaggregated_footprint_path: Path = (
    config.data_path
    / "02_interim"
    / "household survey footprints"
    / str(config.survey_year)
)
disaggregated_footprint_path.mkdir(exist_ok=True, parents=True)
# %% [markdown]
# # Load and prepare data
# %% [markdown]
# Load mapping between survey and IO classification and prepare it.  
# %%
survey_io_concordance = (
    pd.read_excel(
        (
            config.data_path
            / "00_auxiliary"
            / "concordances"
            / "survey_and_io.xlsx"
        ),
        index_col=[0, 1, 2, 3, 4]
    )
    .rename_axis(["sector"], axis=1)
)
survey_io_concordance

level_mappers = dict()
for x in [4, 3, 2]:
    for y in range(x-1):
        level_mappers[f"map_{x}_to_{y+1}"] = (
            survey_io_concordance
            .reset_index()
            .set_index([f"level {x}"])
            .loc[:, f"level {y+1}"]
            .to_dict()
        )

level_mappers["map_4_to_code"] = (
    survey_io_concordance
    .reset_index()
    .set_index(["level 4"])
    .loc[:, "code"]
    .to_dict()
)

concordance_matrix = (
    survey_io_concordance
    .droplevel([0, 1, 2, 4])
)


# %% [markdown]
# Load surveys. 
# Note: Survey 14157 (Expenditure per household per year, by commodity and service group, type of household, income quartile, contents and year)
# is not used further as matching population statistics is missing, 
# so it is not possible to scale up properly.
# %%
survey_14100 = load_survey(surveys_path, name="14100", n_index=2)
#survey_14152 = load_survey(surveys_path, name="14152", n_index=4)
#survey_14156 = load_survey(surveys_path, name="14156", n_index=4)
#survey_14157 = load_survey(surveys_path, name="14157", n_index=6)
#survey_14161 = load_survey(surveys_path, name="14161", n_index=4)

# %%
survey_names_to_code = (
    survey_14100
    .reset_index()
    .set_index("sector")["code"]
    .to_dict()
)
# %%
# %% [markdown]
# Load general statistics used to scale up base survey (14100).

# %%
household_statistics = (
    pd.read_excel(
        household_statistics_paths / "10986.xlsx",
        index_col=[0, 1, 2, 3, 4, 5],
        skiprows=2,
    )
    .dropna(how="all", axis=0)
    .rename_axis([
        "variable_NO",
        "variable",
        "region_code",
        "region",
        "category_code",
        "category"
    ], axis=0)
    .rename_axis(["year"], axis=1)
    .loc[:, str(config.survey_year)]
    .droplevel([
        "variable_NO",
        "region_code",
        "region",
        "category_code",
    ])
    .unstack("variable")
)
household_counts = household_statistics.loc[:, "Private households"]
household_persons = household_statistics.loc[:, "Persons in private households"]

# %%
Y_domestic = pd.read_excel(
    norwegian_IO_path / "Y_domestic.xlsx",
    index_col=[0]
).fillna(0)

Y_import = pd.read_excel(
    norwegian_IO_path / "Y_import.xlsx",
    index_col=[0]
).fillna(0)

valuation_layers = pd.read_excel(
    norwegian_IO_path / "valuation_layers.xlsx",
    index_col=[0]
).fillna(0)

# %%
# Get household direct emissions 
household_direct_emissions = (
    pd.read_csv(
        (
            cafean_result_path
            / "household_emissions.tsv"
        ),
        sep="\t",
        index_col=[0, 1, 2]
    )
)
household_direct_emissions = (
    household_direct_emissions.loc[
        ~household_direct_emissions.index.isin(["Households, totals"], level="hhld_component")
    ]
    .loc[config.survey_year]
    .loc["GHG"]  # Note: Emission type can be changed here!
    .loc[:, "value"]
    .rename(lambda x: f"Households; {x}", axis=0)
)

# %% [markdown]
# Get household footprints calculated by the Hybrid SNAC model.
# %%
# _bp = basic prices, _pp = purchaser prices
footprint_bp = (
    pd.read_csv(
        cafean_result_path 
        / "full_household_footprint.tsv",
        sep="\t",
    )
)
footprint_bp = (
    footprint_bp[
        (footprint_bp["year"] == config.survey_year)
        & (footprint_bp["emission"] == "GHG")  # Note: Emission type can be changed here!
        & ~(footprint_bp["sector"] == "Final consumption expenditure by households")
    ].groupby(["sector"])["value"].sum()
)

# %%
# * Get margin multipliers
margin_share_multiplier = (
    valuation_layers.loc[:, "Trade and transport margins"]
    .div(
        valuation_layers.loc[:, "Total supply at basic prices"]
    )
)
taxes_share_multiplier = (
    valuation_layers.loc[:, "Taxes less subsidies on products"]
    .div(
        valuation_layers.loc[:, "Total supply at basic prices"]
    )
)

# %% [markdown]
# Re-allocate footprint from basic prices to purchaser prices
# %%
# Note: "Y" refers to only final consumption by households.
# * Calculate margin layers and purchaser price household consumption
Y_bp = (
    Y_domestic["Final consumption expenditure by households"]
    + Y_import["Final consumption expenditure by households"]
)
Y_margin = margin_share_multiplier.mul(Y_bp)
Y_taxes = taxes_share_multiplier.mul(Y_bp)
Y_pp = Y_bp + Y_margin + Y_taxes

# %%
# * Identify margin and non-margin products
margin_products = Y_margin[Y_margin.fillna(0) < 0].index
non_margin_products = Y_margin[Y_margin.fillna(0) >= 0].index

# %%
# * Calculate share of household consumption that needs to be reallocated due to margins
Y_margin_reallocate = (
    -Y_margin.loc[margin_products]
    .div(Y_bp.loc[margin_products])
)

# * Calculate footprint associated with margin that needs to be reallocated
footprint_reallocate_margin = (
    footprint_bp.loc[margin_products]
    .mul(Y_margin_reallocate)
)

# * Calculate the the non-margin products share of non-margin total
Y_non_margin_share_of_margins = (
    Y_margin.loc[non_margin_products]
    .div(Y_margin.loc[non_margin_products].sum())
)

# * Allocate margin footprint using non-margin products share of non-margin total
footprint_reallocate_non_margin = (
    Y_non_margin_share_of_margins
    .mul(footprint_reallocate_margin.sum())
)

# * Concatenate the two vectors to get the total reallocated footprint
footprint_reallocate = pd.concat(
    [-footprint_reallocate_margin, footprint_reallocate_non_margin],
    axis=0
)

# * Calculate footprint in purchaser prices
footprint_pp = (footprint_bp + footprint_reallocate).fillna(0)

# %%
# * Combine with direct emissions for later use
footprint_io = (
    pd.concat([
        footprint_pp, 
        household_direct_emissions
    ])
)

# %%
# %% [markdown]
# Use GRAS to balance most detailed survey (14100) with the household consumption in purchaser prices.
# %%
# * Use survey as row totals in GRAS. 
row_total = (
    survey_14100
    .loc[survey_14100.index.get_level_values("code_level") == 4, :]
    .droplevel([0, 2, 3, 4, 5, 6])
    .loc[concordance_matrix.index, "Expenditure (NOK)"]
    .replace("..", 0)
    .mul(household_counts.sum())  # Transform from per household to national total
    .div(1e6)  # Transform to Million NOK to match IO units
)

# * Use household consumption in purchaser price as column totals for the balancing.
column_total = Y_pp.fillna(0)

# * Rescale household survey total to household consumption total (so it matches the total footprint)
row_total_rescaled = (
    row_total
    .mul(
        column_total.sum() / row_total.sum()
    )
)

# * Use concordance matrix as starting point for balancing.
post_balancing = gras(
    concordance_matrix.loc[row_total_rescaled.index, column_total.index].values,
    coltot=column_total.values,
    rowtot=row_total_rescaled.values,
    iter_in=20
)
print(f"Table total:", post_balancing.sum().sum())
print("IO total:", column_total.sum())
print("Survey total:", row_total_rescaled.sum())

# Turn the post balancing results into a pandas DataFrame
foootprint_allocation_matrix = pd.DataFrame(
    post_balancing,
    index=row_total_rescaled.index,
    columns=column_total.index
).fillna(0)


# %%
rescale_factor = (row_total_rescaled / row_total).fillna(1)
rescale_difference = (row_total_rescaled - row_total)
post_balancing_rows = foootprint_allocation_matrix.sum(axis=1)
post_balancing_factor = (post_balancing_rows / row_total_rescaled).fillna(1)
post_balancing_difference = (post_balancing_rows - row_total_rescaled)

survey_adjustments = pd.concat(
    [
        row_total,
        row_total_rescaled,
        rescale_factor,
        rescale_difference,
        post_balancing_rows,
        post_balancing_factor,
        post_balancing_difference,
    ],
    keys=[
        "Original",
        "Rescaled",
        "Rescaled / Original",
        "Rescaled - Original",
        "Post balancing",
        "Post balancing / Rescaled",
        "Post balancing - Rescaled",
    ],
    axis=1
)
# %%
# Allocate direct emissions according to survey expenditure values
direct_emission_allocation_matrix = (
    concordance_matrix
    .loc[:, household_direct_emissions.index]
    .mul(row_total, axis=0)
)

# %%
# * Combine allocation matrices
allocation_matrix = (
    pd.concat([
        foootprint_allocation_matrix,
        direct_emission_allocation_matrix
    ], axis=1)
)

# * Normalise each column
allocation_matrix = (
    allocation_matrix
    .div(allocation_matrix.sum(axis=0), axis=1)
    .replace([np.inf, -np.inf, np.nan], 0)
)

# %%
# Disaggregate footprint to household survey classification
disaggregated_footprint = allocation_matrix.dot(footprint_io)

# And aggregate up to intermediate levels.
disaggregated_footprint_level3 = (
    disaggregated_footprint
    .rename(level_mappers["map_4_to_3"])
    .groupby(level=0)
    .sum()
    .rename_axis(["level 3"])
)

disaggregated_footprint_level2 = (
    disaggregated_footprint
    .rename(level_mappers["map_4_to_2"])
    .groupby(level=0)
    .sum()
    .rename_axis(["level 2"])
)

disaggregated_footprint_level1 = (
    disaggregated_footprint
    .rename(level_mappers["map_4_to_1"])
    .groupby(level=0)
    .sum()
    .rename_axis(["level 1"])
)

# %%
# Add indexes to the allocation matrix
allocation_matrix_levels = (
    allocation_matrix
    .reset_index()
)

allocation_matrix_levels["level 3"] = (
    allocation_matrix_levels["level 4"]
    .map(level_mappers["map_4_to_3"])
)

allocation_matrix_levels["level 2"] = (
    allocation_matrix_levels["level 4"]
    .map(level_mappers["map_4_to_2"])
)

allocation_matrix_levels["level 1"] = (
    allocation_matrix_levels["level 4"]
    .map(level_mappers["map_4_to_1"])
)

allocation_matrix_levels["code"] = (
    allocation_matrix_levels["level 4"]
    .map(level_mappers["map_4_to_code"])
)

allocation_matrix_levels = (
    allocation_matrix_levels
    .set_index(["level 1", "level 2", "level 3", "level 4", "code"])
    .sort_index()
)

# %%
# Add indexes to the post balancing expenditures
post_balancing_expenditures = (
    foootprint_allocation_matrix.sum(axis=1)
    .to_frame("Expenditure (M.NOK.)")
    .reset_index()
)
post_balancing_expenditures["level 3"] = (
    post_balancing_expenditures["level 4"]
    .map(level_mappers["map_4_to_3"])
) 
post_balancing_expenditures["level 2"] = (
    post_balancing_expenditures["level 4"]
    .map(level_mappers["map_4_to_2"])
) 
post_balancing_expenditures["level 1"] = (
    post_balancing_expenditures["level 4"]
    .map(level_mappers["map_4_to_1"])
) 
post_balancing_expenditures["code"] = (
    post_balancing_expenditures["level 4"]
    .map(level_mappers["map_4_to_code"])
) 

post_balancing_expenditures = (
    post_balancing_expenditures
    .set_index(["level 1", "level 2", "level 3", "level 4", "code"])
)

# %%
# Make detailed footprint dataset for further analysis
detailed_footprint = (
    allocation_matrix_levels
    .dot(footprint_io)
    .to_frame("Footprint (Mt CO2eq)")
)
detailed_footprint = (
    detailed_footprint.merge(
        post_balancing_expenditures,
        right_index=True,
        left_index=True,
    )
)

detailed_footprint_level_4 = (
    detailed_footprint
    .groupby(["level 4"]).sum()
    .rename_axis(["Category"])
) 
detailed_footprint_level_4["Emission intensity (t CO2eq/NOK)"] = (
    detailed_footprint_level_4["Footprint (Mt CO2eq)"]
    .div(detailed_footprint_level_4["Expenditure (M.NOK.)"])
    .fillna(0)
)

detailed_footprint_level_3 = (
    detailed_footprint
    .groupby(["level 3"]).sum()
    .rename_axis(["Category"])
) 
detailed_footprint_level_3["Emission intensity (t CO2eq/NOK)"] = (
    detailed_footprint_level_3["Footprint (Mt CO2eq)"]
    .div(detailed_footprint_level_3["Expenditure (M.NOK.)"])
    .fillna(0)
)

detailed_footprint_level_2 = (
    detailed_footprint
    .groupby(["level 2"]).sum()
    .rename_axis(["Category"])
) 
detailed_footprint_level_2["Emission intensity (t CO2eq/NOK)"] = (
    detailed_footprint_level_2["Footprint (Mt CO2eq)"]
    .div(detailed_footprint_level_2["Expenditure (M.NOK.)"])
    .fillna(0)
)

detailed_footprint_level_1 = (
    detailed_footprint
    .groupby(["level 1"]).sum()
    .rename_axis(["Category"])
) 
detailed_footprint_level_1["Emission intensity (t CO2eq/NOK)"] = (
    detailed_footprint_level_1["Footprint (Mt CO2eq)"]
    .div(detailed_footprint_level_1["Expenditure (M.NOK.)"])
    .fillna(0)
)

detailed_footprint_full = (
    pd.concat(
        [
            detailed_footprint_level_1,
            detailed_footprint_level_2,
            detailed_footprint_level_3, 
            detailed_footprint_level_4
        ],
        keys=[1, 2, 3, 4],
        names=["Level"]
    )
    .reset_index()
) 

detailed_footprint_full["Code"] = (
    detailed_footprint_full["Category"]
    .map(survey_names_to_code)
)
detailed_footprint_full = (
    detailed_footprint_full
    .set_index(["Category", "Level", "Code"])
    .sort_index(level="Code")
)

detailed_footprint_allocated = (
    allocation_matrix_levels
    .mul(footprint_io, axis=1)
    .loc[
        survey_io_concordance.index,
        survey_io_concordance.columns
    ]
)

# %%
# Calculate unallocated footprint

unallocated_footprint = (
    footprint_io
    .subtract(detailed_footprint_allocated.sum(axis=0))
    .sort_values(ascending=False)
)

# Save results
# %%
# Save results
disaggregated_footprint.to_excel(
    disaggregated_footprint_path / "footprint_level_4.xlsx"
)

disaggregated_footprint_level3.to_excel(
    disaggregated_footprint_path / "footprint_level_3.xlsx"
)

disaggregated_footprint_level2.to_excel(
    disaggregated_footprint_path / "footprint_level_2.xlsx"
)

disaggregated_footprint_level1.to_excel(
    disaggregated_footprint_path / "footprint_level_1.xlsx"
)

detailed_footprint_full.to_excel(
    disaggregated_footprint_path / "detailed_footprint.xlsx"
)

detailed_footprint_allocated.to_excel(
    disaggregated_footprint_path / "detailed_footprint_allocated.xlsx"
)

unallocated_footprint.to_excel(
    disaggregated_footprint_path / "unallocated_footprint.xlsx"
)

# Save allocation matrix for future reference
allocation_matrix_levels.to_excel(
    disaggregated_footprint_path
    / "allocation matrix.xlsx"
)

# Save adjustment made to survey
survey_adjustments.to_excel(
    disaggregated_footprint_path
    / "survey adjustments.xlsx"
)

