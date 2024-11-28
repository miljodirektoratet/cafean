# -*- coding: utf-8 -*-
# %% [markdown]
#  # Deflator preparation script
#
# This prepares deflators for domestic production and imports from two SSB datasets.

# %% [markdown]
# Copyright (C) 2024 XIO Sustainability Analytics, Inc
#
# Written by
#
# - Kajwan Rasul
#
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
import numpy as np
from functions.general import aggregate_chemical_sectors

# %% [markdown]
# Python internal packages (this don't need to be installed, part of standard python)
# %%
from pathlib import Path


# %% [markdown]
# ### Locations / folder definitions
#
#
# deflators_path: Path = (
#     config.data_path
#     / "01_raw"
#     / "SSB"
#     / "deflators"
# )

# %% [markdown]
# We specify the folder in which we save 
# the deflators
# %%
deflators_interim_path: Path = (
    config.data_path
    / "02_interim" 
    / "deflators"
)
deflators_interim_path.mkdir(exist_ok=True, parents=True)


# %% [markdown]
# # Script

# %% [markdown]
# # Load data

# %%
# Load mappings between deflators and IO data
# * Read deflator mappings
dom_mapping = (
    pd.read_excel(
        (
            config.data_path
            / "00_auxiliary" 
            / "concordances" 
            / "deflators_and_io.xlsx"
        ),
        index_col=[0],
        header=[0],
        sheet_name="domestic"
    )
    .rename_axis(["sector"], axis=1)
    .stack()
)
dom_mapping = (
    dom_mapping[dom_mapping != 0]
    .reset_index()
    .loc[:, ["category", "sector"]]
)

imp_mapping = (
    pd.read_excel(
        (
            config.data_path
            / "00_auxiliary" 
            / "concordances" 
            / "deflators_and_io.xlsx"
        ),
        index_col=[0],
        header=[0],
        sheet_name="import"
    )
    .rename_axis(["sector"], axis=1)
    .rename_axis(["product"], axis=0)
    .stack()
)
imp_mapping = (
    imp_mapping[imp_mapping != 0]
    .reset_index()
    .loc[:, ["product", "sector"]]
)

# %% [markdown]
# ## Load and prepare deflators
# %%
## Load deflators for domestic IO
# Read data
deflator_data_dom = (
    pd.read_excel(
        deflators_path / "09170_A64.xlsx",
        index_col=[0, 1, 2, 3],
        header=[0, 1, 2],
    ).droplevel([0, 1], axis=1)
    .dropna()
    .rename_axis(["indicator_code", "indicator_full", "category_code", "category"], axis=0)
    .rename_axis(["year"], axis=1)
    .droplevel(["indicator_code", "category_code"], axis=0)
    .stack()
    .to_frame("value")
    .reset_index()
)

# Need to split up an index and then further process data
deflator_data_dom = (
    pd.concat([
        deflator_data_dom, 
        (
            deflator_data_dom["indicator_full"]
            .str.split(". ", regex=False, expand=True)
            .rename({0: "indicator", 1: "variable"}, axis=1)
        )
    ], axis=1)
    .set_index(["variable", "indicator", "category", "year"])
    .rename({"Constant 2015 price (NOK million)": "Constant 2015 prices (NOK million)"})
    .loc[["Current prices (NOK million)", "Constant 2015 prices (NOK million)"], "value"]
    .drop(["2022", "2023"], level="year")
    .to_frame("values")
    .reset_index()
    .merge(  # Map to IO classification
        dom_mapping,
        on="category",
        how="outer"
    )
    .groupby(["variable", "indicator", "sector", "year"]).sum()
    .loc[:, "values"]
)

# Aggregate chemical sectors
deflator_data_dom = (
    aggregate_chemical_sectors(
        deflator_data_dom.unstack("sector"),
        axis=1,
        level=0
    ).stack()
)

# Calculate deflators
domestic_deflator_all = (
    deflator_data_dom
    .loc["Constant 2015 prices (NOK million)"]
    .div(
        deflator_data_dom
        .loc["Current prices (NOK million)"]
        .replace(0, np.nan)
    )
    .replace([np.nan, np.inf, -np.inf], 0)
    .reorder_levels([0, 1, 2])
    .sort_index()
)

# Extract the different domestic deflators
domestic_deflator = (
    domestic_deflator_all.loc["Output at basic values"]
)

domestic_intermediate_consumption_deflator = (
    domestic_deflator_all.loc["Intermediate consumption"]
)

value_added_deflator = (
    domestic_deflator_all.loc["Value added at basic prices"]
)


# %%
## Load deflators for import IO
# TODO: Optimise order of commands.
# Read data
deflator_data_imp = (
    pd.read_excel(
        deflators_path / "07337.xlsx",
        index_col=[0, 1, 2, 3],
        header=[0, 1, 2],
    )
    .droplevel([0, 1], axis=1)
    .dropna(how="all", axis=0)
    .rename_axis(["variable_code", "variable", "product_code", "product"], axis=0)
    .rename_axis(["year"], axis=1)
    .droplevel(["variable_code", "product_code"], axis=0)
    .stack()
    .loc[["Current prices (NOK million)", "Constant 2015 prices (NOK million)"]]
)

# Map to IO classification
deflator_data_imp = (
    deflator_data_imp
    .to_frame("values")
    .reset_index()
    .merge(
        imp_mapping,
        on="product",
        how="right"
    )
    .set_index(["variable", "sector", "year"])
    .loc[:, "values"]
)

# Aggregate chemical sectors
deflator_data_imp = (
    aggregate_chemical_sectors(
        deflator_data_imp.unstack("sector"),
        axis=1,
        level=0
    ).stack()
)

# Calculate deflator
import_deflator = (
    deflator_data_imp
    .loc["Constant 2015 prices (NOK million)"]
    .div(
        deflator_data_imp
        .loc["Current prices (NOK million)"]
        .replace(0, np.nan)
    )
    .replace([np.nan, np.inf, -np.inf], 0)
    .sort_index()
)

# %% [markdown]
# Save deflators
# domestic_deflator.to_excel(
#     deflators_interim_path / "domestic.xlsx"
# )
# domestic_intermediate_consumption_deflator.to_excel(
#     deflators_interim_path / "intermediate consumption.xlsx"
# )
# value_added_deflator.to_excel(
#     deflators_interim_path / "value added.xlsx"
# )
# import_deflator.to_excel(
#     deflators_interim_path / "import.xlsx"
# )

# %%
