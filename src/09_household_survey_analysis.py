 # -*- coding: utf-8 -*-
# %% [markdown]
#  # Household survey footprint
#
# This script disaggregates the household footprint calculated by the SNAC
# model using SSB's household consumption surveys.
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
from functions.general import (
    load_survey,
    survey_aggregate_table
)
import matplotlib.pyplot as plt

import seaborn as sns
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
# The path to the disaggregated footprint results
# %%
disaggregated_footprint_path: Path = (
    config.data_path
    / "02_interim"
    / "household survey footprints"
    / str(config.survey_year)
)

# %% [markdown]
# Paths to store results and plots 
# %%
survey_footprint_path: Path = (
    config.data_path
    / "03_results"
    / "Household survey footprints"
    / str(config.survey_year)
)
survey_footprint_path.mkdir(exist_ok=True, parents=True)

survey_footprint_plots: Path = (
    config.data_path
    / "04_plots"
    / "Household survey footprints"
    / str(config.survey_year)
)
survey_footprint_plots.mkdir(exist_ok=True, parents=True)

# %% [markdown]
# The names and classifications doesn't always match across SSB statistics. 
# Here we use a few dictionaries to correct for that.
# %%
household_type_map = (
    pd.read_excel(
        (
            config.data_path
            / "00_auxiliary"
            / "classifications"
            / "survey_categories.xlsx"
        ),
        sheet_name="Household type"
    )
    .set_index(["Categories"])
    .loc[:, "Harmonized"]
    .to_dict()
)

centrality_map = (
    pd.read_excel(
        (
            config.data_path
            / "00_auxiliary"
            / "classifications"
            / "survey_categories.xlsx"
        ),
        sheet_name="Centrality"
    )
    .set_index(["Categories"])
    .loc[:, "Harmonized"]
    .to_dict()
)
# %%
# Orders for showing the results
household_type_order = [
    "Living alone",
    "Mother or father with children, youngest child 0-19 years",
    "Couples with no children",
    "Couples with children, youngest child 0-6 years",
    "Couples with children, youngest child 7-19 years",
    "Other households"
]
centrality_order = [
    'Most central municipalities (index 925-1000)',
    'Second most central municipalities (index 870-924)',
    'Above medium central municipalities (index 775-869)',
    'Medium central municipalities (index 670-774)',
    'Least and second least central municipalities (index 0-669)'
]
income_quartile_order = [
    "Highest income quartile",
    "Third income quartile",
    "Second income quartile",
    "Lowest income quartile"
]


# %%
# Load mapping between survey and IO classification and prepare it.  
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
# %%
# Load disaggregated footprint
disaggregated_footprint_level4 = pd.read_excel(
    disaggregated_footprint_path / "footprint_level_4.xlsx",
    index_col=[0]
).loc[:, 0]

disaggregated_footprint_level3 = pd.read_excel(
    disaggregated_footprint_path / "footprint_level_3.xlsx",
    index_col=[0]
).loc[:, 0]

disaggregated_footprint_level2 = pd.read_excel(
    disaggregated_footprint_path / "footprint_level_2.xlsx",
    index_col=[0]
).loc[:, 0]


# %% [markdown]
# Load surveys. 
# Note: Survey 14157 (Expenditure per household per year, by commodity and service group, type of household, income quartile, contents and year)
# is not used further as matching population statistics is missing, 
# so it is not possible to scale up properly.
# Survey 14100 is only used for disaggregation in the previous step.
# %%
# survey_14100 = load_survey(surveys_path, name="14100", n_index=2)
survey_14152 = load_survey(surveys_path, name="14152", n_index=4)
survey_14156 = load_survey(surveys_path, name="14156", n_index=4)
# survey_14157 = load_survey(surveys_path, name="14157", n_index=6)
survey_14161 = load_survey(surveys_path, name="14161", n_index=4)

# %% [markdown]
# Load general household statistics and make aggregations
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

household_counts_agg = (
    household_counts
    .rename(household_type_map)
    .groupby(level=0).sum()
)
household_persons_agg = (
    household_persons
    .rename(household_type_map)
    .groupby(level=0).sum()
)

# %%
population_centrality = (
    pd.read_excel(
        household_statistics_paths / "07459.xlsx",
        skiprows=2,
        index_col=list(range(8))
    )
    .dropna(how="all", axis=0)
    .rename_axis([
        "centrality_code",
        "centrality",
        "gender_code",
        "gender",
        "age_code",
        "age",
        "year",
        "year_duplicate",
    ], axis=0)
    .droplevel([
        "centrality_code",
        "gender_code",
        "age_code",
        "year_duplicate",
    ], axis=0)
    .loc[:, "Persons"]
    .unstack("year")
    .loc[:, config.survey_year]
    .groupby("centrality").sum()
)

household_centrality = (
    pd.read_excel(
        household_statistics_paths / "06070.xlsx",
        skiprows=2,
        index_col=list(range(6))
    )
    .dropna(how="all", axis=0)
    .rename_axis([
        "variable_NO",
        "variable",
        "centrality_code",
        "centrality",
        "household_type_code",
        "household_type"
    ], axis=0)
    .droplevel([
        "variable_NO",
        "variable",
        "centrality_code",
        "household_type_code",
    ], axis=0)
    .rename_axis(["year"], axis=1)
    .loc[:, str(config.survey_year)]
    .groupby(["centrality"]).sum()
)

population_centrality_agg = (
    population_centrality
    .rename(centrality_map)
    .groupby("centrality").sum()
)

household_centrality_agg = (
    household_centrality
    .rename(centrality_map)
    .groupby("centrality").sum()
)

# %% [markdown]
# # Survey 14152: Expenditure by type of household analysis
# %%
survey_14152_prepared = (
    survey_14152
    .droplevel("category_code")
    .unstack(["category"])
    .loc[3, "Expenditure (NOK)"]
    .drop("All households", axis=1)
)

spending_by_householdtype = (
    survey_14152_prepared.mul(household_counts_agg, axis=1)
    .droplevel([1, 2, 3, 4])
)

agg_spending_by_householdtype = (
    spending_by_householdtype.sum(axis=0)
)

total_spending_and_footprint = (
    pd.concat(
        [spending_by_householdtype.sum(axis=1), disaggregated_footprint_level3],
        keys=["allocation", "footprint"],
         axis=1
    )
)
survey_non_allocated = total_spending_and_footprint[
    (total_spending_and_footprint["allocation"] == 0)
    & (total_spending_and_footprint["footprint"] != 0)
]

for index in survey_non_allocated.index:
    spending_by_householdtype.loc[index, :] = agg_spending_by_householdtype

spending_by_householdtype_normalised = (
    spending_by_householdtype
    .div(
        spending_by_householdtype.sum(axis=1),
        axis=0
    )
    .replace([np.inf, -np.inf, np.nan], 0)
)

footprint_by_householdtype = (
    spending_by_householdtype_normalised
    .mul(disaggregated_footprint_level3, axis=0)
    .loc[:, household_type_order]
)

footprint_by_householdtype_per_household = (
    footprint_by_householdtype
    .div(household_counts_agg, axis=1)
    .mul(1e6)  # Mt to tons? 
    .loc[:, household_type_order]
)

footprint_by_householdtype_per_person = (
    footprint_by_householdtype
    .div(household_persons_agg, axis=1)
    .mul(1e6)  # Mt to tons?
    .loc[:, household_type_order]
)

householdtype_footprint = pd.concat(
    [
        footprint_by_householdtype,
        footprint_by_householdtype_per_household,
        footprint_by_householdtype_per_person
    ],
    keys=[
        "Total",
        "Per household",
        "Per person"
    ],
    names=["variable"],
    axis=1
)
householdtype_footprint.name = "Footprint by household type"

# %% [markdown]
# # Survey 14156: Expenditure by income quartile
# %%
survey_14156_prepared = (
    survey_14156
    .drop(["0", "0EU"], level="category_code", axis=0)
    .droplevel([
        "category_code",
    ])
    .loc[2, "Expenditure (NOK)"]
    .unstack("category")
)
spending_by_income_quartile = (
    survey_14156_prepared.mul(household_counts_agg.sum()/4)
    .droplevel([1, 2, 3])
)

agg_spending_by_income_quartile = (
    spending_by_income_quartile.sum(axis=0)
)

total_spending_and_footprint = (
    pd.concat(
        [spending_by_income_quartile.sum(axis=1), disaggregated_footprint_level2],
        keys=["allocation", "footprint"],
         axis=1
    )
)
survey_non_allocated = total_spending_and_footprint[
    (total_spending_and_footprint["allocation"] == 0)
    & (total_spending_and_footprint["footprint"] != 0)
]

for index in survey_non_allocated.index:
    spending_by_income_quartile.loc[index, :] = agg_spending_by_income_quartile

spending_by_income_quartile_normalised = (
    spending_by_income_quartile
    .div(
        spending_by_income_quartile.sum(axis=1),
        axis=0
    )
    .replace([np.inf, -np.inf, np.nan], 0)
)

footprint_by_income_quartile = (
    spending_by_income_quartile_normalised
    .mul(disaggregated_footprint_level2, axis=0)
    .loc[:, income_quartile_order]
)

footprint_by_income_quartile_per_household = (
    footprint_by_income_quartile
    .div(household_counts_agg.sum()/4, axis=1)
    .mul(1e6)  # Mt to tons? 
    .loc[:, income_quartile_order]
)

income_quartile_footprint = (
    pd.concat(
        [
            footprint_by_income_quartile,
            footprint_by_income_quartile_per_household,
        ],
        keys=[
            "Total",
            "Per household"
        ],
        names=["variable"],
        axis=1
    )
    .rename_axis(["sector"], axis=0)
)
income_quartile_footprint.name = "Footprint by income quartile"

# %% [markdown]
# # Survey 14156: Expenditure by centrality

# %%
survey_14161_prepared = (
    survey_14161
    .droplevel("category_code")
    .unstack(["category"])
    .loc[2, "Expenditure (NOK)"]
    .drop("All municipalities", axis=1)
)

spending_by_centrality = (
    survey_14161_prepared.mul(household_centrality_agg, axis=1)
    .droplevel([1, 2, 3])
)

agg_spending_by_centrality = (
    spending_by_centrality.sum(axis=0)
)

total_spending_and_footprint = (
    pd.concat(
        [spending_by_centrality.sum(axis=1), disaggregated_footprint_level2],
        keys=["allocation", "footprint"],
         axis=1
    )
)
survey_non_allocated = total_spending_and_footprint[
    (total_spending_and_footprint["allocation"] == 0)
    & (total_spending_and_footprint["footprint"] != 0)
]

for index in survey_non_allocated.index:
    spending_by_centrality.loc[index, :] = agg_spending_by_centrality

spending_by_centrality_normalised = (
    spending_by_centrality
    .div(
        spending_by_centrality.sum(axis=1),
        axis=0
    )
    .replace([np.inf, -np.inf, np.nan], 0)
)

footprint_by_centrality = (
    spending_by_centrality_normalised
    .mul(disaggregated_footprint_level2, axis=0)
    .loc[:, centrality_order]
)

footprint_by_centrality_per_household = (
    footprint_by_centrality
    .div(household_centrality_agg, axis=1)
    .mul(1e6)  # Mt to tons? 
    .loc[:, centrality_order]
)

footprint_by_centrality_per_person = (
    footprint_by_centrality
    .div(population_centrality_agg, axis=1)
    .mul(1e6)  # Mt to tons?
    .loc[:, centrality_order]
)

centrality_footprint = pd.concat(
    [
        footprint_by_centrality,
        footprint_by_centrality_per_household,
        footprint_by_centrality_per_person
    ],
    keys=[
        "Total",
        "Per household",
        "Per person"
    ],
    names=["variable"],
    axis=1
)
centrality_footprint.name = "Footprint by centrality"

# %%
# Add column totals
householdtype_footprint.loc["Total", :] = householdtype_footprint.sum(axis=0)
income_quartile_footprint.loc["Total", :] = income_quartile_footprint.sum(axis=0)
centrality_footprint.loc["Total", :] = centrality_footprint.sum(axis=0)

# %% [markdown]
# # Save results
# %%
# Then save survey footprint results
householdtype_footprint.to_excel(
    survey_footprint_path / "Footprint by household type.xlsx"
)
income_quartile_footprint.to_excel(
    survey_footprint_path / "Footprint by income quartile.xlsx"
)
centrality_footprint.to_excel(
    survey_footprint_path / "Footprint by centrality.xlsx"
)


# %% [markdown]
# # Create aggregates tables and save them
# %%
survey_aggregate_table(
    survey=householdtype_footprint,
    variable="Per person",
    mappers=level_mappers,
    column_order=household_type_order,
    show_table=True,
    save_table=True,
    save_path=survey_footprint_path
)
# %%
survey_aggregate_table(
    survey=householdtype_footprint,
    variable="Per household",
    mappers=level_mappers,
    show_table=True,
    save_table=True,
    save_path=survey_footprint_path
)
# %%
survey_aggregate_table(
    survey=income_quartile_footprint,
    variable="Per household",
    mappers=level_mappers,
    column_order=income_quartile_order,
    show_table=True,
    save_table=True,
    save_path=survey_footprint_path
)
# %%
survey_aggregate_table(
    survey=centrality_footprint,
    variable="Per person",
    mappers=level_mappers,
    column_order=centrality_order,
    show_table=True,
    save_table=True,
    save_path=survey_footprint_path
)
# %%
survey_aggregate_table(
    survey=centrality_footprint,
    variable="Per household",
    mappers=level_mappers,
    column_order=centrality_order,
    show_table=True,
    save_table=True,
    save_path=survey_footprint_path
)

# %% [markdown]
# Make plots for documentation
# %%
householdtype_footprint_plot = (
    householdtype_footprint
    .sum(axis=0)
    .unstack("variable")
    .sort_values(by="Per person")
    .round(2)
    .rename(lambda x: x.replace(", ", ",\n"))
)

fig, axes = plt.subplots(
    figsize=(10, 15),
    nrows=3,
    sharex=True,
    gridspec_kw={"hspace": 0.05}
)

sns.barplot(
    data=householdtype_footprint_plot,
    x="category",
    y="Per person",
    ax=axes[0],
    hue="category"
)
axes[0].set_ylabel("Per person [t]")
sns.barplot(
    data=householdtype_footprint_plot,
    x="category",
    y="Per household",
    ax=axes[1],
    hue="category"
)
axes[1].set_ylabel("Per household [t]")
sns.barplot(
    data=householdtype_footprint_plot,
    x="category",
    y="Total",
    ax=axes[2],
    hue="category"
)
axes[2].set_xticks(
    axes[2].get_xticks(),
    axes[2].get_xticklabels(),
    rotation=45,
    ha='right'
)
axes[2].set_ylabel("Total [Mt]")
axes[2].set_xlabel(None)
fig.supylabel("CO2 equivalent")
fig.savefig(
    survey_footprint_plots / "Footprint by household type.png",
    bbox_inches="tight"
)
fig.show()

# %%
income_quartile_footprint_plot = (
    income_quartile_footprint
    .sum(axis=0)
    .unstack("variable")
    #.sort_values(by="Per household")
    .round(2)
)

fig, axes = plt.subplots(
    figsize=(10, 10),
    nrows=2,
    sharex=True, 
    gridspec_kw={"hspace": 0.05}
)

sns.barplot(
    data=income_quartile_footprint_plot,
    x="category",
    y="Per household",
    ax=axes[0],
    hue="category"
)
axes[0].set_ylabel("Per household [t]")
sns.barplot(
    data=income_quartile_footprint_plot,
    x="category",
    y="Total",
    ax=axes[1],
    hue="category",
)
axes[1].set_ylabel("Total [Mt]")
axes[1].set_xticks(
    axes[1].get_xticks(),
    axes[1].get_xticklabels(),
    rotation=45,
    ha='right'
)
axes[1].set_xlabel(None)
fig.supylabel("CO2 equivalent")
fig.savefig(
    survey_footprint_plots / "Footprint by income quartile.png",
    bbox_inches="tight"
)
fig.show()

# %%
centrality_footprint_plot = (
    centrality_footprint
    .sum(axis=0)
    .unstack("variable")
    #.loc[centrality_order, :]
    .rename(lambda x: x.replace(" (", "\n("))
    .round(2)
)

fig, axes = plt.subplots(
    figsize=(10, 15),
    nrows=3,
    sharex=True,
    gridspec_kw={"hspace": 0.05}
)

sns.barplot(
    data=centrality_footprint_plot,
    x="category",
    y="Per person",
    ax=axes[0],
    hue="category"
)
axes[0].set_ylabel("Per person [t]")
sns.barplot(
    data=centrality_footprint_plot,
    x="category",
    y="Per household",
    ax=axes[1],
    hue="category"
)
axes[1].set_ylabel("Per household [t]")
sns.barplot(
    data=centrality_footprint_plot,
    x="category",
    y="Total",
    ax=axes[2],
    hue="category"
)
axes[2].set_xticks(
    axes[2].get_xticks(),
    axes[2].get_xticklabels(),
    rotation=45,
    ha='right'
)
axes[2].set_ylabel("Total [Mt]")
axes[2].set_xlabel(None)
fig.supylabel("CO2 equivalent")
fig.savefig(
    survey_footprint_plots / "Footprint by centrality.png",
    bbox_inches="tight"
)
fig.show()


# %% [markdown]
# Make plots for presentation
# %%
householdtype_footprint_plot = (
    householdtype_footprint
    .sum(axis=0)
    .unstack("variable")
    .sort_values(by="Per person")
    .round(2)
    .rename(lambda x: x.replace(", ", ",\n"))
)

fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=householdtype_footprint_plot,
    x="category",
    y="Per person",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Per person\nCO2 eqv. [t]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by household type (per person).png",
   bbox_inches="tight"
)
fig.show()

fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=householdtype_footprint_plot,
    x="category",
    y="Per household",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Per household\nCO2 eqv. [t]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by household type (per household).png",
   bbox_inches="tight"
)
fig.show()

fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=householdtype_footprint_plot,
    x="category",
    y="Total",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Total\nCO2 eqv. [Mt]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by household type (total).png",
   bbox_inches="tight"
)
fig.show()

# %%
income_quartile_footprint_plot = (
    income_quartile_footprint
    .sum(axis=0)
    .unstack("variable")
    #.sort_values(by="Per household")
    .round(2)
)


fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=income_quartile_footprint_plot,
    x="category",
    y="Per household",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Per household\nCO2 eqv. [t]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by income quartile (per household).png",
   bbox_inches="tight"
)
fig.show()

fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=income_quartile_footprint_plot,
    x="category",
    y="Total",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Total\nCO2 eqv. [Mt]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by income quartile (total).png",
   bbox_inches="tight"
)
fig.show()

# %%
centrality_footprint_plot = (
    centrality_footprint
    .sum(axis=0)
    .unstack("variable")
    #.loc[centrality_order, :]
    .rename(lambda x: x.replace(" (", "\n("))
    .round(2)
)

fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=centrality_footprint_plot,
    x="category",
    y="Per person",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Per person\nCO2 eqv. [t]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by centrality (per person).png",
   bbox_inches="tight"
)
fig.show()

fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=centrality_footprint_plot,
    x="category",
    y="Per household",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Per household\nCO2 eqv. [t]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by centrality (per household).png",
   bbox_inches="tight"
)
fig.show()

fig, ax = plt.subplots(
    figsize=(10, 3),
)
sns.barplot(
    data=centrality_footprint_plot,
    x="category",
    y="Total",
    ax=ax,
    hue="category"
)
ax.set_ylabel("Total\nCO2 eqv. [Mt]")
ax.set_xticks(
    ax.get_xticks(),
    ax.get_xticklabels(),
    rotation=45,
    ha='right'
)
ax.set_xlabel(None)
fig.savefig(
   survey_footprint_plots / "Footprint by centrality (total).png",
   bbox_inches="tight"
)
fig.show()
# %%
