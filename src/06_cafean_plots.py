# -*- coding: utf-8 -*-
# %% [markdown]
#  # Figures and tables for the CaFEAN report
#
# This script creates all the figures and tables used in the original CaFEAN report,
# both for the original (legacy) data, and the new updated report.
# All functions are kept in the script to make it easy to configure.

# %% [markdown]
# Copyright (C) 2023  XIO Sustainability Analytics, Inc
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
# Python internal packages (this don't need to be installed, part of standard python)
# %%
import os
from pathlib import Path

# %% [markdown]
# Import of required external packages. These are all available through pip/conda/mamba install.
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %% [markdown]
# ### Locations / folder definitions
# %%

cafean_result_path: Path = (
    config.data_path
    / "03_results"
    / "CaFEAN"
)

cafean_save_path: Path = (
    config.data_path
    / "04_plots"
    / "CaFEAN"
)

legacy_cafean_result_path: Path = (
    config.data_path
    / "legacy"
    / "CaFEAN results"
)

legacy_cafean_save_path: Path = (
    legacy_cafean_result_path
    / "figures"
)

# %%
# Charactarization factors
charactarization_factors_ghg = pd.read_excel(
    (
        config.data_path
        / "00_auxiliary" 
        / "other"
        / "emis_conv_char.xlsx"
    ), sheet_name="nor_char"
)

# %%
set(plt.style.available)
plt.style.use('seaborn-v0_8-whitegrid')


color_category = {
    "Production": "tab:olive",
    "Import": "tab:blue",
    "Export": "tab:orange",
    "Footprint": "tab:green",
    "Domestic": "tab:olive",
    "Imported": "tab:blue"
}

color_emission = {
    "CO2": "tab:blue",
    "CH4": "tab:orange",
    "N2O": "tab:green",
    "HFC": "tab:red",
    "SF6_NF3": "tab:purple",
    "PFC": "tab:brown"
}

color_source = {
    "Norway": "tab:green",
    "Europe": "tab:blue",
    "Developing": "tab:red",
    "Advanced": "tab:purple"
}

# %%
# Simple function to shorten names in plots.
def shorten_name(x, max_len=45):
    if len(x) > max_len:
        x = f"{x[:max_len]}..."
    return x


# %%
def figure5(
        result_path, 
        year=config.plot_year, 
        save_plot=True, 
        save_path=None,
        save_name="Figure5",
        save_extension="png",
        return_data=True
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / "GHG_totals.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        ).loc[
            [
                "production_account_with_households",
                "footprint_gross_imports",
                "footprint_gross_exports",
                "total_footprint_with_households"    
            ],
            f"{year}"
        ]
        .rename(
            {
                "production_account_with_households": "Production",
                "footprint_gross_imports": "Import",
                "footprint_gross_exports": "Export",
                "total_footprint_with_households": "Footprint"    
            }
        )
        .to_frame("values")
    )
    df.loc["Export"] = df.loc[["Export"]]*(-1)
    df.loc[:, "bottom"] = df.cumsum().shift(1)
    df = df.round(1).reset_index()
    df.loc[0, "bottom"] = 0
    df.loc[len(df)-1, "bottom"] = 0
    
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(
        x=df["component"],
        height=df["values"],
        bottom=df["bottom"],
        color=list(color_category.values()),
        width=0.9,
        edgecolor="k"
    )
    ax.set_ylabel("CO2-eq [Mt]")
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax*1.1)
    ax.bar_label(ax.containers[0], label_type="center")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_title(f"GHG emissions in {year}")
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")

    fig.show()

    if return_data:
        return df.set_index(["component", "emission", "unit"])["values"]
figure5(cafean_result_path, save_plot=False)

# %%
def figure6(
        result_path, 
        year=config.plot_year, 
        save_plot=True, 
        save_path=None,
        save_name="Figure6",
        save_extension="png",
        return_data=True
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / "GHG_totals.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        )
        .droplevel([1, 2])
        .loc[:, f"{year}"]
    )

    df_plot = pd.DataFrame(
        None,
        index=("Production", "Footprint", "Import"),
        columns=("Domestic", "Export", "Import")
    )

    df_plot.loc["Production", "Domestic"] = (
        df["domestic_footprint"]
        + df["household_totals"]
    )
    df_plot.loc["Footprint", "Domestic"] = (
        df["domestic_footprint"]
        + df["household_totals"]
    )
    df_plot.loc["Production", "Export"] = (
        df["production_account_with_households"]
        - df_plot.loc["Footprint", "Domestic"]
    )
    df_plot.loc["Import", "Export"] = (
        df["footprint_gross_imports"]
        - (
            df["total_footprint_with_households"]
            - df_plot.loc["Production", "Domestic"]
        )
    )
    df_plot.loc["Footprint", "Import"] = (
        df["import_footprint"]
    )

    df_plot = df_plot.round(1)
    
    fig, ax = plt.subplots(figsize=(8, 4))
    df_plot.plot.bar(
        stacked=True,
        color=color_category,
        ax=ax,
        width=0.9,
        edgecolor="k"
    )
    ax.set_ylabel("CO2-eq [Mt]")
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax*1.1)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_title(f"GHG emissions in {year}")
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")

    fig.show()

    if return_data:
        return df_plot
figure6(cafean_result_path, save_plot=False)

# %%
def figure7(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure7",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "total_accounts.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        )
        .loc[year, :]
        .loc[
            [
                "production_account_with_households",
                "footprint_gross_imports",
                "footprint_gross_exports",
                "total_footprint_with_households"    
            ],
            :
        ]
        .rename(
            {
                "production_account_with_households": "Production",
                "footprint_gross_imports": "Import",
                "footprint_gross_exports": "Export",
                "total_footprint_with_households": "Footprint"    
            }
        )
        .reset_index()
        .merge(
            charactarization_factors_ghg,
            how="left",
            left_on=["emission", "unit"],
            right_on=["stressor", "stressor_unit"]
        )
    )
    df["GHG eq"] = df["value"].mul(df["factor"])
    df = (
        df
        .set_index(["component", "emission"])
        .loc[:, "GHG eq"]
        .replace(0, np.nan)
        .dropna(how="all")
        .unstack("emission")
    )
    category_order = ["Production", "Import", "Export", "Footprint"]
    gas_order = df.sum(axis=0).sort_values(ascending=False).index
    df = df.loc[category_order, gas_order]
    df.loc["Export", :] = df.loc["Export", :]*(-1)

    fig, ax = plt.subplots(figsize=(8, 4))
    (
        df
        .plot.bar(
            stacked=True,
            ax=ax,
            color=color_emission,
            width=0.9,
            edgecolor="k"
        )
    )
    ax.set_title(f"GHG by gas type in {year}")
    ax.axhline(color="k")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_xlabel(None)
    ax.legend(title=None, ncols=2, loc="lower left")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df
figure7(cafean_result_path, save_plot=False)

# %%
# Figure 8 and 9
def figure8_and_9(
        result_path, 
        save_plot=True, 
        save_path=None,
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / "GHG_totals.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        ).loc[
            [
                "production_account_with_households",
                "footprint_gross_imports",
                "footprint_gross_exports",
                "total_footprint_with_households"    
            ],
            :
        ]
        .rename(
            {
                "production_account_with_households": "Production",
                "footprint_gross_imports": "Import",
                "footprint_gross_exports": "Export",
                "total_footprint_with_households": "Footprint"    
            }
        )
        .rename_axis(["year"], axis=1)
    )

    df_index = (
        df.div(df.loc[:, "2012"], axis=0)
        .stack()
        .to_frame("values")
        .reset_index()
    )

    df = (
        df
        .stack()
        .to_frame("values")
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.lineplot(data=df, x="year", y="values", hue="component", palette=color_category, ax=ax)
    ax.legend(title=None, ncols=4)
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_xlabel(None)
    ax.set_title("Trend in GHG emissions")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"Figure8.{save_extension}", bbox_inches="tight")
    fig.show()

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.lineplot(data=df_index, x="year", y="values", hue="component", palette=color_category, ax=ax)
    ax.legend(title=None, ncols=4)
    ax.set_title("Trend in GHG emissions (Index 2012)")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_xlabel(None)
    ax.set_ylabel(None) 
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"Figure9.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return (
            df.set_index(["year", "component", "emission", "unit"]),
            df_index.set_index(["year", "component", "emission", "unit"])
        )
figure8_and_9(cafean_result_path, save_plot=False)
# %%
def figure10(
        result_path,
        save_plot=True, 
        save_path=None,
        save_name="Figure10",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / "GHG_totals.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        ).loc[
            [
                "domestic_footprint",
                "household_totals",
                "import_footprint",    
            ],
            :
        ]
        .rename(
            {
                "domestic_footprint": "Domestic",
                "household_totals": "Domestic",
                "import_footprint": "Imported",
            }
        )
        .rename_axis(["year"], axis=1)
        .stack()
        
    )
    df = df.groupby(df.index.names).sum()
    df_plot = df.unstack("component").droplevel([0, 1])

    fig, ax = plt.subplots(figsize=(8, 4))
    df_plot.plot.bar(stacked=True, ax=ax, color=color_category, width=0.9, edgecolor="k")
    ax.legend(title=None, ncols=3)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_title("GHG footprint")
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_xlabel(None)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df
figure10(cafean_result_path, save_plot=False)
# %%

def figure11_to_13_and_27(
        result_path, 
        ghg_gas="CO2",
        save_plot=True, 
        save_path=None,
        save_name="Figure11",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "total_accounts.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        )
        .rename(
            {
                "production_account_with_households": "Production",
                "footprint_gross_imports": "Import",
                "footprint_gross_exports": "Export",
                "total_footprint_with_households": "Footprint",
                "domestic_footprint": "Domestic",
                "household_totals": "Domestic",
                "import_footprint": "Imported",
            }, level="component"
        )
        .reset_index()
    )
    if ghg_gas in ["CO2", "Biomass CO2"]:
        df["value"] = df["value"]*1e-3
        df["unit"] = "Mt"
    elif ghg_gas in ["CH4", "N2O"]:
        df["value"] = df["value"]*1e-3
        df["unit"] = "kt"

    df = df[df["emission"] == ghg_gas]
    unit = df["unit"].iloc[0]

    df_left = df[
        df["component"].isin(["Production", "Import", "Export", "Footprint"])
    ]
    df_right = (
        df[df["component"].isin(["Domestic", "Imported"])]
        .set_index(["year", "component"])["value"]
        .groupby(["year", "component"]).sum()
        .unstack("component")
    )

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
    sns.lineplot(data=df_left, x="year", y="value", hue="component", palette=color_category, ax=axes[0])
    axes[0].legend(title=None)
    ymin, ymax = axes[0].get_ylim()
    axes[0].set_ylim(0, ymax)
    axes[0].set_title(f"Trend for {ghg_gas}")
    axes[0].set_ylabel(f"{ghg_gas} [{unit}]")
    axes[0].set_xlabel(None)
    plt.setp(axes[0].spines.values(), color="k", linewidth=1)
    df_right.plot.bar(stacked=True, ax=axes[1], color=color_category, width=0.9, edgecolor="k")
    axes[1].legend(title=None, ncols=1)
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45)
    axes[1].set_title(f"{ghg_gas} footprint source")
    axes[1].set_ylabel(f"{ghg_gas} [{unit}]")
    axes[1].set_xlabel(None)
    plt.setp(axes[1].spines.values(), color="k", linewidth=1)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")

    fig.show()

    if return_data:
        return df.set_index(["emission", "component", "year", "unit"])["value"]

figure11_to_13_and_27(cafean_result_path, save_plot=False)
# %%
def figure14(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure14",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_{year}_production_account.tsv",
            sep="\t",
            index_col=[0]
        ).loc[:, ["GHG"]]
        .sort_values(by="GHG", ascending=False)
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(data=df, x="sector", height="GHG", width=0.9, edgecolor="k")
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_title("Production account, GHG emissions")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")

    fig.show()

    if return_data:
        return df.set_index("sector")
figure14(cafean_result_path, save_plot=False)
# %%
def figure15(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure15",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_{year}_production_account.tsv",
            sep="\t",
            index_col=[0]
        )
        .drop("GHG", axis=1)
        .rename_axis(["emission"], axis=1)
        .stack()
        .to_frame("values")
        .reset_index()
        .merge(
            charactarization_factors_ghg,
            how="left",
            left_on=["emission"],
            right_on=["stressor"]
        )
        #.loc[:, ["GHG"]]
        #.sort_values(by="GHG", ascending=False)
        #.reset_index()
    )
    df["GHG eq"] = df["values"].mul(df["factor"])

    df = (
        df
        .set_index(["sector", "emission"])["GHG eq"]
        .unstack("emission")
        .drop("Biomass CO2", axis=1)
    )
    sector_order = df.sum(axis=1).sort_values(ascending=False).index
    emission_order = df.sum(axis=0).sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(8, 4))
    (
        df
        .loc[sector_order, emission_order]
        .plot.bar(
            stacked=True,
            ax=ax,
            color=color_emission,
            width=0.9,
            edgecolor="k"
        )
    )
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_xlabel(None)
    ax.set_title("Production account, GHG emissions by gas")
    ax.legend(title=None)
    plt.setp(ax.spines.values(), color="k", linewidth=1)

    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")

    fig.show()

    if return_data:
        return df
figure15(cafean_result_path, save_plot=False) 
# %%
def figure16(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure16",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_{year}_total_footprint.tsv",
            sep="\t",
            index_col=[0]
        ).loc[:, ["GHG"]]
        .sort_values(by="GHG", ascending=False)
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(data=df, x="sector", height="GHG", width=0.9, edgecolor="k")
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_title("Footprint account by final goods and services")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")

    fig.show()

    if return_data:
        return df.set_index("sector")
figure16(cafean_result_path, save_plot=False) 

# %%
def figure17(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure17",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_{year}_total_footprint.tsv",
            sep="\t",
            index_col=[0]
        )
        .drop("GHG", axis=1)
        .rename_axis(["emission"], axis=1)
        .stack()
        .to_frame("values")
        .reset_index()
        .merge(
            charactarization_factors_ghg,
            how="left",
            left_on=["emission"],
            right_on=["stressor"]
        )
        #.loc[:, ["GHG"]]
        #.sort_values(by="GHG", ascending=False)
        #.reset_index()
    )
    df["GHG eq"] = df["values"].mul(df["factor"])

    df = (
        df
        .set_index(["sector", "emission"])["GHG eq"]
        .unstack("emission")
        .drop("Biomass CO2", axis=1)
    )
    sector_order = df.sum(axis=1).sort_values(ascending=False).index
    emission_order = df.sum(axis=0).sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(8, 4))
    (
        df
        .loc[sector_order, emission_order]
        .plot.bar(
            stacked=True,
            ax=ax,
            color=color_emission,
            width=0.9,
            edgecolor="k"
        )
    )
    ax.legend(title=None)
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_title("Footprint account by final goods and services, and by gas")
    ax.set_xlabel(None)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")

    fig.show()

    if return_data:
        return df
figure17(cafean_result_path, save_plot=False) 

# %%
def figure18(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure18",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.concat([
            pd.read_csv(
                result_path / "excel_report_format" / f"sector_accounts_agg_{year}_domestic_footprint.tsv",
                sep="\t",
                index_col=[0]
            ),
            pd.read_csv(
                result_path / "excel_report_format" / f"sector_accounts_agg_{year}_import_footprint.tsv",
                sep="\t",
                index_col=[0]
            ),
        ], axis=0, keys=["Domestic", "Imported"], names=["Category"])
        .loc[:, "GHG"]
        .unstack("Category")
        #.sort_values(by="GHG", ascending=False)
        #.reset_index()
    )
    sector_order = df.sum(axis=1).sort_values(ascending=False).index
    category_order = df.sum(axis=0).sort_values(ascending=False).index
    fig, ax = plt.subplots(figsize=(8, 4))
    (
        df
        .loc[sector_order, category_order]
        .plot.bar(
            stacked=True,
            color=color_category,
            ax=ax,
            width=0.9,
            edgecolor="k"
        )
    )
    ax.set_xlabel(None)
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_title("Domestic and imported GHG footprint account")
    ax.legend(title=None)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()
    if return_data:
        return df
figure18(cafean_result_path, save_plot=False) 

# %%
def figure19(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure19",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_{year}_footprint_gross_imports.tsv",
            sep="\t",
            index_col=[0]
        ).loc[:, ["GHG"]]
        .sort_values(by="GHG", ascending=False)
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(data=df, x="sector", height="GHG", width=0.9, edgecolor="k")
    ax.set_title("GHG emissions embodied in imports")
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df.set_index("sector")
figure19(cafean_result_path, save_plot=False) 
# %%
def figure20(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure20",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_{year}_footprint_gross_exports.tsv",
            sep="\t",
            index_col=[0]
        ).loc[:, ["GHG"]]
        .sort_values(by="GHG", ascending=False)
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set_title("GHG emissions embodied in exports")
    ax.set_ylabel("CO2-eq [Mt]")
    ax.bar(data=df, x="sector", height="GHG", width=0.9, edgecolor="k")
    plt.setp(ax.spines.values(), color="k", linewidth=1)    
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df.set_index("sector")
figure20(cafean_result_path, save_plot=False) 

# %%
def table5_to_7(
        result_path,
        ghg_gas = "CO2",
        year_start=2012,
        year_end=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Table5",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "sector_accounts.tsv",
            sep="\t",
            index_col=[1, 2, 0, 3]
        ).loc["total_footprint", "value"]
        .loc[ghg_gas, :]
        .loc[[year_start, year_end], :]
        .unstack("year")
        .diff(axis=1)
        .loc[:, [year_end]]
        .sort_values(by=year_end, ascending=False)
        .rename(shorten_name, axis=0)
    )
    if ghg_gas in ["CO2", "CH4", "N2O"]:
        df = df.div(1e3)
    df_full = df.copy()
    df_bottom = df.tail(5)
    df_bottom["color"] = "tab:green"
    df_top = df.head(5)
    df_top["color"] = "tab:red"

    df = pd.concat([df_top, df_bottom], axis=0).round(2)

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.barplot(
        data=df.reset_index(),
        y="sector",
        x=year_end,
        hue="color",
        palette=["tab:red", "tab:green"],
        ax=ax,
        legend=False
    )
    ax.set_title(f"Changes in {ghg_gas} from {year_start} to {year_end} - Top 5 negative and positive changes")
    ax.bar_label(ax.containers[0], label_type="edge")
    ax.bar_label(ax.containers[1], label_type="edge")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.axvline(0, color="k", lw=1)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin*1.1, xmax*1.1)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df_full
table5_to_7(cafean_result_path, save_plot=False)
table5_to_7(cafean_result_path, year_start=2012, year_end=config.plot_year, save_plot=True, save_name=f"2012 to {config.plot_year} - Table5", save_path=cafean_save_path)
table5_to_7(cafean_result_path, year_start=2019, year_end=config.plot_year, save_plot=True, save_name=f"2019 to {config.plot_year} - Table5", save_path=cafean_save_path)
table5_to_7(cafean_result_path, year_start=2020, year_end=config.plot_year, save_plot=True, save_name=f"2020 to {config.plot_year} - Table5", save_path=cafean_save_path)
# %%
def figure21(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure21",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"emis_sources_totals_with_households_{year}.tsv",
            sep="\t",
            index_col=[0, 1]
        )
        .drop("Biomass CO2", axis=0)
        .rename_axis(["source"], axis=1)
        .stack()
        .to_frame("value")
        .reset_index()
        .merge(
            charactarization_factors_ghg,
            how="left",
            left_on=["emission"],
            right_on=["stressor"]
        )
        .fillna(1)
    )
    df["GHG eq"] = df["value"].mul(df["factor"])
    df = (
        df
        .set_index(["emission", "source"])
        .loc[:, "GHG eq"]
        .unstack("source")
    )
    gas_order = df.sum(axis=1).sort_values(ascending=False).index
    source_order = df.sum(axis=0).sort_values(ascending=False).index

    df = df.loc[gas_order, source_order]

    fig, ax = plt.subplots(figsize=(8, 4))
    df.plot.bar(stacked=True, color=color_source, ax=ax, width=0.9, edgecolor="k")
    ax.set_title("Emission by origin of emission and emission type")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_xlabel(None)
    ax.legend(title=None)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df
t1 = figure21(cafean_result_path, save_plot=False, return_data=True)

# %%
def figure22(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure22",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"ghg_source_reg_agg_{year}.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        )
        .droplevel(["emission", "unit"])
    )
    sector_order = df.sum(axis=1).sort_values(ascending=False).index
    source_order = df.sum(axis=0).sort_values(ascending=False).index

    df = df.loc[sector_order, source_order]

    fig, ax = plt.subplots(figsize=(8, 4))
    df.plot.bar(stacked=True, color=color_source, ax=ax, width=0.9, edgecolor="k")
    ax.set_title("Emission in final product by origin of emission")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_xlabel(None)
    ax.legend(title=None)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df
figure22(cafean_result_path, save_plot=False)

# %%
def figure23(
        result_path,
        year=config.plot_year,
        save_plot=True, 
        save_path=None,
        save_name="Figure23",
        save_extension="png",
        return_data=False
    ):
    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_GHG_final_demand_breakdown_{year}.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        )
        .droplevel(["emission", "unit"])
        .drop(["Exports fob (2)"], axis=1)
        .rename({
            'Changes in inventories': "Changes in\nInventories",
            'Final consumption expenditure by government': "Government",
            'Final consumption expenditure by households': "Households",
            'Final consumption expenditure by non-profit organisations serving households (NPISH)': "Non-profit\nOrganisations",
            'Gross fixed capital formation': "Capital\nFormation"
        }, axis=1)
        .sum(axis=0)
        .sort_values(ascending=False)
        .replace(0, np.nan)
        .dropna()
    )
    household_direct = (
        pd.read_csv(
            result_path / "household_emissions.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        )
        .loc[year]
        .loc["GHG", "value"]
        .drop("Households, totals")
        .sum(axis=0)
    )
    df["Households"] = df["Households"] + household_direct

    fig, ax = plt.subplots(figsize=(8, 4))
    df.plot.bar(stacked=True, ax=ax, width=0.9, edgecolor="k")
    ax.set_xlabel(None)
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_ylabel("CO2-eq [Mt]")
    ax.set_title("GHG footprints by final demand category")
    ax.legend(title=None)
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df

# %%
def figure24_to_26(
        result_path,
        year=config.plot_year,
        final_demand="Final consumption expenditure by government",
        save_plot=True, 
        save_path=None,
        save_name="Figure24",
        save_extension="png",
        return_data=False
    ):

    df = (
        pd.read_csv(
            result_path / "excel_report_format" / f"sector_accounts_agg_GHG_final_demand_breakdown_{year}.tsv",
            sep="\t",
            index_col=[0, 1, 2]
        )
        .loc["GHG"]
        .loc["Mt", final_demand]
    )

    if final_demand == "Final consumption expenditure by households":
        df_direct = (
            pd.read_csv(
                result_path / "household_emissions.tsv",
                sep="\t",
                index_col=[0, 1, 2]
            )
            .loc[year]
            .loc["GHG", "value"]
            .drop("Households, totals")
            .rename(lambda x: x + " (Direct)")
            .rename_axis(["sector"])
        )
        df = pd.concat(
            [df_direct, df]
        )
    mapper = {
        'Changes in inventories': "Changes in inventories",
        'Final consumption expenditure by government': "Government consumption",
        'Final consumption expenditure by households': "Household consumption",
        'Final consumption expenditure by non-profit organisations serving households (NPISH)': "Non-profit organisations",
        'Gross fixed capital formation': "Capital formation"
    }
    name = mapper[final_demand]
    df_plot = (
        df
        .sort_values()
        .replace(0, np.nan)
        .dropna()
        .rename(shorten_name, axis=0)
        .round(2)
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    df_plot.plot.barh(ax=ax, width=0.8, edgecolor="k")
    ax.set_ylabel(None)
    ax.set_title(f"{name} - GHG emissions")
    plt.setp(ax.spines.values(), color="k", linewidth=1)
    ax.set_xlabel("CO2-eq [Mt]")
    if save_plot:
        if save_path == None:
            save_path = result_path / "figures"
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path / f"{save_name}.{save_extension}", bbox_inches="tight")
    fig.show()

    if return_data:
        return df

figure24_to_26(
    cafean_result_path,
    save_plot=False,
    final_demand="Final consumption expenditure by households",
    return_data=True
).sort_values()

# %%
# Wrapper function to make all plots.
def make_all_plots(result_path, year=config.plot_year, save_path=cafean_save_path):
    figure5(result_path, year=year, save_path=save_path)
    figure6(result_path, year=year, save_path=save_path)
    figure7(result_path, year=year, save_path=save_path)
    figure8_and_9(result_path, save_path=save_path)
    figure10(result_path, save_path=save_path)
    figure11_to_13_and_27(result_path, ghg_gas="CO2", save_name="Figure11", save_path=save_path)
    figure11_to_13_and_27(result_path, ghg_gas="CH4", save_name="Figure12", save_path=save_path)
    figure11_to_13_and_27(result_path, ghg_gas="N2O", save_name="Figure13", save_path=save_path)
    figure11_to_13_and_27(result_path, ghg_gas="Biomass CO2", save_name="Figure27", save_path=save_path)
    figure14(result_path, year=year, save_path=save_path)
    figure15(result_path, year=year, save_path=save_path)
    figure16(result_path, year=year, save_path=save_path)
    figure17(result_path, year=year, save_path=save_path)
    figure18(result_path, year=year, save_path=save_path)
    figure19(result_path, year=year, save_path=save_path)
    figure20(result_path, year=year, save_path=save_path)
    table5_to_7(result_path, year_end=year, ghg_gas="CO2", save_name="Table5", save_path=save_path)
    table5_to_7(result_path, year_end=year, ghg_gas="CH4", save_name="Table6", save_path=save_path)
    table5_to_7(result_path, year_end=year, ghg_gas="N2O", save_name="Table7", save_path=save_path)
    figure21(result_path, year=year, save_path=save_path)
    figure22(result_path, year=year, save_path=save_path)
    figure23(result_path, year=year, save_path=save_path)
    figure24_to_26(result_path, year=year, final_demand="Final consumption expenditure by government", save_name="Figure24", save_path=save_path)
    figure24_to_26(result_path, year=year, final_demand="Final consumption expenditure by households", save_name="Figure25", save_path=save_path)
    figure24_to_26(result_path, year=year, final_demand="Gross fixed capital formation", save_name="Figure26", save_path=save_path)


# %%
# Make plots for current version
make_all_plots(cafean_result_path, year=config.plot_year, save_path=cafean_save_path)

# %%
# Make plots for 2020 for comparison with original CaFEAN report
make_all_plots(cafean_result_path, year=2020, save_path=(cafean_save_path / "2020 Comparison"))
#%%
# Make plots for legacy version
make_all_plots(legacy_cafean_result_path, save_path=None, year=2020)

# %%
