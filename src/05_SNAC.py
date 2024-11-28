# -*- coding: utf-8 -*-
# %% [markdown]
#  # Hybrid SNAC model adopted for Norway SSB data
#
# Top level script for calculation Norwegian carbon footprints
# based on the hybrid-SNAC approach, coupling EXIOBASE 3
# with official data from SSB Norway
#

# %% [markdown]
# Copyright (C) 2024 XIO Sustainability Analytics, Inc
#
# Written by
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
# Import own functions
# %%
from functions.general import (
    aggregate_chemical_sectors,
    aggregate_results
)
from functions.snac import (
    snac_findem_breakdown, 
    calc_snac_sector_accounts,
    snac_coupling,
    stack_and_unit
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

exiobase_interim_path: Path = (
    config.data_path 
    / "02_interim"
    / "exiobase"
)

cafean_result_path: Path = (
    config.data_path
    / "03_results"
    / "CaFEAN"
)
cafean_result_path.mkdir(parents=True, exist_ok=True)

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
industry_emissions = (
    aggregate_chemical_sectors(
        industry_emissions,
        axis=1, 
        level=0
    )
)

# %%
col_dom_io = {}
col_Q_imp_exio = {}
col_snac = {}
col_findem_breakdown = {}
col_sector_accounts = {}
col_grand_totals = {}
col_sources = {}
col_sources_tot = {}
col_sources_tot_with_hhld = {}
col_household_footprint = {}

# %%

for year in config.years_io:
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

    x_import = (
        Z_import.sum(axis=1)
        + Y_import_breakdown.sum(axis=1)
    )

    A_domestic = pymrio.calc_A(Z_domestic, x_domestic)
    A_import = pymrio.calc_A(Z_import, x_domestic)
    L_domestic = pymrio.calc_L(A_domestic)

    F_domestic = industry_emissions.loc[year]
    S_domestic = pymrio.calc_S(F_domestic, x_domestic)

    SL_domestic = pymrio.calc_M(S=S_domestic, L=L_domestic) 

    dom_io = dict()
    dom_io["Zdom"] = Z_domestic
    dom_io["Ydom_exp"] = Y_domestic_export
    dom_io["Ydom_br"] = Y_domestic_breakdown
    dom_io["Ydom"] = Y_domestic
    dom_io["Zimp"] = Z_import
    dom_io["Yimp_exp"] = Y_import_export
    dom_io["Yimp_br"] = Y_import_breakdown
    dom_io["Yimp"] = Y_import
    dom_io["x"] = x_domestic
    dom_io["Adom"] = A_domestic
    dom_io["Aimp"] = A_import
    dom_io["Ldom"] = L_domestic
    dom_io["Sdom"] = S_domestic
    dom_io["SLdom"] = SL_domestic

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

    # Fix emission mismatch
    # Due to confidentiality of data for the three sectors below in both Norwegian emissions and economic data,
    # the data is aggregated to one sector, to be consistent with the Norwegian IO data
    import_data = (
        aggregate_chemical_sectors(import_data, axis=0, level=0)
    )

    x_import_exiobase = (
        import_data.loc[:, "imports"]
    )

    F_import = (
        import_data
        .drop("imports", axis=1, level=0)
    )
    F_import_unit = (
        F_import
        .columns
        .to_frame()
        .reset_index(drop=True)
        .rename({"variable": "emission"}, axis=1)
    )
    # We have Q_ something as abbreviation for the (full-supply-chain) multipliers
    # throughout the script. Note, that these are called .M in the pymrio package.
    Q_import = (
        F_import
        .div(
            x_import_exiobase.loc[:, "M.NOK"],
            axis=0
        )
        .replace([np.nan, np.inf, -np.inf], 0)
        .rename(lambda x: f"{x}/M.NOK", axis=1, level=1)
        .rename_axis(["emission", "unit"], axis=1)
    )
    Q_import_unit = (
        Q_import.columns
        .to_frame()
        .reset_index(drop=True)
        .rename({"account": "emission"}, axis=1)
    )
    Q_import = (
        Q_import
        .droplevel("unit", axis=1)
        .T
    )

    snac = snac_coupling(
        Aimp=dom_io["Aimp"],
        Ldom=dom_io["Ldom"],
        Qimp_mrio=Q_import
    )

    findem_breakdown_gather = dict()
    for emission in config.emission_types + ["GHG"]:
        SLdom_diag = pd.DataFrame(
            np.diag(dom_io["SLdom"].loc[emission, :]),
            index=dom_io["SLdom"].columns,
            columns=dom_io["SLdom"].columns,
        )

        QAimp_Ldom_diag = pd.DataFrame(
            np.diag(snac["QAimpLdom"].loc[emission, :]),
            index=snac["QAimpLdom"].columns,
            columns=snac["QAimpLdom"].columns,
        )
        Q_imp_exio_diag = pd.DataFrame(
            np.diag(Q_import.loc[emission, :]),
            index=Q_import.columns,
            columns=Q_import.columns,
        )
        findem_breakdown_gather[emission] = snac_findem_breakdown(
            SLdom_diag=SLdom_diag,
            QAimp_Ldom_diag=QAimp_Ldom_diag,
            Q_imp_exio_diag=Q_imp_exio_diag,
            Ydom_br=dom_io["Ydom_br"],
            Yimp_br=dom_io["Yimp_br"],
        )

    _df = pd.concat(findem_breakdown_gather, axis=1)
    _df.columns.names = ["emission"]

    fp_findem_breakdown = _df.stack()
    del _df

    # footprint of imports (goods traded at the border)
    total_imports = dom_io["Zimp"].sum(axis=1) + dom_io["Yimp"] + dom_io["Yimp_exp"]

    sector_accounts = calc_snac_sector_accounts(
        SLdom=dom_io["SLdom"],
        QAimp_Ldom=snac["QAimpLdom"],
        Ydom_diag=pd.DataFrame(
            data=np.diag(dom_io["Ydom"]),
            index=dom_io["Ydom"].index,
            columns=dom_io["Ydom"].index,
        ),
        Yimp_diag=pd.DataFrame(
            data=np.diag(dom_io["Yimp"]),
            index=dom_io["Yimp"].index,
            columns=dom_io["Yimp"].index,
        ),
        Ydom_exp_diag=pd.DataFrame(
            data=np.diag(dom_io["Ydom_exp"]),
            index=dom_io["Ydom_exp"].index,
            columns=dom_io["Ydom_exp"].index,
        ),
        Yimp_exp_diag=pd.DataFrame(
            data=np.diag(dom_io["Yimp_exp"]),
            index=dom_io["Yimp_exp"].index,
            columns=dom_io["Yimp_exp"].index,
        ),
        Qimp_mrio=Q_import,
        total_imports=total_imports,
        pba=industry_emissions.loc[year],
    )

    # fp-src calc: get source of footprint from exiobase
    _source = pd.read_excel(
        exiobase_interim_path / "exio_no_bp_raw.xlsx",
        sheet_name=f"imp_source_{year}",
        index_col=[0, 1],
        header=[0],
    )
    _source.columns.name = "region"
    _source = (
        _source
        .rename(exiobase_indudstry_mapper, axis=0, level="sector")
        .stack("region")
        .unstack("sector")
    )

    # fp-src calc: convert to nor class and do the chem sector adjustment
    source_nor_class = (
        _source.dot(weighted_concordance_normalised)
    )
    source_nor_class = aggregate_chemical_sectors(source_nor_class)

    # fp-src: get src scaled import mulitpliers for each region
    Q_src_imp = (
        source_nor_class.div(x_import_exiobase.loc[:, "M.NOK"], axis=1)
        .fillna(0)
        .replace(np.inf, 0)
    )

    # fp-src: multiply by imports and dom demand
    Q_Aimp_src = Q_src_imp @ dom_io["Aimp"]
    Q_Aimp_Ldom_src = Q_Aimp_src @ dom_io["Ldom"]

    # fp-src: multiply by final demand and add norwegian domestic demand
    src_all = Q_Aimp_Ldom_src * dom_io["Ydom"] + Q_src_imp * dom_io["Yimp"]

    src_all = (
        src_all
        .rename_axis(["sector"], axis=1)
        .stack("sector")
        .unstack("region")
    )
    src_all.loc[:, "NO"] = (
        src_all.loc[:, "NO"]
        + (
            (dom_io["SLdom"] * dom_io["Ydom"])
            .rename_axis(["sector"], axis=1)
            .stack("sector")
        )
    )
    src_all = src_all.stack()
    src_all.name = "value"
    src_total = src_all.groupby(level=["region", "emission"]).sum()


    # household footprint
    src_hh_dom = (
        dom_io["SLdom"].mul(
            dom_io["Ydom_br"]
            .loc[:, "Final consumption expenditure by households"]
        )
        .stack()
        .to_frame("value")
        .assign(region="NO")
        .set_index("region", append=True)
        .rename_axis(["emission", "sector", "region"])
    )
    src_hh_imp_dom = (
        Q_Aimp_Ldom_src.mul(
            dom_io["Ydom_br"]
            .loc[:, "Final consumption expenditure by households"]
        )
        .stack()
        .to_frame("value")
        .rename_axis(["emission", "region", "sector"])
    )
    src_hh_imp = (
        Q_src_imp.mul(
            dom_io["Yimp_br"]
            .loc[:, "Final consumption expenditure by households"]
        )
        .stack()
        .to_frame("value")
        .rename_axis(["emission", "region", "sector"])
    )
    src_hh_direct = (
        final_demand_emissions.loc[year]
        .loc[:, "Households, totals"]
        .to_frame("value")
        .assign(sector="Final consumption expenditure by households")
        .assign(region="NO")
        .set_index(["sector", "region"], append=True)
        .rename_axis(["emission", "sector", "region"])
    )
    src_hh_fp = (
        pd.concat(
            [
                src_hh_dom,
                src_hh_imp_dom.reorder_levels(src_hh_dom.index.names),
                src_hh_imp.reorder_levels(src_hh_dom.index.names),
                src_hh_direct.reorder_levels(src_hh_dom.index.names)
            ],
            keys=[
                "domestic",
                "import_domestic",
                "import",
                "direct"
            ],
            names=["component"]
        )
     )

    # calculate some convenient totals
    total_accounts = pd.DataFrame(
        sector_accounts.groupby(["component", "emission"]).sum(),
        columns=["value"],
    )
    total_hhld = (
        pd.DataFrame(
            final_demand_emissions.loc[year, "Households, totals"]
        )
        .rename_axis(["emission"])
    )
    total_hhld.columns = ["value"]
    tot_fp_with_hhld = (
        (total_accounts.loc["total_footprint"] + total_hhld)
        .assign(component="total_footprint_with_households")
        .set_index("component", append=True)
        .reorder_levels(["component", "emission"])
    )
    tot_pba_with_hhld = (
        (total_accounts.loc["production_account"] + total_hhld)
        .assign(component="production_account_with_households")
        .set_index("component", append=True)
        .reorder_levels(["component", "emission"])
    )
    total_hhld = (
        total_hhld.assign(component="household_totals")
        .set_index("component", append=True)
        .reorder_levels(["component", "emission"])
    )
    total_accounts = pd.concat(
        [total_accounts, tot_fp_with_hhld, tot_pba_with_hhld, total_hhld],
        axis=0,
    ).squeeze()

    src_tot_hhld = src_total.copy().unstack('region')
    src_tot_hhld.loc[:, 'NO'] = (src_tot_hhld.loc[:, 'NO'] + 
                                 total_hhld.loc['household_totals', 'value'])
    src_total_with_hhld = src_tot_hhld.stack('region').reorder_levels(['region', 'emission']).sort_index()

    col_dom_io[year] = dom_io
    col_Q_imp_exio[year] = Q_import
    col_snac[year] = snac
    col_findem_breakdown[year] = fp_findem_breakdown.rename_axis(["component", "sector", "category", "emission"])
    col_sector_accounts[year] = sector_accounts.rename_axis(["component", "emission", "sector"])
    col_grand_totals[year] = total_accounts
    col_sources[year] = src_all
    col_sources_tot[year] = src_total
    col_sources_tot_with_hhld[year] = src_total_with_hhld
    col_household_footprint[year] = src_hh_fp.loc[:, "value"]


# %%
# Getting the household emissions directly from the interim data.
emis_nor_hhld = (
        final_demand_emissions
        .rename_axis(["year", "emission"])
        .rename_axis(["hhld_component"], axis=1)
        .stack()
        .to_frame("value")
)

# %% [markdown]
# Stacking the results into a long list format

# %%
sort_sector = dom_io["x"].index.tolist()
emission_unit_mapper = F_import_unit.set_index("emission")["unit"].to_dict()

# %%
all_sector_accounts = stack_and_unit(
    col_sector_accounts, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)
all_findem_breakdown = stack_and_unit(
    col_findem_breakdown, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)
all_hhld = stack_and_unit(
    emis_nor_hhld, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)
all_sources = stack_and_unit(
    col_sources, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)
all_sources_tot = stack_and_unit(
    col_sources_tot, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)

all_sources_tot_with_households = stack_and_unit(
    col_sources_tot_with_hhld, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)

all_totals = stack_and_unit(
    col_grand_totals, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)
all_household_footprint = stack_and_unit(
    col_household_footprint, emission_unit_mapper=emission_unit_mapper, sector_order=sort_sector
)

# %% [markdown]
# Finally we store the results
# %%
all_sector_accounts.to_csv(cafean_result_path / "sector_accounts.tsv", sep="\t")
all_findem_breakdown.to_csv(cafean_result_path / "footprints_final_demand_breakdown.tsv", sep="\t")
all_hhld.to_csv(cafean_result_path / "household_emissions.tsv", sep="\t")
all_totals.to_csv(cafean_result_path / "total_accounts.tsv", sep="\t")
all_sources.to_csv(cafean_result_path / "footprint_sources.tsv", sep="\t")
all_sources_tot.to_csv(cafean_result_path / "footprint_sources_totals.tsv", sep="\t")
all_sources_tot_with_households.to_csv(cafean_result_path / "footprint_sources_totals_with_households.tsv", sep="\t")
all_household_footprint.to_csv(cafean_result_path / "full_household_footprint.tsv", sep="\t")

# %% [markdown]
# File(s) to aggregate.
# The script can also be used to just aggregate a specific file,
# just pass the path to that file within a list.

# %%
src_files = list(cafean_result_path.glob("*.tsv"))

# %% [markdown]
# Suffix for the aggregated file.
# This will be appended to the filename before the extension.
# For example: result.tsv -> result_agg.tsv
# Note: result files will be stored in the same folder as the source file and
# overwrite already existing files with the same name.

# %%
agg_suffix = "_agg"

# %% [markdown]
# Aggregation specification
# This must be a pandas data-frame with columns `src` and `agg`.
# This can be stored anywhere, we provide the ones used in the
# report in the `data` folder.
#
# Set reg_agg or sec_agg to None if aggregation along this dimension is not required.
# One can also manually define the aggregation specification as
# a dictionary (or a pandas data-frame and convert to dict),
# with src as key and agg as value.

# %%
# reg_agg = None
reg_agg = (
    pd.read_csv((
        config.data_path 
        / "00_auxiliary" 
        / "classifications"
        / "region_agg_spec.tsv"
        ), sep="\t")
    .apply(lambda x: x.str.strip())
    .set_index("src")["agg"]
    .to_dict()
)

# %%
sec_agg = (
    pd.read_csv((
        config.data_path 
        / "00_auxiliary" 
        / "classifications"
        / "sector_agg_spec.tsv"
        ), sep="\t")
    .apply(lambda x: x.str.strip())
    .set_index("src")["agg"]
    .to_dict()
)

# %% [markdown]
# Loop over all files and aggregate them

# %%
for fl in src_files:
    df = pd.read_csv(fl, sep="\t")

    data_changed = False
    if (
        "sector" in df.columns
        and sec_agg
        and set(sec_agg.keys()).issubset(set(df.sector))
    ):
        print(f" Aggregating {fl.name} along sector dimension")
        df = aggregate(df, sec_agg, "sector")
        data_changed = True
    if (
        "region" in df.columns
        and reg_agg
        and set(reg_agg.keys()).issubset(set(df.region))
    ):
        print(f" Aggregating {fl.name} along region dimension")
        df = aggregate(df, reg_agg, "region")
        data_changed = True

    if data_changed:
        result_file = fl.parent / (fl.stem + agg_suffix + fl.suffix)
        print(f" Writing aggregated file to {result_file}")
        df.to_csv(result_file, sep="\t", index=False)
    else:
        print(f" No aggregation done for {fl.name}")


# %%
# Folder definitions
report_folder = cafean_result_path / "excel_report_format"
report_folder.mkdir(exist_ok=True)


# %% [markdown]
# Read result files and convert to table format

# %%
_tot_all = pd.read_csv(cafean_result_path / "total_accounts.tsv", sep="\t")
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

all_sources_agg = pd.read_csv(cafean_result_path / "footprint_sources_agg.tsv", sep="\t")
ghg_source_reg_agg = (
    all_sources_agg.query("emission=='GHG'").query(f"year=={config.year_to_extract}").drop("year", axis=1)
)
ghg_source_reg_agg.set_index(
    ["emission", "unit", "sector", "region"], inplace=True, append=True
)
ghg_source_reg_agg = ghg_source_reg_agg.droplevel(0, axis=0)
ghg_source_reg_agg_pivot = ghg_source_reg_agg.unstack("region")
ghg_source_reg_agg_pivot = ghg_source_reg_agg_pivot.droplevel(0, axis=1)
ghg_source_reg_agg_pivot.to_csv(report_folder / f"ghg_source_reg_agg_{config.year_to_extract}.tsv", sep="\t")

emis_sources_totals = (
    pd.read_csv(cafean_result_path / "footprint_sources_totals_with_households_agg.tsv", sep="\t")
)


src_tot_year = (
    emis_sources_totals
    .loc[emis_sources_totals.year == config.year_to_extract, :]
    .set_index(['year', 'emission', 'region', 'unit'])
    .loc[config.year_to_extract, 'value'].unstack('region')
)
src_tot_year.to_csv(report_folder / f"emis_sources_totals_with_households_{config.year_to_extract}.tsv", sep="\t")

# Alternative method to save in excel in different sheets
# with pd.ExcelWriter(report_folder / "emis_sources_totals_with_households_all_years.xlsx") as writer:
#     for year in emis_sources_totals.year.unique():
#         df = emis_sources_totals.loc[emis_sources_totals.year == year, :].set_index(['year', 'emission', 'region', 'unit']).loc[year, 'value'].unstack('region')
#         df.to_excel(writer, sheet_name=str(year))

all_sector_breakdown_agg = pd.read_csv(
    cafean_result_path / "sector_accounts_agg.tsv", sep="\t"
)
for loop_component in all_sector_breakdown_agg.component.unique():
    ghg_sector = (
        all_sector_breakdown_agg.query(f"year=={config.year_to_extract}")
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
        f"sector_accounts_agg_{config.year_to_extract}_" + str(loop_component) + ".tsv"
    )
    ghg_sector_pivot.to_csv(output_file_path, sep="\t")


all_findem_breakdown_agg = pd.read_csv(
    cafean_result_path / "footprints_final_demand_breakdown_agg.tsv", sep="\t"
)
ghg_findem = (
    all_findem_breakdown_agg
    .query("emission=='GHG'")
    .query(f"year=={config.year_to_extract}")
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
    f"sector_accounts_agg_GHG_final_demand_breakdown_{config.year_to_extract}.tsv"
)
ghg_findem_pivot.to_csv(output_file_path, sep="\t")
# %%
print(f"Done - Results stored at {cafean_result_path}")

# %%
