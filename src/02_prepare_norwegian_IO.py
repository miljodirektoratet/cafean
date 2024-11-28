# -*- coding: utf-8 -*-
# %% [markdown]
#  # Norwegian IO preparation script
#
# This script prepares the Norwegian IO data for general use.
#
# It processes the raw excel files and stores the data in easy machine readable format
# for further processing. 

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

# %% [markdown]
# Here we define the path to the Norwegian IO data
# %%
norwegian_IO_raw_path: Path = (
    config.data_path
    / "01_raw"
    / "SSB"
    / "norwegian IO"
)

domestic_IO_raw_path: Path = (
    norwegian_IO_raw_path
    / "domestic"
)

import_IO_raw_path: Path = (
    norwegian_IO_raw_path
    / "import"
)

supply_tables_raw_path: Path = (
    norwegian_IO_raw_path
    / "supply tables"
)

use_tables_raw_path: Path = (
    norwegian_IO_raw_path
    / "use tables"
)

# %% [markdown]
# Finally, we specify the folder in which we save 
# the processed Norwegian IO data.

# %%
norwegian_IO_interim_path: Path = (
    config.data_path
    / "02_interim" 
    / "norwegian IO"
)

norwegian_IO_interim_path.mkdir(exist_ok=True, parents=True)

# %%
# "Changes in valuables" are for some reason written in different ways in the different tables.
# Here we make a simple lambda function to correct for it.
change_in_valuables_name_fix = (
    lambda x: "Changes in valuables" if "Changes in valuables" in x else x
)

real_estate_sector_name_fix = {            
    "Real estate activities (excluding imputed rents)": "Real estate services (excluding imputed rents)"
}


# %% [markdown]
# # Script

# %%
# We iterate over each year, to then load, extract,
# and store the different parts of the datasets.

for year in config.years_io:
    # Create folder for each year
    save_path: Path = (
        norwegian_IO_interim_path 
        / str(year)
    )
    save_path.mkdir(exist_ok=True, parents=True)

    # IO domestic data
    domestic_IO = (
        pd.read_excel(
            domestic_IO_raw_path/ f"iot_1850_{year}.xlsx",
            index_col=[0, 1, 2, 3],
            skiprows=23,
            skipfooter=30,
            header=[0, 1, 2, 3, 4]
        )
        .droplevel(1)
        .rename(lambda x: x.strip(), level=2, axis=0)
        .rename(lambda x: x.strip(), level=3, axis=1)
        .rename(
            real_estate_sector_name_fix,
            axis=1,
            level=3
        )
    )

    # Aggregate chemical sectors
    domestic_IO = aggregate_chemical_sectors(domestic_IO, axis=1, level=3)
    domestic_IO = aggregate_chemical_sectors(domestic_IO, axis=0, level=2)

    # Extract row and column levels we want to use for filtering
    transaction_level = domestic_IO.columns.get_level_values("Transaction")
    row_code_0 = domestic_IO.index.get_level_values(0)
    row_code_1 = domestic_IO.index.get_level_values(1)
    row_code_2 = domestic_IO.index.get_level_values(2)

    Y_dom_all = (
        domestic_IO.loc[
            (row_code_0 == "Transaction\n(3. quadrant)")
            & ~(row_code_1.isin(config.value_added_rows))
            ,
            ~(transaction_level == "P2")
            & ~(transaction_level.isin(config.agg_final_demand_columns))
        ]
        .droplevel([0, 1], axis=0)
        .droplevel([0, 1, 2, 4], axis=1)
        .rename_axis(["sector"], axis=0)
        .rename_axis(["category"], axis=1)
        .drop("Total", axis=0)
        .rename(change_in_valuables_name_fix, axis=1)
    )

    Y_dom = (
        Y_dom_all
        .loc[:, config.final_demand_categories]
        .replace(np.nan, 0)
    )

    Y_dom.to_excel(
        save_path / "Y_domestic.xlsx"
    )

    Z_dom = (
        domestic_IO
        .loc[
            (row_code_0.isin(["Transaction\n(3. quadrant)"], level=0))
            & ~(row_code_1.isin(['R', 'RNAM', 'RNTS', 'RZ', "RADJ"])), 
            transaction_level.isin(["P2"], level="Transaction")
        ]
        .droplevel([0, 1], axis=0)
        .droplevel([0, 1, 2, 4], axis=1)
        .drop("Total", axis=1)
        .replace(np.nan, 0)
        .rename_axis(["sector"], axis=0)
        .rename_axis(["sector"], axis=1)
    )

    Z_dom.to_excel(
        save_path / "Z_domestic.xlsx"
    )

    value_added = (
        domestic_IO
        .loc[
            row_code_2.isin(["Value added at basic prices"], level=0),
            transaction_level.isin(["P2"], level="Transaction")
        ]
        .droplevel([0, 1], axis=0)
        .droplevel([0, 1, 2, 4], axis=1)
        .drop("Total", axis=1)
        .replace(np.nan, 0)
        .rename_axis(["sector"], axis=0)
        .rename_axis(["sector"], axis=1)
        .iloc[0, :]
        .T
    )

    value_added.to_excel(
        save_path / "value_added.xlsx"
    )

    value_added_detailed = (
        domestic_IO
        .loc[
            row_code_1.isin(config.value_added_rows + ["RZ"]),
            transaction_level.isin(["P2"], level="Transaction")
        ]
        .droplevel([0, 1], axis=0)
        .droplevel([0, 1, 2, 4], axis=1)
        .drop("Total", axis=1)
        .replace(np.nan, 0)
        .rename_axis(["sector"], axis=0)
        .rename_axis(["sector"], axis=1)
    )

    value_added_detailed.to_excel(
        save_path / "value_added_detailed.xlsx"
    )

    # IO import data
    imp = (
        pd.read_excel(
            import_IO_raw_path / f"iot_1950_{year}.xlsx",
            index_col=[0, 1, 2, 3],
            skiprows=23,
            skipfooter=43,
            header=[0, 1, 2, 3, 4]
        )
        .droplevel([0, 1, 2])
        .rename_axis(["sector"])
        .rename(lambda x: x.strip(), level=0, axis=0)
        .rename(lambda x: x.strip(), level=3, axis=1)
        .rename(
            real_estate_sector_name_fix,
            axis=1,
            level=3
        )
    )
    imp = aggregate_chemical_sectors(imp, axis=1, level=3)
    imp = aggregate_chemical_sectors(imp, axis=0, level=0)

    transaction_level = imp.columns.get_level_values("Transaction")

    Y_imp_all = (
        imp
        .loc[
            :, 
            ~(transaction_level == "P2")
            & ~(transaction_level.isin(config.agg_final_demand_columns))
        ]
        .droplevel([0, 1, 2, 4], axis=1)
        .rename_axis(["category"], axis=1)
        .drop("Total", axis=0)
        .rename(change_in_valuables_name_fix, axis=1)
    )

    Y_imp = (
        Y_imp_all
        .loc[:, config.final_demand_categories]
        .replace(np.nan, 0)
    )

    Y_imp.to_excel(
        save_path / "Y_import.xlsx"
    )

    Z_imp = (
        imp
        .loc[
            (imp.index != "Total"), 
            (transaction_level == "P2")
        ]
        .droplevel([0, 1, 2, 4], axis=1)
        .drop("Total", axis=1)
        .replace(np.nan, 0)
        .rename_axis(["sector"], axis=0)
        .rename_axis(["sector"], axis=1)
    )

    Z_imp.to_excel(
        save_path / "Z_import.xlsx"
    )

    # Downloads SUTs for subset of years that are needed for more advanced analysis.
    if year in config.years_sut:
        # Supply data
        supply_data = (
            pd.read_excel(
                supply_tables_raw_path / f"iot_1500_{year}.xlsx",
                index_col=[0, 1, 2],
                skiprows=23,
                skipfooter=28,
                header=[0, 1, 2, 3, 4]
            )
            .droplevel([0])
            .rename_axis(["code", "sector"])
            .rename(lambda x: x.strip(), level=1, axis=0)
            .rename(lambda x: x.strip(), level=3, axis=1)
        )

        supply_data = aggregate_chemical_sectors(supply_data, axis=0, level=1)

        transaction_level = supply_data.columns.get_level_values("Transaction")
        row_code_level = supply_data.index.get_level_values("code")

        intermediate_supply = (
            supply_data
            .loc[
                ~row_code_level.isin(config.supply_tables_extra_row_codes),
                transaction_level == "P1"
            ]
            .droplevel(0, axis=0)
            .droplevel([0, 1, 2, 4], axis=1)
            .drop("Total", axis=1)
        )

        intermediate_supply.to_excel(
            save_path / "supply.xlsx"
        )

        extra_supply_data = (
            supply_data
            .loc[
                row_code_level.isin(config.supply_tables_extra_row_codes),
                transaction_level != "P7"
            ]
            .droplevel(0, axis=0)
            .droplevel([0, 1, 2, 4], axis=1)
            .drop("Total", axis=1)
            .rename({"": "... of which"})
        )
        extra_supply_data.to_excel(
            save_path / "extra_supply_data.xlsx"
        )

        valuation_layers = (
            supply_data
            .loc[
                ~row_code_level.isin(config.supply_tables_extra_row_codes),
                ~(transaction_level.isin(["P1", "P7"]))
            ]
            .droplevel(0, axis=0)
            .droplevel([0, 1, 2, 4], axis=1)
            .rename_axis(["sector"], axis=1)
        )

        valuation_layers.to_excel(
            save_path / "valuation_layers.xlsx"
        )

        # Use data
        use_data = (
            pd.read_excel(
                use_tables_raw_path / f"iot_1600_{year}.xlsx",
                index_col=[0, 1, 2, 3, 4, 5, 6],
                skiprows=23,
                skipfooter=18,
                header=[0, 1, 2, 3, 4]
            )
            .droplevel([0, 1, 2, 3, 4])
            .rename_axis(["code", "sector"])
            .rename(lambda x: x.strip(), level=1, axis=0)
            .rename(lambda x: x.strip(), level=3, axis=1)
        )

        transaction_level = use_data.columns.get_level_values("Transaction")
        row_code_level = use_data.index.get_level_values("code")

        intermediate_use = (
            use_data
            .loc[
                ~row_code_level.isin(config.use_tables_extra_row_codes),
                transaction_level == "P2"
            ]
            .droplevel(0, axis=0)
            .droplevel([0, 1, 2, 4], axis=1)
            .drop("Total", axis=1)
        )

        intermediate_use.to_excel(
            save_path / "intermediate_use.xlsx"
        )

        final_demand_use = (
            use_data
            .loc[
                ~row_code_level.isin(config.use_tables_extra_row_codes),
                ~(transaction_level == "P2")
            ]
            .droplevel(0, axis=0)
            .droplevel([0, 1, 2, 4], axis=1)
            .rename(change_in_valuables_name_fix, axis=1)
            .loc[:, config.final_demand_categories]
        )

        final_demand_use.to_excel(
            save_path / "final_demand_use.xlsx"
        )

        total_uses = (
            use_data
            .loc[
                ~row_code_level.isin(config.use_tables_extra_row_codes),
                transaction_level.isin(["TFUPR", "TUPR"])
            ]
            .droplevel(0, axis=0)
            .droplevel([0, 1, 3, 4], axis=1)
        )
        total_uses.to_excel(
            save_path / "total_use.xlsx"
        )

        # TODO: Extract the extra use data available in dataset.
        # extra_use_data = (
        #     use_data
        #     .loc[
        #         row_code_level.isin(use_tables_extra_row_codes),
        #         : 
        #     ]
        #     #.droplevel(0, axis=0)
        #     .droplevel([0, 1, 2, 4], axis=1)
        #     #.rename(change_in_valuables_name_fix, axis=1)
        #     #.loc[:, final_demand_categories]
        # )


# %%
