# -*- coding: utf-8 -*-
# %% [markdown]
#  # Functions used across the scripts in the CaFEAN project.
#
# %% [markdown]
# Copyright (C) 2024 XIO Sustainability Analytics, Inc
#
# Written by
#
# - Richard Wood
# - Konstantin Stadler
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

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl

def aggregate_chemical_sectors(
        df,
        axis: int = 1,
        level: int =0
    ) -> pd.DataFrame | pd.Series :
    """Aggregates the three chemical sectors to one sector.

    Args:
        df (pd.DataFrame | pd.Series): Accounts in Norwegian classification in either the columns or rows.
        axis (int, optional): Whether the classification is in the rows (=0) or columns (=1). Defaults to 1.
        level (int, optional): The level in the row or column index. Defaults to 0.

    Returns:
        pd.DataFrame | pd.Series: Account with the three chemical sectors aggregated to one sector. 
                        Rest is left unchanged.
    """
    chem_sectors_to_sum = [
        "Coke and refined petroleum products",
        "Chemicals and chemical products",
        "Basic pharmaceutical products and pharmaceutical preparations",
    ]

    chem_sector = ["Basic pharmaceutical products and pharmaceutical preparations"]
    if axis == 1:
        sum_sectors_pos = (
            df
            .columns
            .isin(chem_sectors_to_sum, level=level)
        )
        agg_sector_pos = (
            df
            .columns
            .isin(chem_sector, level=level)
        )
        chem_sum = df.loc[:, sum_sectors_pos].sum(axis=1)
        
        df.loc[:, sum_sectors_pos] = 0
        # ! Due to weird behaviour in multiindex this conditional statement is needed.
        if level == 0:
            df.loc[:, agg_sector_pos] = pd.DataFrame(chem_sum, columns=chem_sector)
        else: 
            df.loc[:, agg_sector_pos] = chem_sum
    if axis == 0:
        sum_sectors_pos = (
            df
            .index
            .isin(chem_sectors_to_sum, level=level)
        )
        agg_sector_pos = (
            df
            .index
            .isin(chem_sector, level=level)
        )
        if type(df) == pd.DataFrame:
            chem_sum = df.loc[sum_sectors_pos, :].sum(axis=0)
            df.loc[sum_sectors_pos, :] = 0
            df.loc[agg_sector_pos, :] = chem_sum.values      
        elif type(df) == pd.Series:
            chem_sum = df.loc[sum_sectors_pos].sum(axis=0)
            df.loc[sum_sectors_pos] = 0
            df.loc[agg_sector_pos] = chem_sum   
    return df


def load_survey(
        folder_path: Path,
        name: str,
        n_index: int,
        index_names: list = ["code", "sector", "category_code", "category", "subcategory_code", "subcategory"],
        file_type: str = "xlsx"
    ) -> pd.DataFrame:
    """Loads SSB consumer expenditure surveys.

    Args:
        folder_path (Path): Path to the folder of the survey.
        name (str): Name of the survey.
        n_index (int): Number of index levels that is present in survey.
        index_names (list, optional): Names of the indexes. Defaults to ["code", "sector", "category_code", "category", "subcategory_code", "subcategory"].
        file_type (str, optional): Whether the file type is xlsx or csv. Defaults to "xlsx".

    Returns:
        pd.DataFrame: SSB consumer expenditure survey.
    """
    index_levels = list(range(n_index))
    index_names = index_names[:n_index]
    survey = (
        pd.read_excel(
            folder_path / f"{name}.{file_type}",
            index_col=index_levels,
            header=[0, 1, 2, 3],
        )
    )
    survey.name = survey.columns.get_level_values(0)[0]
    print(f"Loading survey:\n{survey.name}")

    survey = (
        survey
        .droplevel([0, 1, 3], axis=1)
        .dropna()
        .rename_axis(index_names, axis=0)
        .rename_axis(["variable"], axis=1)
    )

    survey_codes = (
        survey.index.get_level_values("code")
        .str
        .split(".", expand=True)
        .to_frame()
        .rename(lambda x: f"code_level_{x+1}", axis=1)
        .reset_index(drop=True)
        .fillna(0)
    )

    survey["code_level"] = (
        survey
        .index.get_level_values("code")
        .map(lambda x: len(x.split(".")))
    )
    survey = (
        survey
        .set_index(["code_level"], append=True)
    )

    survey  = (
        pd.concat(
            [survey.reset_index(), survey_codes],
            axis=1
        ).set_index(
            list(survey.index.names[::-1]) 
            + list(survey_codes.columns)
        )
        .replace("..", 0)
    )
    return survey


def aggregate_results(
        df: pd.DataFrame,
        agg_dict: dict,
        level_to_agg: str
    ) -> pd.DataFrame:
    """Aggregates results based on a specific column in dataset.

    Args:
        df (pd.DataFrame): Result dataset to be aggregated.
        agg_dict (dict): Dictionary used to aggregate the column.
        level_to_agg (str): Name of column to aggregate.

    Returns:
        pd.DataFrame: Aggregated results.
    """
    assert (
        level_to_agg in df.columns
    ), f"{level_to_agg} aggregation requested but no {level_to_agg} column in source file"
    
    assert (
        "value" in df.columns
    ), "Sector aggregation requested but no value column in source file"

    col = [c for c in df.columns if c != "value"]

    agg_df = (
        df.set_index(col)
        .value.unstack(level_to_agg)
        .T.groupby(agg_dict, sort=False)
        .sum()
        .T.stack(level_to_agg, sort=False)
    )

    agg_df.name = "value"
    agg_df = agg_df.reset_index()[df.columns]

    return agg_df


def diagonalise_series(x: pd.Series) -> pd.DataFrame:
    """Simple function to diagonalise a Pandas series.

    Args:
        x (pd.Series): Pandas series to be diagonalised.

    Returns:
        pd.DataFrame: Diagonal version of pandas series.
    """
    x_diagonal = pd.DataFrame(
        data=np.diag(x),
        index=x.index,
        columns=x.index
    )
    return x_diagonal


def analyse_sector(
        year: int, 
        dataset: pd.DataFrame,
        final_sector: int=1, 
        source_share_limit: int=5, 
        variable: str=None,
    ) -> pd.DataFrame:
    """Provides a detailed analysis of specific sector.

    Args:
        year (int): The year for which to analyse the sector.
        dataset (pd.DataFrame): The result dataset. Either aggregated or not.
        final_sector (int, optional): The name of the sector to be analysed. If integer n is given, it takes the largest n'th sector. Defaults to 1.
        source_share_limit (int, optional): The minimum share (in percentage) of a source sectors share of the total footprint for it to be included in the results. Defaults to 5.
        variable (str, optional): Whether to only consider domestic, imports, or both (None). Defaults to None.

    Raises:
        TypeError: Returned if final_sector keyword is not an suitable integer or string.

    Returns:
        pd.DataFrame: Subset of the input dataset for detailed analysis.
    """
    if variable:
        row_filter = (year, variable)
    else: 
        row_filter = year

    largest_footprints = (
        dataset
        .loc[row_filter, "Footprint"]
        .groupby(["Final sector"])
        .sum()
        .sort_values(ascending=False)
        .round(2)
    )
    if type(final_sector) is int:
        final_sector_name = (
            largest_footprints.index[final_sector-1]
        )
        final_sector_footprint = (
            largest_footprints.values[final_sector-1]
        )
    elif type(final_sector) is str:
        final_sector_name = final_sector
        final_sector_footprint = (
            largest_footprints.loc[final_sector]
        )
    else:
        print("final_sector argument must be a integer or string")
        raise TypeError

    print(f"Sector: {final_sector_name}")
    if variable:
        print(f"{variable} footprint of sector in year {year}: {final_sector_footprint}")
    else:
        print(f"Total footprint of sector in year {year}: {final_sector_footprint}")
    
    df = (
        dataset
        .loc[
            (
                dataset.index.isin(
                    [final_sector_name],
                    level="Final sector"
                )
            ) 
            & (
                dataset["Source sector footprint share"]
                > source_share_limit
            )
        ]
        .loc[year]
        .reorder_levels(["Final sector", "Variable", "Source sector"])
        .sort_values("Source sector footprint share", ascending=False)
        .sort_index(level=["Final sector", "Variable"], sort_remaining=False)
    )
    return df


def make_results(
        year: int,
        domestic_footprint_ts: pd.DataFrame,
        import_footprint_ts: pd.DataFrame,
        S_results: pd.DataFrame,
        source_sector_data: pd.DataFrame,
        source_sector_change_data: pd.DataFrame,
        sector_aggregator: dict,
        aggregate_sectors: bool=True
    ) -> pd.DataFrame:
    """Combines the different data to create combined result dataset.

    Args:
        year (int): Base year to compare against.
        domestic_footprint_ts (pd.DataFrame): Timeseries dataset of domestic footprint.
        import_footprint_ts (pd.DataFrame): Timeseries of emissions emboddied in imports.
        S_results (pd.DataFrame): Timeseries of stressor data.
        source_sector_data (pd.DataFrame): Data on source sectors.
        source_sector_change_data (pd.DataFrame): Data on the changes of the source sectors.
        sector_aggregator (dict): Keys to aggregate.
        aggregate_sectors (bool, optional): Whether or not to aggregate the sectors. Defaults to True.

    Returns:
        pd.DataFrame: Combined result dataset.
    """
        
    footprint = (
        pd.concat(
            [
                domestic_footprint_ts,
                import_footprint_ts
            ],
            keys=["Domestic", "Import"],
            names=["Variable"],
            axis=0
        )
    )
    if aggregate_sectors:
        footprint = (
            footprint
            .rename(sector_aggregator, axis=1)
        )
    footprint = (
        footprint
        .T.groupby(level="Final sector").sum().T
        .stack()
        .reorder_levels(["Year", "Variable", "Final sector", "Source sector"])
    )

    footprint_pivot = (
        footprint
        .unstack("Year")
    )

    footprint_difference = (
        footprint_pivot
        .subtract(
            footprint_pivot.loc[:, year],
            axis=0)
        .unstack(["Variable", "Final sector", "Source sector"])
    )
    footprint_percentage_change = (
        footprint_difference
        .unstack("Year")
        .div(
            footprint_pivot.loc[:, year],
            axis=0
        )
        .replace([np.nan, np.inf, -np.inf], 0)
        .unstack(["Variable", "Final sector", "Source sector"])
        .mul(100)
        .round(2)
    )

    footprint_shares = (
        footprint
        .div(
            footprint
            .groupby(["Year", "Final sector"])
            .transform("sum")
        )
        .replace([np.nan, np.inf, -np.inf], 0)
        .mul(100)
        .round(2)
    )

    footprint_shares_difference = (
        footprint_shares
        .unstack("Year")
        .subtract(
            footprint_shares
            .unstack("Year")
            .loc[:, year], axis=0
        )
        .unstack(["Variable", "Final sector", "Source sector"])
    )

    results = (
        pd.concat(
            [
                footprint,
                footprint_shares,
                footprint_shares_difference,
                footprint_difference,
                footprint_percentage_change
            ],
            keys=[
                "Footprint",
                "Source sector footprint share",
                f"Absolute difference in source sector footprint share from {year}",
                f"Absolute difference in footprint from {year}",
                f"Percentage change in footprint from {year}"
            ],
            names=["Indicator"],
            axis=1
        )
        .reset_index()
        .merge(
            S_results.reset_index(),
            on=["Year", "Variable", "Source sector"],
            how="left"
        )
        .merge(
            source_sector_data.reset_index(),
            on=["Year", "Variable", "Source sector"],
            how="left"
        )
        .merge(
            source_sector_change_data.reset_index(),
            on=["Year", "Variable", "Source sector"],
            how="left"
        )
        .set_index(footprint.index.names)
    )
    return results


def survey_aggregate_table(
        survey: pd.DataFrame,
        variable: str,
        mappers: dict,
        column_order: list=None,
        show_table: bool=True,
        save_table: bool=True,
        save_path: Path | str=None
    ) -> pd.DataFrame:
    """Aggregates household footprint into a easy interpretable format with highlighted backgrounds.

    Args:
        survey (pd.DataFrame): Survey footprint result.
        variable (str): _description_
        mappers (dict): Keys used to aggregate survey.
        column_order (list, optional): Order of the columns. Defaults to None.
        show_table (bool, optional): Whether or not to return table. Defaults to True.
        save_table (bool, optional): Whether or not to save the table. Defaults to True.
        save_path (Path | str, optional): Where to save the table. Defaults to None.

    Returns:
        pd.DataFrame: Aggregates household footprint table. Only returned if show_table=True.
    """
    if len(survey) < 100:
        mapper = mappers["map_2_to_1"].copy()
    elif len(survey) < 200:
        mapper = mappers["map_3_to_1"].copy()

    table = (
        survey
        .rename(mapper, axis=0)
        .groupby("sector").sum()
        .stack("variable")
        .reorder_levels([1, 0])
        .sort_index()
        .loc[variable]
    )
    table_row_order = (
        table
        .drop("Total", axis=0)
        .sum(axis=1)
        .sort_values(ascending=False)
        .index
        .values
    )
    table_row_order = np.concatenate((table_row_order, ["Total"]))

    min_color = (11/255, 187/255, 183/255)
    max_color = (6/255, 98/255, 97/255)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "Miljødir",
        [min_color, max_color]
    )
    if column_order is None:
        table_column_order = (
            table.sum(axis=0)
            .sort_values(ascending=False)
            .index
        )
    else:
        table_column_order = column_order

    table_sorted = (
        table
        .loc[table_row_order, table_column_order]
        .style
        .background_gradient(axis=1, cmap=cmap)
        .format(precision=2)
    )

    if save_table:
        table_sorted.to_excel(
            save_path 
            / f"{survey.name} ({variable} - Aggregated).xlsx",
        )

    if show_table:
        return table_sorted
