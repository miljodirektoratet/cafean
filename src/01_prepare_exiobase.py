# -*- coding: utf-8 -*-
# %% [markdown]
#  # EXIOBASE download and preparation script
#
# This script prepares the EXIOBASE3 data for use in the SNAC coupling.
#
# It downloads the EXIOBASE3 data from Zenodo, extracts the data and prepares
# it for use in the SNAC coupling. The preparation includes adding the biogenic
# CO2 data and the characterization of the emissions to the Norwegian emission
# data. The script also calculates the Norwegian imports and the Norwegian
# final demand for each year. Finally, it calculates the Norwegian emissions
# based on EXIOBASE for each year, both in normal and diagonalized/per source format.
#
# In the default setting, the EXIOBASE data is downloaded to ./data/exiobase within
# the project repository. This folder is not tracked by version control (added
# to .gitignore). The script works with both, the original compressed zip
# downloads or extracted folders.
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
# Import of required external packages.
# These are all available through pip/conda/mamba install.
# %%
import pandas as pd
import numpy as np
import pymrio
import pymrio.tools.iomath as iomath

# %% [markdown]
# Import Python internal packages 
# (These don't need to be installed, part of standard python)
# %%
from pathlib import Path

# %% [markdown]
# ### Locations / folder definitions

# %% [markdown]
# Here we define the folder for the EXIOBASE raw data.
# This can be anywhere on the disk (in case EXIOBASE is available locally already).
# By default, we set it to ./exiobase within the project repository. This folder
# is not tracked by version control (added to .gitignore).
# The script works with both, the original
# compressed zip downloads or extracted folders.
#
# Update: The new version of EXIOBASE comes with extensions in separate folders.
# Hence we need to specify paths to each of these. 
# These include biogenic CO2.
# %%
exio_raw_path: Path = (
    config.data_path
    / "01_raw"
    / "exiobase"
)

exio_core_path: Path = (
    config.exiobase_raw_path
    / "core"
)

exio_non_combustion_ghg_path: Path = (
    config.exiobase_raw_path
    / "extension"
    / "air_emissions_non_fuel_combustion"
)

exio_combustion_ghg_path: Path = (
    config.exiobase_raw_path
    / "extension"
    / "air_emissions_fuel_combustion"
)

# %% [markdown]
# Finally we specify the folder and file for storing the extracted EXIOBASE data.
# This needs to match the EXIOBASE folder in the main SNAC modelling script, the
# same file name must be specified there.

# %%
exio_interim_path: Path = (
    config.data_path
    / "02_interim" 
    / "exiobase"
) 
exio_interim_path.mkdir(exist_ok=True, parents=True)
exio_extract_file: Path = exio_interim_path / "exio_no_bp_raw.xlsx"

# %% [markdown]
# # Script
# %% [markdown]
# Download EXIOBASE if needed.
if config.download_exiobase:
    for year in config.years_io:
        print(f"Downloading EXIOBASE for year {year}")
        exio_core_path.mkdir(exist_ok=True, parents=True)
        pymrio.download_exiobase3(
            storage_folder=exio_core_path,
            doi=config.exiobase_doi,
            system="ixi",
            years=exio_core_path,
            overwrite_existing=False,
        )

# %% [markdown]
# Load the characterization factors/exiobase stressor conversion.
# %%
exio_stressor_to_nor = pd.read_excel(
    (
        config.data_path
        / "00_auxiliary"
        / "other"
        / "emis_conv_char.xlsx"
    ),
    sheet_name="exio_nor_conv",
    index_col=0
)
charact_ghg = pd.read_excel(
    (
        config.data_path
        / "00_auxiliary"
        / "other"
        / "emis_conv_char.xlsx"
    ),
    sheet_name="nor_char",
    index_col=0
)

# %% [markdown]
# ## Emission and import account generation

# %%
yearly_totals = dict()
yearly_imports = dict()
yearly_diag_imp = dict()
yearly_diag_fd = dict()

for year in config.years_io:
    print(f"Processing EXIOBASE files for year {year}")

    exio_files = exio_core_path.glob(f"ixi/IOT_{year}_ixi*")

    if not exio_files:
        print(f"Warning: No EXIOBASE files found for year {year}.")
        continue
    else:
        # We just take the first found, in case we have zip and extracted folders
        exio_file = next(exio_files)

    exio3 = pymrio.parse_exiobase3(exio_file)

    exio3.reset_extensions()

    # Read combustion and non-combustion emissions accounts
    F_combustion = (
        pd.read_csv(
            (
                exio_combustion_ghg_path
                / "ixi"
                / f"IOT_{year}_ixi"
                / "F.tsv"
            ),
            sep="\t",
            index_col=[0],
            header=[0, 1]
        )
        .rename_axis(exio3.satellite.F.columns.names, axis=1)
        .rename_axis([exio3.satellite.F.index.name], axis=0)
    )
    F_Y_combustion = (
        pd.read_csv(
            (
                exio_combustion_ghg_path
                / "ixi"
                / f"IOT_{year}_ixi"
                / "F.tsv"
            ),
            sep="\t",
            index_col=[0],
            header=[0, 1]
        )
        .rename_axis(exio3.satellite.F_Y.columns.names, axis=1)
        .rename_axis([exio3.satellite.F_Y.index.name], axis=0)
    )
    F_combustion_unit = pd.DataFrame(
        index=F_combustion.index,
        columns=["unit"],
        data="kg"
    )

    # Read combustion and non-combustion emissions accounts
    F_non_combustion = (
        pd.read_csv(
            (
                exio_non_combustion_ghg_path
                / "ixi"
                / f"IOT_{year}_ixi"
                / "F.txt"
            ),
            sep="\t",
            index_col=[0],
            header=[0, 1]
        )
        .rename_axis(exio3.satellite.F.columns.names, axis=1)
        .rename_axis([exio3.satellite.F.index.name], axis=0)
    )
    F_Y_non_combustion = (
        pd.read_csv(
            (
                exio_non_combustion_ghg_path
                / "ixi"
                / f"IOT_{year}_ixi"
                / "F.txt"
            ),
            sep="\t",
            index_col=[0],
            header=[0, 1]
        )
        .rename_axis(exio3.satellite.F_Y.columns.names, axis=1)
        .rename_axis([exio3.satellite.F_Y.index.name], axis=0)
    )
    F_non_combustion_unit = pd.DataFrame(
        index=F_non_combustion.index,
        columns=["unit"],
        data="kg"
    )


    exio3.satellite = pymrio.concate_extension(
        exio3.satellite,
        pymrio.Extension(
            F=F_combustion,
            F_Y=F_Y_combustion,
            unit=F_combustion_unit,
            name="non_combustion"
        ),
        pymrio.Extension(
            F=F_non_combustion,
            F_Y=F_Y_non_combustion,
            unit=F_non_combustion_unit,
            name="non_combustion"
        ),
        name=exio3.satellite.name,
    )

    # HACK: on some systems the concatenate extension sets the index name
    # to "indicator". This happens when the index name of both are not the same.
    # They are the same, so probably some undocumented behaviour in some pandas/python
    # version. To resolve this, we just set the index name again.
    for df in exio3.satellite.get_DataFrame(data=True, with_population=False):
        df.index.name = "stressor"

    # Make a satellite account mirroring the Norwegian emission data
    emis_like_nor = exio3.satellite.characterize(
        factors=exio_stressor_to_nor.reset_index(),
        characterized_name_column="ghg_type",
        characterization_factors_column="factor",
        characterized_unit_column="ghg_type_unit",
        name="emissions",
    )

    emis_like_nor.F.index.name = "stressor"
    emis_like_nor.F_Y.index.name = "stressor"

    emis_ghg = emis_like_nor.characterize(
        factors=charact_ghg.reset_index(),
        characterized_name_column="impact",
        characterization_factors_column="factor",
        characterized_unit_column="impact_unit",
        name="ghg",
    )

    emissions = pymrio.concate_extension(
        emis_ghg,
        emis_like_nor,
        name="emissions"
    )

    exio3.emissions = emissions
    exio3.satellite = None
    exio3.impacts = None

    exio3.calc_all()

    gross_trade = iomath.calc_gross_trade(exio3.Z, exio3.Y)

    imp_emis = (
        exio3.emissions.M
        .mul(
            gross_trade.bilat_flows["NO"],
            axis=1
        )
        .T.groupby("sector", sort=False).sum()
    )

    fd_emis = (
        exio3
        .emissions
        .D_cba.loc[:, "NO"]
        .T.groupby("sector", sort=False).sum()
    )

    imp_monetary = gross_trade.totals.imports["NO"]

    fd_monetary = (
        exio3.Y.loc[:, "NO"]
        .sum(axis=1)
        .T.groupby("sector", sort=False)
        .sum()
    )
    fd_monetary.name = "Final demand"

    # Prep output
    sec = pymrio.get_classification("exio3_ixi").sectors
    rename_dict = (
        sec
        .loc[:, ("ExioName", "ExioCode")]
        .set_index("ExioName")
        .to_dict()["ExioCode"]
    )

    _tot_fd = pd.concat([fd_monetary, fd_emis], axis=1).T
    _tot_fd.loc[:, "unit"] = exio3.emissions.unit
    _tot_fd.loc["Final demand", "unit"] = exio3.unit.iloc[0, 0]
    total_fd = (
        _tot_fd
        .set_index("unit", append=True)
        .T
        .rename(index=rename_dict)
    )

    _tot_imp = pd.concat([imp_monetary, imp_emis], axis=1).T
    _tot_imp.loc[:, "unit"] = exio3.emissions.unit
    _tot_imp.loc["imports", "unit"] = exio3.unit.iloc[0, 0]
    total_imp = (
        _tot_imp
        .set_index("unit", append=True)
        .T
        .rename(index=rename_dict)
    )

    yearly_totals[year] = total_fd
    yearly_imports[year] = total_imp

    # Diagonalization of the emissions
    diag_coll_imp = dict()
    diag_coll_fd = dict()

    for emis in exio3.emissions.get_rows():
        print(f"Diagonalizing {emis}")
        eclean = emis.replace(" ", "_") + "_diag"
        diag_emis = exio3.emissions.diag_stressor(emis, name=emis)
        for df_name in diag_emis.get_DataFrame(
            data=False, with_population=False, with_unit=False
        ):
            df = getattr(diag_emis, df_name)
            setattr(diag_emis, df_name, df.groupby("region").sum())
        diag_emis.unit = diag_emis.unit.groupby("region").first()
        setattr(exio3, eclean, diag_emis)
        exio3.calc_all()

        S_regsum = getattr(exio3, eclean).S

        diag_gross_imports = pd.DataFrame(
            np.diag(gross_trade.bilat_flows.NO),
            index=gross_trade.bilat_flows.NO.index,
            columns=gross_trade.bilat_flows.NO.index,
        )

        fp_src_imp = S_regsum @ exio3.L @ diag_gross_imports

        fp_src_imp_agg = (
            fp_src_imp
            .T.groupby("sector", sort=False).sum()
            .rename(index=rename_dict)
        )

        D_cba = getattr(exio3, eclean).D_cba.NO

        diag_coll_imp[emis] = fp_src_imp_agg
        diag_coll_fd[emis] = D_cba.T.rename(index=rename_dict)

        delattr(exio3, eclean)

    yearly_diag_imp[year] = pd.concat(diag_coll_imp, axis=0, names=["emission"])
    yearly_diag_fd[year] = pd.concat(diag_coll_fd, axis=0, names=["emission"])

# %% [markdown]
# # Saving results
# Save the totals and imports in excel, with one sheet per year/account

# %%
with pd.ExcelWriter(exio_extract_file) as writer:
    for year in config.years_io:
        yearly_totals[year].to_excel(writer, sheet_name=f"fd_{str(year)}")
        yearly_diag_imp[year].to_excel(writer, sheet_name=f"imp_source_{str(year)}")
        yearly_diag_fd[year].to_excel(writer, sheet_name=f"fd_source_{str(year)}")
        yearly_imports[year].to_excel(writer, sheet_name=f"imp_{str(year)}")

# %%
print(f"Done and results stored in {exio_extract_file}.")

# %%
