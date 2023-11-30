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
# In the default setting, the EXIOBASE data is downloaded to ./exiobase within
# the project repository. This folder is not tracked by version control (added
# to .gitignore). The script works with both, the original compressed zip
# downloads or extracted folders.
#

# %% [markdown]
# Copyright (C) 2023  XIO Sustainability Analytics, Inc
# 
# Written by
# 
# - Richard Wood
# - Konstantin Stadler
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
import pandas as pd
import numpy as np
import pymrio
import pymrio.tools.iomath as iomath

# %% [markdown]
# Python internal packages (this don't need to be installed, part of standard python)

# %%
from pathlib import Path

# %% [markdown]
# ## Settings

# %% [markdown]
# ### Year parameter
# First, the years for which the EXIOBASE data should be prepared are defined.
# These need to cover the years for the SNAC coupling.

# %%
years: list[int] = list(range(2012, 2021))
# years : list[int] = [2020]


# %% [markdown]
# ### EXIOBASE version
# The doi parameter points the version 3.8.2 (10.5281/zenodo.5589597), you can also
# just download the latest available version by pointing to 10.5281/zenodo.3583070 .
# However, we not recommend this:
#   - we did develop specifically for 3.8.2
#   - by just pointing to the latest version you
#       might not be aware of the actual version used

# %%
exio_doi = "10.5281/zenodo.5589597"

# %% [markdown]
# ### Locations / folder definitions

# %% [markdown]
# First we define the general working directory.

# %%
# Set the work path to the directory where this script is located
# and if this is not available to the current working directory
try:
    work_path = Path(__file__).parent.absolute()  # when running as script
except NameError:
    work_path = Path.cwd()

data_path: Path = work_path / "data"


# %% [markdown]
# Then the file with the characterization factors/exiobase stressor conversion.

# %%
exio_stressor_to_nor = pd.read_excel(
    data_path / "emis_conv_char.xlsx", sheet_name="exio_nor_conv", index_col=0
)
charact_ghg = pd.read_excel(
    data_path / "emis_conv_char.xlsx", sheet_name="nor_char", index_col=0
)


# %% [markdown]
# Here we define the folder for the EXIOBASE raw data.
# This can be anywhere on the disk (in case EXIOBASE is available locally already).
# By default, we set it to ./exiobase within the project repository. This folder
# is not tracked by version control (added to .gitignore).
# The script works with both, the original
# compressed zip downloads or extracted folders.

# %%
exio_raw_path: Path = work_path / "exiobase"

# %% [markdown]
# Then the two files containing the additional biogenic CO2 dataset for EXIOBASE

# %%
F_CO2_bio_file = data_path / "F_ixi_CO2biogenic.csv"
F_Y_CO2_bio_file = data_path / "F_Y_ixi_CO2biogenic.csv"

# %% [markdown]
# Finally we specify the folder and file for storing the extracted EXIOBASE data.
# This needs to match the EXIOBASE folder in the main SNAC modelling script, the
# same file name must be specified there.

# %%
exio_prep_path: Path = data_path
exio_prep_path.mkdir(exist_ok=True)
exio_extract_file: Path = exio_prep_path / "exio_no_bp_raw.xlsx"

# %% [markdown]
# ## EXIOBASE download
# NOTE: Comment this cell if the EXIOBASE files are already downloaded and unpacked.
# Update in pymrio coming soon to handle that case as well.
# If you have downloaded the EXIOBASE files and did not unpack/zip them (as is the case
# when using the script as is), you don’t have to do anything.

# %%
for year in years:
    print(f"Downloading EXIOBASE for year {year}")

    pymrio.download_exiobase3(
        storage_folder=exio_raw_path,
        doi=exio_doi,
        system="ixi",
        years=years,
        overwrite_existing=False,
    )

# %% [markdown]
# ## Emission and import account generation

# %%
yearly_totals = dict()
yearly_imports = dict()
yearly_diag_imp = dict()
yearly_diag_fd = dict()

for year in years:
    print(f"Processing EXIOBASE files for year {year}")

    exio_files = exio_raw_path.glob(f"IOT_{year}_ixi*")

    if not exio_files:
        print(f"Warning: No EXIOBASE files found for year {year}.")
        continue
    else:
        # We just take the first found, in case we have zip and extracted folders
        exio_file = next(exio_files)

    exio3 = pymrio.parse_exiobase3(exio_file)

    exio3.reset_extensions()

    # Add the biogenic CO2 data
    F_CO2_bio = pd.read_csv(F_CO2_bio_file, index_col=0, header=[0, 1], sep="\t").loc[
        year, :
    ]
    F_CO2_bio.name = "CO2 – biogenic - air"
    F_CO2_bio = pd.DataFrame(F_CO2_bio).T
    F_CO2_bio.index.name = exio3.satellite.F.index.name

    F_Y_CO2_bio = pd.read_csv(
        F_Y_CO2_bio_file, index_col=0, header=[0, 1], sep="\t"
    ).loc[year, :]
    F_Y_CO2_bio.name = "CO2 – biogenic - air"
    F_Y_CO2_bio = pd.DataFrame(F_Y_CO2_bio).T
    F_Y_CO2_bio.index.name = exio3.satellite.F_Y.index.name

    F_CO2_bio_unit = pd.DataFrame(index=F_CO2_bio.index, columns=["unit"], data="kg")

    exio3.satellite = pymrio.concate_extension(
        exio3.satellite,
        pymrio.Extension(
            F=F_CO2_bio, F_Y=F_Y_CO2_bio, unit=F_CO2_bio_unit, name="CO2bio"
        ),
        name=exio3.satellite.name,
    )

    # HACK: on some systems the concatenate extension sets the index name
    # to "indicator". This happens when the index name of both are not the same.
    # They are the same, so probably some undocumented behaviour in some pandas/python
    # version. To resolve this, we just set the index name again.
    for df in exio3.satellite.get_DataFrame(data=True, with_population=False):
        df.index.name = "stressor"

    # Make a satellite account mirriroing the Norwegian emission data
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

    emissions = pymrio.concate_extension(emis_ghg, emis_like_nor, name="emissions")

    exio3.emissions = emissions
    exio3.satellite = None
    exio3.impacts = None

    exio3.calc_all()

    gross_trade = iomath.calc_gross_trade(exio3.Z, exio3.Y)

    imp_emis = (
        (exio3.emissions.M * gross_trade.bilat_flows["NO"])
        .T.groupby("sector", sort=False)
        .sum()
    )

    fd_emis = exio3.emissions.D_cba.loc[:, "NO"].T.groupby("sector", sort=False).sum()

    imp_monetary = gross_trade.totals.imports["NO"]

    fd_monetary = exio3.Y.loc[:, "NO"].sum(axis=1).T.groupby("sector", sort=False).sum()
    fd_monetary.name = "Final demand"

    # Prep output
    sec = pymrio.get_classification("exio3_ixi").sectors
    rename_dict = (
        sec.loc[:, ("ExioName", "ExioCode")].set_index("ExioName").to_dict()["ExioCode"]
    )

    _tot_fd = pd.concat([fd_monetary, fd_emis], axis=1).T
    _tot_fd.loc[:, "unit"] = exio3.emissions.unit
    _tot_fd.loc["Final demand", "unit"] = exio3.unit.iloc[0, 0]
    total_fd = _tot_fd.set_index("unit", append=True).T.rename(index=rename_dict)

    _tot_imp = pd.concat([imp_monetary, imp_emis], axis=1).T
    _tot_imp.loc[:, "unit"] = exio3.emissions.unit
    _tot_imp.loc["imports", "unit"] = exio3.unit.iloc[0, 0]
    total_imp = _tot_imp.set_index("unit", append=True).T.rename(index=rename_dict)

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
            fp_src_imp.T.groupby("sector", sort=False).sum().rename(index=rename_dict)
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
    for year in years:
        yearly_totals[year].to_excel(writer, sheet_name=f"fd_{str(year)}")
        yearly_diag_imp[year].to_excel(writer, sheet_name=f"imp_source_{str(year)}")
        yearly_diag_fd[year].to_excel(writer, sheet_name=f"fd_source_{str(year)}")
        yearly_imports[year].to_excel(writer, sheet_name=f"imp_{str(year)}")

# %%
print(f"Done and results stored in {exio_extract_file}.")
