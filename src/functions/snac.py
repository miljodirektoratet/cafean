# -*- coding: utf-8 -*-
# %% [markdown]
#  # Functions used specifically for the SNAC coupling.
#
# %% [markdown]
# Copyright (C) 2024 XIO Sustainability Analytics, Inc
#
# Written by
#
# - Richard Wood
# - Konstantin Stadler
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

# %%
def snac_coupling(Aimp, Ldom, Qimp_mrio):
    """Coupling of MRIO based multipliers with domestic IO tables

    Parameters
    ----------

    Aimp: pd.DataFrame
        Import A matrix
    Ldom: pd.DataFrame
        Domestic L matrix
    Qimp_mrio: pd.DataFrame
        Import multipliers from MRIO

    Returns
    -------
    dict with
        QAimp: pd.DataFrame
            QAimp matrix
        QAimpLdom: pd.DataFrame
            QAimpLdom matrix
    """
    QAimp = Qimp_mrio @ Aimp
    QAimpLdom = QAimp @ Ldom

    return dict(QAimp=QAimp, QAimpLdom=QAimpLdom)


# %%
def snac_findem_breakdown(
    SLdom_diag, QAimp_Ldom_diag, Q_imp_exio_diag, Ydom_br, Yimp_br
):
    """Breakdown of the total fp into final demand categories

    Components: domestic_footprint, import_footprint, total

    Parameters
    ----------
    SLdom_diag: pd.DataFrame
        Diagonalized SLdom matrix for one stressor/gas/impact
    QAimp_Ldom_diag: pd.DataFrame
        Diagonalized QAimpLdom matrix for one stressor/gas/impact
    Q_imp_exio_diag: pd.DataFrame
        Diagonalized Q_imp_exio matrix for one stressor/gas/impact
    Ydom_br: pd.DataFrame
        domestic_footprint final demand matrix without aggregates
    Yimp_br: pd.DataFrame
        Import final demand matrix without aggregates

    Returns
    --------
    pd.Series

        index: pd.MultiIndex

            level 0: component
                domestic_footprint, import_footprint, total_footprint
            level 1: sector
                sector names
            level 2: category
                final demand category


    """

    dom_breakdown = SLdom_diag @ Ydom_br
    imp_breakdown = QAimp_Ldom_diag @ Ydom_br + Q_imp_exio_diag @ Yimp_br
    tot_footprint_breakdown = dom_breakdown + imp_breakdown

    index_names = ["component"] + list(dom_breakdown.index.names)

    ds = pd.concat(
        [dom_breakdown, imp_breakdown, tot_footprint_breakdown],
        keys=["domestic_footprint", "import_footprint", "total_footprint"],
        names=index_names,
        axis=0,
    ).stack()

    return ds


# %%
def calc_snac_sector_accounts(
    SLdom,
    QAimp_Ldom,
    Ydom_diag,
    Qimp_mrio,
    Yimp_diag,
    Ydom_exp_diag,
    Yimp_exp_diag,
    total_imports,
    pba,
):
    """Calculate and merge sector accounts

    Parameters
    ----------
    SLdom: pd.DataFrame

    QAimp_Ldom: pd.DataFrame

    Qimp_mrio: pd.DataFrame

    Ydom_diag: pd.DataFrame

    Yimp_diag: pd.DataFrame

    Ydom_exp_diag: pd.DataFrame

    Yimp_exp_diag: pd.DataFrame

    total_imports: pd.DataFrame

    pba: pd.DataFrame

    Returns
    --------
    pd.Series

        index: pd.MultiIndex

            level 0: component
                domestic_footprint, import_footprint, total_footprint
            level 1: emission
                sector names
            level 2: sector
                sector names
    """

    total_footprint = SLdom @ Ydom_diag + QAimp_Ldom @ Ydom_diag + Qimp_mrio @ Yimp_diag

    domestic_footprint = SLdom * Ydom_diag.sum(axis=1)
    import_footprint = QAimp_Ldom @ Ydom_diag + Qimp_mrio @ Yimp_diag
    export_footprint = (
        SLdom @ Ydom_exp_diag + QAimp_Ldom @ Ydom_exp_diag + Qimp_mrio @ Yimp_exp_diag
    )

    footprint_gross_imports = Qimp_mrio * total_imports
    footprint_gross_imports.columns.name = "sector"

    index_names = ["component", "emission"]

    ds = pd.concat(
        [
            total_footprint,
            domestic_footprint,
            import_footprint,
            export_footprint,
            footprint_gross_imports,
            pba,
        ],
        keys=[
            "total_footprint",
            "domestic_footprint",
            "import_footprint",
            "footprint_gross_exports",
            "footprint_gross_imports",
            "production_account",
        ],
        names=index_names,
        axis=0,
    ).stack()

    return ds



# %%
def exio_nor_converter(to_conv, sector_matching, new_col_names="name"):
    """Converts classification from EXIOBASE to Norwegian sectors

    This checks if the columns of to_conv matches on of
    the index of sector_matching.

    Parameters
    ----------
    to_conv: pd.DataFrame
        Accounts to convert to Norwegian sectors.
        EXIOBASE classification in the columns.

    sector_matching: pd.DataFrame
        Matching between EXIOBASE and Norwegian sectors.
        Read from sector matching file.

    new_col_names: str
        Can be 'name', 'code' or 'nr'.
        Taken from sector_matching columns and used as the new
        sector names.

    Returns
    -------

    pd.DataFrame
        Accounts converted to Norwegian sectors
    """
    # As we do an concordance based conversion, we check for the correct order
    # for at least on of the index levels of sector_matching
    for ind_number in range(0, 100):
        try:
            sec_match_order = sector_matching.index.get_level_values(ind_number)
            if sec_match_order.equals(to_conv.columns):
                break
        except IndexError:
            raise ValueError("Could not find matching index value do not match")

    ret = to_conv @ sector_matching.values
    ret.columns = sector_matching.columns.get_level_values(new_col_names).str.strip()
    return ret

def stack_and_unit(dd, emission_unit_mapper, sector_order):
    """Stacks dataframes for final save and assign units

    The results are ordered
        - emission: order given in emis_units
        - sectors: order given in sector_order
        - years: ascending

    Parameters
    ----------
    dd : dict
        With years as keys and data-frames as values
    emis_units : pd.Series
        Emission units with index of emission names.
        The order of the index is used for the order of
        the result.
    sector_order : list
        Orders list of sectors for reordering the results

    Returns
    -------
    pd.DataFrame
        index: Multiindex with year and the remaining ones from the dd values
        columns: value, unit

    """
    if type(dd) is dict:
        merged = (
            pd.concat(dd, axis=0, names=["year"])
            .to_frame("value")
        )
    else:
        merged = dd.copy()

    merged["unit"] = merged.index.get_level_values("emission").map(emission_unit_mapper)
    merged = (
        merged
        .sort_index(level="year", sort_remaining=False)
        .reindex(list(emission_unit_mapper.keys()), level="emission")
    )
    if "sector" in merged.index.names:
        merged = merged.reindex(sector_order, level="sector")

    return merged