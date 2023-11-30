# -*- coding: utf-8 -*-
"""
Visualization app for the results.

This script provides an interactive way to explore the results of the CaFEAN project.
To run the app, use the command `streamlit run vis_app.py` in the terminal.


Copyright (C) 2023  XIO Sustainability Analytics, Inc

Written by

- Richard Wood
- Konstantin Stadler

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

from pathlib import Path

import pandas as pd
import altair as alt
import streamlit as st

# %%
# Set the work path to the directory where this script is located
try:
    work_path = Path(__file__).parent.absolute()  # when running as script
except NameError:
    raise NameError("Please run this script with 'streamlit run vis_app.py'")

RESULTS_PATH = work_path / "results"
if not RESULTS_PATH.exists():
    raise FileNotFoundError(
        f"Results folder not found at {RESULT_PATH}. Please run the aggregator script first."
    )

EMISSION_ORDER = ["GHG", "CO2", "CH4", "HFC", "N2O", "PFC", "SF6_NF3", "Biomass CO2"]
REQUIRED_COLS = ["year", "value", "unit", "emission"]

def generic_line_chart(df):
    """ The standard line chart for most plots.

    Works with a preselected dataframe.
    """
    if df.emission.nunique() > 1:
        st.write("Please select only one emission to print graph.")
        return
    else:
        selected_emission = df.emission.unique()[0]
        unit = df.unit.unique()[0]
    non_unique_col = []
    for col in df.columns:
        if col in REQUIRED_COLS:
            continue
        if df[col].nunique() > 1:
            non_unique_col.append(col)
    if len(non_unique_col) > 1:
        st.write(
            f"Only one of the following columns can have multiple selections: {non_unique_col}"
        )
        return
    elif len(non_unique_col) == 0:
        color_graph = "emission"
    else:
        color_graph = non_unique_col[0]
        # nr_legend_col = int(len(df[color_graph].unique()) / 10) + 1  # for vertical legends
        if int(len(df[color_graph].unique()) < 20):
            nr_legend_col = int(len(df[color_graph].unique()) / 3) + 1
        else: 
            # guard for many sector in the sector graph
            nr_legend_col = int(len(df[color_graph].unique()) / 20) + 1

    common_x = alt.X("year", 
                     axis=alt.Axis(format=".0f",
                                   tickCount=len(df.year.unique())))
    common_y = alt.Y("value", title=f"{selected_emission} [{unit}]")
    common_tooltip = [
        alt.Tooltip("value", title=f"{selected_emission} in {unit}", format=".2f"),
        alt.Tooltip("year", title="Year"),
        alt.Tooltip(f"{color_graph}", title=f"{color_graph.title()}"),
    ]
    selection = alt.selection_multi(fields=[color_graph], bind="legend")

    line = (
        alt.Chart(df)
        .mark_line()
        .encode(
            x=common_x,
            y=common_y,
            color=color_graph,
            opacity=alt.condition(selection, alt.value(1.0), alt.value(0.1)),
        )
        .add_selection(selection)
    )

    dots = (
        alt.Chart(df)
        .mark_circle(size=200)
        .encode(
            x=common_x,
            y=common_y,
            color=alt.Color(color_graph, title=f"{color_graph.title()}"),
            tooltip=common_tooltip,
            opacity=alt.condition(selection, alt.value(1.0), alt.value(0)),
        )
        .add_selection(selection)
    )

    chart = (
        (line + dots)
        .properties(height=900)
        .interactive()
        .configure_axis(labelFontSize=20, titleFontSize=25)
        .configure_legend(
            labelLimit=400,
            symbolLimit=0,
            orient="bottom",
            titleFontSize=25,
            labelFontSize=20,
            labelColor="silver",
            columns=nr_legend_col,
        )
    )

    st.altair_chart(chart, use_container_width=True)
    st.divider()
    st.write("Tip: (Shift) + Click on the legend to highlight one graph-line.")


def show_home():
    st.header("CafEAN - Result Visualization", divider="rainbow")
    st.subheader("Interactive visualization demo")

    st.markdown(
        """
    The Visualization App offers an interactive way to 
    explore the results of the CaFEAN project.
   
    The CaFEAN project linked
    the official Norwegian (SSB) IO and greenhouse gas emissions data 
    with the Multi-Regional Input-Output (MRIO) model EXIOBASE
    in order to produce carbon footprint results for Norway. 
    These results are consistent with the official Norwegian statistics.

    To use the app, select the perspective you want to visualize in the sidebar.

    Most perspectives depict the aggregated results of the CaFEAN model.
    However, the *Result File Explorer* perspective allows you to upload 
    any tsv file in long format and visualize it.

    See the perspective 'help' for more details.
    """
    )

    st.divider()

    st.markdown(
        """
    Developed by 

    - Vector Sustainability Pty Ltd 
    - XIO Sustainability Analytics A/S 

    for the Norwegian Environment Protection Agency."""
    )


def show_top_level_results():
    st.header("Top Level Results")

    df = pd.read_csv(
        RESULTS_PATH / Path("total_accounts.tsv"), sep="\t"
    )  # Replace "your_data.csv" with the actual file path

    absolute, distribution = st.tabs(["Absolute", "Distribution"])
    st.session_state.top_sidebar_select_disabled = True

    with st.sidebar:
        selected_emission = st.sidebar.selectbox("Emission", EMISSION_ORDER)

    to_show = df[(df["emission"] == selected_emission)]

    with absolute:
        generic_line_chart(df=to_show)

    with distribution:
        min_y = df.year.min() - 0.5
        max_y = df.year.max() + 0.5

        to_show_total_stack = to_show[
            to_show["component"].isin(
                ["domestic_footprint", "import_footprint", "household_totals"]
            )
        ]
        unit = to_show_total_stack.unit.unique()[0]
        common_x = alt.X(
            "year",
            axis=alt.Axis(
                format=".0f", tickCount=len(to_show_total_stack.year.unique())
            ),
        ).scale(domain=(min_y, max_y))
        common_y = alt.Y("value", title=f"{selected_emission} [% of total]").stack("normalize")
        common_tooltip = [
            alt.Tooltip("value", title=f"{selected_emission} in {unit}", format=".2f"),
            alt.Tooltip("year", title="Year"),
            alt.Tooltip("component", title="Component"),
        ]

        bar = (
            alt.Chart(to_show_total_stack)
            .mark_bar(size=50)
            .encode(
                x=common_x,
                y=common_y,
                color="component",
                tooltip=common_tooltip,
            )
        )
        chart = (
            (bar)
            .properties(height=900)
            .configure_axis(labelFontSize=20, titleFontSize=25)
            .interactive()
            .configure_legend(
                labelLimit=400,
                symbolLimit=800,
                titleFontSize=25,
                labelFontSize=20,
                labelColor="silver",
        )

        )

        st.altair_chart(chart, use_container_width=True)


def show_sources():
    st.header("Source of Emissions")

    df_detail = pd.read_csv(
        RESULTS_PATH / Path("footprint_sources_agg.tsv"), sep="\t"
    )  # Replace "your_data.csv" with the actual file path
    df_tot = pd.read_csv(
        RESULTS_PATH / Path("footprint_sources_totals_agg.tsv"), sep="\t"
    )  # Replace "your_data.csv" with the actual file path
    df_tot.loc[:, "sector"] = "TOTAL"

    df = pd.concat([df_tot, df_detail])

    with st.sidebar:
        selected_emission = st.sidebar.selectbox("Emission", EMISSION_ORDER)
        selected_sector = st.sidebar.selectbox("Sector", df["sector"].unique())
        to_show = df[
            (df["emission"] == selected_emission) & (df["sector"] == selected_sector)
        ]
    generic_line_chart(df=to_show)


def show_final_demand_categories():
    st.header("Final Demand Categories Details")

    df = pd.read_csv(
        RESULTS_PATH / Path("footprints_final_demand_breakdown_agg.tsv"), sep="\t"
    )

    with st.sidebar:
        selected_emission = st.sidebar.selectbox("Emission", EMISSION_ORDER)
        selected_sector = st.sidebar.selectbox("Sector", df["sector"].unique())
        selected_component = st.sidebar.selectbox("Component", df["component"].unique())
        to_show = df[
            (df["emission"] == selected_emission)
            & (df["sector"] == selected_sector)
            & (df["component"] == selected_component)
        ]
    generic_line_chart(df=to_show)


def show_details():
    st.header("Generic Result File Explorer")

    with st.sidebar:
        uploaded_file = st.file_uploader("Upload CSV file", type="tsv")

        col_to_select = []
        if uploaded_file is not None:
            df = pd.read_csv(uploaded_file, sep="\t")

            if not all([col in df.columns for col in REQUIRED_COLS]):
                st.write(f"File most contain columns: {REQUIRED_COLS}")
                return

            for col in df.columns:
                if col == "value":
                    df[col] = df[col].astype(float)
                    continue
                elif col == "unit":
                    continue
                if col == "year":
                    df[col] = df[col].astype(int)
                col_to_select.append(col)

            selection = dict()
            for col in col_to_select:
                if col == "year":
                    year_select = st.sidebar.slider(
                        "Years",
                        min_value=df[col].min(),
                        max_value=df[col].max(),
                        value=(df[col].min(), df[col].max()),
                        step=1,
                    )

                    selection[col] = range(year_select[0], year_select[1] + 1)
                else:
                    selection[col] = st.multiselect(f"Select {col}", df[col].unique())
            for col, values in selection.items():
                if len(values) > 0:
                    df = df[df.loc[:, col].isin(values)]

    graph_tab, table_tab = st.tabs(["Graph", "Table"])

    def write_empty():
        st.write("Please select file to show.")
        st.markdown(
            """ 
                The detail explorer works with any tsv file in long format.\n

                In requires the columns "value" to be present and the file
                separator to be a tab.

                It works "out of the box" with the results generated 
                by the *nor_exio_snac.py* script.
                """
        )

    with table_tab:
        try:
            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "value": st.column_config.NumberColumn(format="%.2f"),
                    "year": st.column_config.NumberColumn(format="%.0f"),
                },
            )
        except UnboundLocalError:
            write_empty()

    with graph_tab:
        try:
            generic_line_chart(df)
        except UnboundLocalError:
            write_empty()


def show_help():
    st.markdown(
        """
    ## Variable Descriptions

    ### Accounts in the Top Level Results

    - **total_footprint**: Total footprint, sum of domestic and import footprint
    - **domestic_footprint**: domestic component of the footprint
    - **import_footprint**: import component of the footprint
    - **footprint_gross_imports**: embodied emissions in imported goods and services as they enter the Norwegian border
    - **footprint_gross_exports**: embodied emissions in exports as they leave the Norwegian border
    - **production_account**: production based account, official Norwegian statistics
    - **household_totals**: total direct household emissions (same as Household totals in household_emissions.tsv / official statistics)
    - **total_footprint_with_households**: total footprint including household emissions 
    - **production_account_with_households**: total production based account including household emissions 

    ### Sources

    These show the sources (in which country the emissions actually occur) of the Norwegian footprint.

    ### Final Demand Categories

    These show the footprint split into the official Norwegian final demand categories.

    ### Result File Explorer

    This perspective allows you to upload any result tsv file in long format and visualize it.

    """
    )


def main():
    st.set_page_config(
        page_title="CaFEAN Vis App",
        page_icon="📈",
        layout="wide",
        initial_sidebar_state="expanded",
        menu_items={
            "About": """ Developed by _Vector Sustainability Pty Ltd_ and
                        _XIO Sustainability Analytics A/S_ 
                        for the Norwegian Environment Protection Agency."""
        },
    )

    st.markdown(
        r"""
        <style>
        .stDeployButton {
                visibility: hidden;
            }
        </style>
        """, unsafe_allow_html=True
    )
    # Create sub-pages
    subpages = [
        "Home",
        "Top Level Results",
        "Sources",
        "Final Demand Categories",
        "Result File Explorer",
        "Help",
    ]
    with st.sidebar:
        st.header("CafEAN - Vis App", divider="rainbow")
        subpage = st.sidebar.selectbox("Perspective: ", subpages)

    if subpage == "Home":
        show_home()
    elif subpage == "Top Level Results":
        show_top_level_results()
    elif subpage == "Sources":
        show_sources()
    elif subpage == "Final Demand Categories":
        show_final_demand_categories()
    elif subpage == "Result File Explorer":
        show_details()
    elif subpage == "Help":
        show_help()


if __name__ == "__main__":
    main()
