# -*- coding: utf-8 -*-
# 
from pathlib import Path

# %%
# Paths
## Working directory
try:
    work_path = Path(__file__).parent.absolute() / ".."  # when running as script
except NameError:
    work_path = Path.cwd() / ".."

data_path: Path = work_path / "data"

# Change if you don't have EXIOBASE in local folder.
exiobase_raw_path = (
    data_path
    / "01_raw"
    / "exiobase"
)

# Whether or not to download EXIOBASE
download_exiobase = True
exiobase_doi = "10.5281/zenodo.14614930"


## General input parameters
# The start and end year for the SNAC coupling.
start_year_io: int = 2012
end_year_io: int = 2021

# Year to extract for the SNAC coupling CaFEAN excel_report_format
year_to_extract = 2021

# Year to use for the plotting. 
# 06_cafean_plots.py will only work if 05_SNAC.py has been run
# with the above year_to_extract variable equal to the same same year.
plot_year = 2021

# The start and end year for SUTs needed for the household footprint analysis
start_year_sut: int = 2021
end_year_sut: int = 2021

# IO year for household survey analysis
survey_year: int = 2021

# Base year, substance, and unit of substance for detailed sectoral analysis
base_year = 2012
substance = "GHG"  # Options are "CO2", "CH4", "N2O", "GHG".
substance_unit = "Mt"  # Must match substance unit in result files.

## Derived parameters
years_io = list(range(start_year_io, end_year_io+1))
years_sut = list(range(start_year_sut, end_year_sut+1))


# Norwegian IO variables
final_demand_categories = [
    'Final consumption expenditure by households',
    'Final consumption expenditure by non-profit organisations serving households (NPISH)',
    'Final consumption expenditure by government',
    # 'Final consumption expenditure',
    'Gross fixed capital formation',
    'Changes in valuables',
    'Changes in inventories',
    # 'Changes in inventories and valuables', 
    # 'Gross capital formation',
    # 'Exports intra EU fob (1)',
    # '\nExports fob to members of the euro area (1)',
    # '\nExports fob to non-members of the euro area (1)',
    # 'Exports extra EU fob (1)',
    'Exports fob (2)'
]

# TODO: Consider adding more to extract more information for potential future uses.
agg_final_demand_columns = ["TFUBS", "TUBS"]
value_added_rows = ["RNAM", "RNTS", "RADJ"]
supply_tables_extra_row_codes = ["R", "RNA2", "RN33", "RADJ", "of which:"]
use_tables_extra_row_codes = ["R", "RNA1", "RN33", "RN34", "RADJ", "RZ"]

# Norwegian territorial emissions
emission_types = [
    "CO2",
    "Biomass CO2",
    "N2O",
    "SF6_NF3",
    "CH4",
    "HFC",
    "PFC"
]
