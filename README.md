# This is not up to date!
# CaFEAN: Code for the coupled footprint/simplified-SNAC model
This repository contains the Python code for the CaFEAN model, linking the Multi-Regional Input-Output (MRIO) model EXIOBASE with the official Norwegian (SSB) IO and greenhouse gas emissions data in order to produce carbon footprint results for Norway.

The README is structured as follows:

- [Repository Structure](#repository-structure)
    - [Project folder](#project-folder)
    - [Code](#code)
    - [Documents](#documents)
    - [Data](#data)
    - [CaFEAN results](#cafean-results)
    - [Detailed sector analysis results](#detailed-sector-analysis-results)
    - [Household survey footprint results](#household-survey-footprint-results)
- [Running the code](#running-the-code)
    - [Setting up the environment](#setting-up-the-environment)
    - [Use cases](#use-cases)
- [Development guide](#development-guide)

The CaFEAN model is licensed under the AGPLv3 license. See the [LICENSE](./LICENCE.md) file for details.

Copyright (C) <2024>  <XIO Sustainability Analytics, Inc>

## Repository Structure

### Project folder

The folder contains the following files:

- `README.md`: This file.
- `LICENSE.md`: License which is applicable to this repository.
- `environment.yml`: Environment file used in Anaconda (or Mamba) to install the Python packages required to run the code.
- `.gitignore`: Git ignore file.

and the following folders:

- `src/`: All code used in the project.
- `data/`: All the data used in the project.
- `docs/`: Various documents relevant for the CaFEAN project.

<br>

### Code

The `src/` folder contains a config file, 3 function files, and 9 scripts. All script files can be run in interactive mode in Python.

The purpose of each file is:

- `./config.py`: This file contains various configurations options that are used in the scripts. User needs to ensure to keep this up to date.

<br>

- `./functions/general.py`: Provides functions that are used across the scripts.
- `./functions/gras.py`: Provides an implementation of the GRAS algorithm.
- `./functions/snac.py`: Provides functions used for the SNAC coupling.

<br>

- `./01_prepare_exiobase.py`: This script prepares EXIOBASE 3.9 data for the use in the SNAC coupling
- `./02_prepare_norwegian_IO.py`: This script prepares the Norwegian IO data for general use.
- `./03_prepare_deflators.py`: This script prepares the defaltors for deomstic production and imports from two SSB datasets.
- `./04_prepare_emissions.py`: This script prepares the Norwegian territorial emissions dataset as reported by SSB for the SNAC coupling.
- `./05_SNAC.py`: This script calcualted the Norwegian consumption-based emisions accounts based on the hybrid-SNAC approach.
- `./06_cafean_plots.py`: This script all figures and tables shown in the original CaFEAN report.
- `./07_detailed_sectoral_analysis.py`: This script makes the results files used in the detailed sector analysis.
- `./08_prepare_household_surveys.py`: This script disaggregates the household footprint calculated by the SNAC model using SSB's consumer expenditure surveys.
- `./09_household_survey_analysis.py`: This script makes the tables and figures used for the household footprint analysis.

<br>

### Documents

The `docs/` folder is available to upload any relevant documetation, such as the CaFEAN report and its appendices.

<br>

### Data

The `data/` folder contains all the data required to run the code. It is structured as follows:
```
├───00_auxiliary
│   ├───classifications
│   │       region_agg_spec.tsv         # Aggregation used for regions in the report.
│   │       sector_agg_spec.tsv         # Aggregation used for sectors in the report.
│   │       sector_list.xlsx            # Sector codes and names.
│   │       survey_categories.xlsx      # Aggregation used to match household statistics with consumer surveys.
│   ├───concordances
│   │       deflators_and_io.xlsx       # Maps between deflator datasets and the Norwegian IO accounts.
│   │       exiobase_and_io.xlsx        # Maps between EXIOBASE and the Norwegian IO accounts.
│   │       survey_and_io.xlsx          # Maps between the consumer surveys and the Norwegian IO accounts.
│   └───other
│           emis_conv_char.xlsx         # Contains conversions factors for GHG emissions.
│           exchange_rates.csv          # EUROSTAT exchange rates (https://doi.org/10.2908/TEC00033)
├───01_raw
│   ├───exiobase                        # EXIOBASE 3.9.4 data which is not included in the repository. 
│   │   ├───core
│   │   │   └───ixi
│   │   │       ├───IOT_1995_ixi
│   │   │       │   ├───impacts
│   │   │       │   └───satellite
│   │   │       └─── ...
│   │   └───extension
│   │       ├───air_emissions_fuel_combustion
│   │       │   └───ixi
│   │       │       ├───IOT_1995_ixi
│   │       │       └───...
│   │       └───air_emissions_non_fuel_combustion
│   │           └───ixi
│   │               ├───IOT_1995_ixi
│   │               └───...
│   └───ssb                             # SSB data as downloaded from their website.
│       ├───deflators
│       │       07337.xlsx
│       │       09170_A64.xlsx
│       ├───household consumption surveys
│       │       14100.xlsx
│       │       14152.xlsx
│       │       14156.xlsx
│       │       14157.xlsx
│       │       14161.xlsx
│       ├───household statistics
│       │       06070.xlsx
│       │       07459.xlsx
│       │       10986.xlsx
│       │       12563.xlsx
│       ├───norwegian IO
│       │   ├───domestic
│       │   │       iot_1850_[year].xlsx
│       │   ├───import
│       │   │       iot_1950_[year].xlsx
│       │   ├───supply tables
│       │   │       iot_1500_2021.xlsx
│       │   └───use tables
│       │           iot_1600_2021.xlsx
│       └───territorial emissions
│             AEA Questionnaire_2023_Norway2023_til footprint.xlsm
├───02_interim
│   ├───deflators                       # Output of script 03_prepare_deflators.py
│   ├───exiobase                        # Output of script 01_prepare_exiobase.py
│   ├───household survey footprints     
│   │   └───2021                        # Output of 08_prepare_household_surveys.py - contains useful intermediate files.      
│   ├───norwegian IO                    # Output of 02_prepare_norwegian_IO.py - useful for any future IO analysis.
│   │   ├───2012
│   │   └───...
│   └───territorial emissions           # Output of 04_prepare_emissions.py - Norwegian emissions data in a better format.
├───03_results
│   ├───CaFEAN                          # Output made for original CaFEAN deliverable
│   │   └───excel_report_format         # Output in Excel format used for original plots.
│   ├───Detailed analysis               # Output of 07_detailed_sectoral analysis.
│   └───Household survey footprints     # Output of 09_household_survey_analysis.py
│       └───2021
├───04_plots
│   ├───CaFEAN                          # Output of 06_cafean_plots.py
│   └───Household survey footprints     # Output of 09_household_survey analysis.
│       └───2021
└───legacy
    └───CaFEAN results                  # Original CaFEAN results.
        ├───excel_report_format
        └───figures
```

<br>

###  CaFEAN results

Result files are stored in a long table format. This can easily be converted to
any data format, including SQL. In excel, the list can easily be pivoted to a
wide format or filtered for specific years. In Python/Pandas, this can be
extended to a wide format with .unstack(index_name). The table also contains
the unit for each data entry. For further reprocessing in Python this might be
cumbersome and can be dropped with .drop('unit', axis=1). 

The files are:

- `./data/03_results/CaFEAN/sector_accounts.tsv`: Sector accounts with the following columns:
    - year
    - component: 
        - total_footprint: total footprint of a sector
        - domestic_footprint: domestic component of the footprint
        - import_footprint: import component of the footprint
        - footprint_gross_imports: embodied emissions in imported goods and services as they enter the Norwegian border
        - footprint_gross_exports: embodied emissions in exports as they leave the Norwegian border
        - production_account: production based account, official Norwegian statistics
    - emission (GHG, CO2, etc)
    - sector: Norewegian sector classification
    - value
    - unit

- `./data/03_results/CaFEAN/footprints_final_demand_breakdown.tsv`: Final demand breakdown of the footprint accounts with the following columns:
    - year
    - component:
        - domestic_footprint: domestic component of the footprint
        - import_footprint: import component of the footprint
        - total_footprint: total footprint of a sector
    - category (final demand category, including exports)
    - emission (GHG, CO2, etc)
    - sector
    - value
    - unit

- `./data/03_results/CaFEAN/household_emissions.tsv`: Household emissions (official Norwegian data) with the following columns:

    - year
    - emission
    - hhld_component (Household totals, Transport, Heating/cooling, Other, as per official Norwegian statistics)
    - value
    - unit

- `./data/03_results/CaFEAN/total_accounts.tsv`: These are total sums of the sector_accounts, also including household emissions. The columns are:
    - year
    - component: 
        - totals of all components of sector_accounts.tsv 
        - household_totals: total direct household emissions (same as Household totals in household_emissions.tsv / official statistics)
        - total_footprint_with_households: total footprint including household emissions 
        - production_account_with_households: total production based account including household emissions 
    - emission
    - value
    - unit

- `./data/03_results/CaFEAN/footprint_sources.tsv`: Source of imports
    - year
    - emission
    - region (region where the emissions occur/imported from)
    - sector
    - value
    - unit

- `./data/03_results/CaFEAN/footprint_sources_totals.tsv`: Total source of imports (all sectors aggregated)
    - year
    - emission
    - region (region where the emissions occur/imported from)
    - value
    - unit

- `./data/03_results/CaFEAN/footprint_sources_totals_with_households.tsv`: Total source of imports incl. households emissions (all sectors aggregated)
    - year
    - emission
    - region (region where the emissions occur/imported from)
    - value
    - unit

- `./data/03_results/CaFEAN/full_household_footprint.tsv`: Detailed footprint of households.
    - year
    - component
        - domestic: emissions occur in Norway. 
        - import_domestic: emissions embodied embodied in imported goods/services that are first used by Norwegian producers before delivering final good to Norwegian households.
        - import: emissions embodied in imported goods/services that are consumed directly by households.
    - emission (GHG gas type.)
    - sector (final good/service consumed by households)
    - region (region where the emissions occur/imported from)
    - value
    - unit

- `./data/03_results/CaFEAN/excel_report_format/` : Wide table format for excel import of the data underlying the figures of the report. 

Besides these data files, the script also creates aggregated versions
of these files, indicated by the suffix "_agg". These files have the same 
structure as explained above, with sectors/regions aggregated.

<br>

###  Detailed sector analysis results

The `sectoral_analysis.xlsx` file contains the following columns:
- <b>Year</b>: The year at which the economic activity and emission of greenhouse gas takes place. Ranges from 2012 to 2021. 
- <b>Variable</b>: Whether the emission occurs domestically (Domestic) or if the emission is embodied in imports (Import).
- <b>Final sector</b>: The sector that delivers the good or service to final demand. Either at the full 64 sector IO classification or the aggregate 17 sector classification (version with suffix “_agg.tsv”).
- <b>Source sector</b>: The sector that causes the emissions in the supply chain of the final sector. Provided at the full 64 sector classification. 
- <b>Footprint</b>: The amount of GHG emissions emitted in source sector to satisfy the demand on the good or service provided by final sector. In megatons of CO2 equivalent.
- <b>Source sector footprint share</b>: The source sector’s share (in percentages) of the total footprint of Final sector.
- <b>Absolute difference in source sector footprint share from 2012</b>: Difference in Source sector footprint share from 2012 to year.
- <b>Absolute difference in footprint from 2012</b>: Absolute difference in footprint from 2012 to year occurring in source sector due to final demand on final sector. In megatons of CO2 equivalent.
- <b>Percentage change in footprint from 201</b>: Percentage change in footprint from 2012 to year in source sector due to final demand on final sector.
- <b>Percentage change in stressor (current prices) from 2012</b>: Percentage change in non-deflated stressor for source sector from 2012 to year. 
- <b>Percentage change in stressor (constant prices) from 2012</b>: Percentage change in deflated stressor for source sector from 2012 to year.
- <b>Source sector: F</b>: Total emissions in source sector in year. In megatons of CO2 equivalent.
- <b>Source sector: x (current)</b>: Total output of source sector in year. In millions of Norwegian kroner. 
- <b>Source sector: x (constant)</b>: Total output of source sector in year. In millions of Norwegian kroner in 2015 prices.
- <b>Percentage change in source sector: F</b>: Percentage change in total emissions in source sector from year 2012 to year. 
- <b>Percentage change in source sector: x (current)</b>: Percentage change in total output of source sector from year 2012 to year.
- <b>Percentage change in source sector: x (constant)</b>: Percentage change in total output of source sector from year 2012 to year.

The `sectoral_analysis_agg.xlsx` has the same structure but with final sector aggregated. 

<br>

### Household survey footprint results

A table is provided for each survey with survey product level along the rows and two-level column levels:

1.	Variable: Total, per household, or per person.
2.	Survey category: Centrality, household type, or income quartile.

The cells then contain the footprint. 
The total columns are in megatons CO2-eqv., while the others are in tons CO2-eqv. 

Aggregated product group versions (i.e., COICOP level one) of the tables are also provided for per person/household results.
These tables have formatted background to highlight large/small values.

Intermediate results of the household survey, which is likely to be of interest, can be found in `/data/02_interim/household survey footprints/[year]/`.

<br>

## Running the code

The code is setup to run with version 3.9.4 of EXIOBASE. 

Before running the scripts, one needs to: 

- [Build the Anaconda environment (first time)](#setting-up-the-environment) - only necessary the first time 
- Activate the environment:

        conda activate cafean

- Configure the `./src/config.py` file

- Run the 9 scripts in sequence as indicated by their prefix. 

Then one can run the 9 scripts in the sequence indicated by their prefix.

It is **not necessary to run the full sequence**. See [section use cases below](#use-cases) for some standard work flows you might encounter.

<br>

### Setting up the environment

We recommend to run the code within a designated environment.
During the development we used the [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) specified in the file *environment.yml*.

To build this environment, you need to have [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) installed. 
We actually recommend mamba, because of the faster environment building.
However, mamba is drop-in replacement for conda, and both work interchangeably.

For mamba specifically on Windows you can [download the installer here.](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Windows-x86_64.exe).
For all other OS follow the link above to the installation instructions.

If you leave all at the default settings (but please make sure to not have any special characters - e.g. Norwegian letters - in the install path), you should get a "Miniforge Prompt" in the start menu.

After opening the prompt, navigate to the folder you have downloaded/cloned the code to.
You can then install the environment by running

        mamba env create -f environment.yml

(replace "mamba" with "conda" if you use conda instead of mamba).

This might take some time! It only need to run once/when the environment changes.

Afterwards run

     mamba activate cafean

and then

     jupyter lab

If that works, you have successfully establish a coding environment called "cafean" (the create command from above), activated it and then started Jupyter lab within the environment.

<br>

### Use cases

Here we give some example use cases and how to run the code for them. We go from simple cases to more complex ones.

**WARNING:** Running the code overwrite previous results! Make sure to backup the results before running the code again. 
Alternatively, all scripts define the output folder in the first lines of the script.

<br>

#### Changing the config.py file

The `config.py` can be changed to adjust the scripts for run for future years when data becomes available. 
Comments are provided in the file that explains what the purpose of each relevant variable. 

<br>

#### Changing the aggregation of the results

For the report we aggregated the results to 17 sectors and 4 regions.
The aggregations are defined in `./data/00_auxiliary/classifications/region_agg_spec.tsv` and `./data/00_auxiliary/classifications/sector_agg_spec.tsv`.

To change the aggregation, simply change the aggregation specification in these files and rerun.
* `05_SNAC.py` for new datasets.
* `06_cafean_plots.py` for new plots.
* `07_detailed_sectoral_analysis` for new detailed sector analysis.

<br>

#### Changing the exchange rates

Exchange rates are given in the file `./data/00_auxiliary/other/exchange_rates.csv`. 
These are average Eurostat exchange rates. In case other data sources for these should be used, just change them in this file.
Then rerun: 

* `05_SNAC.py` for new datasets.
* `06_cafean_plots.py` for new plots.
* `07_detailed_sectoral_analysis.py` for new detailed sector analysis.
* `08_prepare_household_surveys.py` for new (intermediate) household survey results.
* `09_household_survey_analysis.py` for new household survey results.

<br>

#### Updating the Norwegian Emission or IO data

As long as the updated data keeps the same format/excel structure this can be done by simply replacing the old data with the new one in the `./data/01_raw/ssb/territorial emissions/` folder. 
Then rerun: 

* `04_prepare_emissions.py` for new intermediate emissions datasets.
* `05_SNAC.py` for new datasets.
* `06_cafean_plots.py` for new plots.
* `07_detailed_sectoral_analysis.py` for new detailed sector analysis.
* `08_prepare_household_surveys.py` for new (intermediate) household survey results.
* `09_household_survey_analysis.py` for new household survey results.

In case the structure of the data changed, the `04_prepare_emissions.py` script needs to be adapted accordingly.

<br>

#### Updating the GHG characterisation factors 

The characterisation factor from GHG emissions to GHG totals are based on the GWP100 factors of the IPCC Fifth Assessment Report.
These are given in the file `./data/00_auxiliary/other/emis_conv_char.xlsx`. The files contains a sheet "meta" given some information on the structure.
In short, the characterisation factors are given in the sheet "nor_char" and can be updated there.

Then rerun: 

* `01_prepare_exiobase.py` for new intermediate EXIOBASE datasets. 
* `04_prepare_emissions.py` for new intermediate emissions datasets.
* `05_SNAC.py` for new datasets.
* `06_cafean_plots.py` for new plots.
* `07_detailed_sectoral_analysis.py` for new detailed sector analysis.
* `08_prepare_household_surveys.py` for new (intermediate) household survey results.
* `09_household_survey_analysis.py` for new household survey results.

Running the `01_prepare_exiobase.py` can take some time and use a decent amount of RAM memory (around 20 GB). See [Preparing the EXIOBASE raw data](#preparing-the-exiobase-raw-data) for more details.

<br>

#### Running for a new EXIOBASE version

The code is setup to run with version 3.9.4 of EXIOBASE which at the time of writing is not publically available on Zenodo.
Hence the data needs to downloaded manually for now. 
Once new versions of EXIOBASE data is available on Zenodo, one needs to update the `download_exiobase` and `exiobase_doi` variables in `src/config.py`, and then run:

* `01_prepare_exiobase.py` for new intermediate EXIOBASE datasets. 
* `05_SNAC.py` for new datasets.
* `06_cafean_plots.py` for new plots.
* `07_detailed_sectoral_analysis.py` for new detailed sector analysis.
* `08_prepare_household_surveys.py` for new (intermediate) household survey results.
* `09_household_survey_analysis.py` for new household survey results.

Note: if SSB data is also updated, then all scripts needs to be re-run.

Running the `01_prepare_exiobase.py` can take some time and use a decent amount of RAM memory (around 20 GB). See [Preparing the EXIOBASE raw data](#preparing-the-exiobase-raw-data) for more details.

<br>

### Preparing the EXIOBASE raw data

**This step is optional!**

We include the prepared GHG and monetary imports for Norway in the repository, so you can skip this step if you want to run the model with the prepared data. 

In case you want to rerun the pre-processing of the EXIOBASE data with the `01_prepare_exiobase.py` script: 
1. Manually download the data and place the data correctly in the folder `/data/01_raw/exiobase/`
2. Update the `config.py` file by setting `download_exiobase = True` and `exiobase_doi` to the new version of EXIOBASE.

In case the EXIOBASE data changes folder structure, one has two options:
1. Place the data manually in the correct folders and subfolders. See the `placeholder.txt` files in the `data/01_raw/exiobase/` subfolders for a guide. 
2. Update the `01_prepare_exiobase.py` script to properly read the data from the new folder structure.

<br>

### Linking the Norwegian IO data with EXIOBASE

This is done with the `05_SNAC.py` script. It requires the EXIOBASE data to be pre-processed (see above). 

<br>

## Development guide

The code is provided as percent formatted Python scripts (.py), which can be run interactively in most modern editors ([Visual Studio Code recommended](https://code.visualstudio.com/docs/python/jupyter-support-py)). 
