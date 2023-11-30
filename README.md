# CaFEAN: Code for the coupled footprint/simplified-SNAC model

This repository contains the Python code for the CaFEAN model, linking the Multi-Regional Input-Output (MRIO) model EXIOBASE with the official Norwegian (SSB) IO and greenhouse gas emissions data in order to produce carbon footprint results for Norway.

The readme is structured as follows:

- [Repository Structure](#repository-structure)
    - [Scripts and environment](#scripts-and-environment)
    - [Data folder](#data-folder)
    - [Result format](#result-format)
- [Running the code](#running-the-code)
    - [Setting up the environment](#setting-up-the-environment)
    - [Visualizing the results](#visualizing-the-results)
    - [Use cases](#use-cases)
- [Development guide](#development-guide)

The CaFEAN model is licensed under the AGPLv3 license. See the [LICENSE](./LICENCE.md) file for details.

Copyright (C) <2023>  <XIO Sustainability Analytics, Inc>

## Repository Structure

### Scripts and environment

There are 5 main python files in the root folder, except for the final visualization all in a python script (.py) and a Jupyter notebook (.ipynb) format:

- ./exio_prepper.py|.ipynb : This file contains the code for the pre-processing of the EXIOBASE data which extracts Norwegian monetary and GHG imports from EXIOBASE. 

- ./nor_exio_snac.py|.ipynb : This file contains the code for the simplified-SNAC (or coupled) model, linking the Norwegian IO data with the EXIOBASE data.

- ./aggregator.py|.ipynb : This file contains the code for aggregating the results of the simplified-SNAC model

- ./excel_converter.py|.ipynb : This file contains the code for converting the results of the simplified-SNAC model/aggregation to pre-prepared table formats to copy to Excel for some basic visualisations.

- ./vis_app.py : An interactive visualization app based on streamlit. 

In addition, there are

- ./environment.yml : This file contains the environment for the development of the code.
- ./jupytext.toml : Defines the pairing of the .py and .ipynb files [see jupytext readme for details](https://jupytext.readthedocs.io/en/latest/config.html).

### Data folder

The folder ./data contains all the data required to run the code. It is structured as follows:

- ./data/emis_conv_char.xlsx : This file contains the conversion factors for the GHG emissions:

    - from EXIOBASE to Norwegian emission data (unit and name conversion)
    - characterisation factors for emission data in the Norwegian format

    All unit and name definitions are defined in that file.

- ./data/exchange_rates.csv : The exchange rates for the conversion from EXIOBASE monetary units to Norwegian monetary units. This data was obtained from Eurostat. 

- ./data/iot_XXXX_XXXX.xlsx : The official Norwegian IO data. 

- ./data/F_xxx_CO2biogenic.csv : The biogenic CO2 emissions for the EXIBOASE 3.8.2 dataset. This gets linked to EXIOBASE in the pre-processing step within the *exio_prepper* script.

- ./data/AEA Questionnaire_2023_Norway2023_til footprint.xlsm : The official Norwegian GHG emissions data.

- ./data/sector_list.xlsx Names and codes of sectors in EXIOBASE and the official Norwegian IO data.

- ./data/sector_matching.xlsx : This file contains the matching of the sectors in the official Norwegian IO data with the sectors in EXIOBASE.

- ./data/exio_no_bp_raw.xlsx : Data extracted from EXIOBASE for the snac-coupling. This file is created by the *exio_prepper* script. With:

    - imp_xxxx: Monetary and emission imports based on EXIOBASE data for year xxxx
    - imp_source_xxxx : Source of import footprints based on EXIOBASE data for year xxxx (units as in imp_xxxx, totals add up to the same value)
    - fd_xxxx: Final demand and total footprints based on EXIBASE data for year xxxx. This is just for EXIOBASE reporting and not used futher.
    - fp_source_xxxx : Source of footprints based on EXIOBASE data for year xxxx (units as in fd_xxxx, totals add up to the same value). This is just for EXIOBASE reporting and not used futher.

- ./data/region_agg_spec.tsv : Aggregation specification for the regions. This file is used by the *aggregator* script.

- ./data/sector_agg_spec.tsv : Aggregation specification for the sectors. This file is used by the *aggregator* script.

### Result format

Result files are stored in a long table format. This can easily be converted to
any data format, including SQL. In excel, the list can easily be pivoted to a
wide format or filtered for specific years. In Python/Pandas, this can be
extended to a wide format with .unstack(index_name). The table also contains
the unit for each data entry. For further reprocessing in Python this might be
cumbersome and can be dropped with .drop('unit', axis=1). 

- ./results/sector_accounts.tsv : Sector accounts with the following columns:

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

- ./results/footprints_final_demand_breakdown.tsv : Final demand breakdown of the footprint accounts with the following columns:

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

- ./results/household_emissions.tsv : Household emissions (official Norwegian data) with the following columns:

    - year
    - emission
    - hhld_component (Household totals, Transport, Heating/cooling, Other, as per official Norwegian statistics)
    - value
    - unit

- ./results/total_accounts.tsv : These are total sums of the sector_accounts, also including household emissions. The columns are:

    - year
    - component: 
        - totals of all components of sector_accounts.tsv 
        - household_totals: total direct household emissions (same as Household totals in household_emissions.tsv / official statistics)
        - total_footprint_with_households: total footprint including household emissions 
        - production_account_with_households: total production based account including household emissions 
    - emission
    - value
    - unit

- ./results/household_emissions.tsv : Household emissions (official Norwegian data) with the following columns:

    - year
    - emission
    - hhld_component (Household totals, Transport, Heating/cooling, Other, as per official Norwegian statistics)
    - value
    - unit

- ./results/footprint_sources.tsv : Source of imports

    - year
    - emission
    - region (region where the emissions occur/imported from)
    - sector
    - value
    - unit

- ./results/footprint_sources_totals.tsv : Total source of imports (all sectors aggregated)

    - year
    - emission
    - region (region where the emissions occur/imported from)
    - value
    - unit

- ./results/excel_report_format : Wide table format for excel import of the data underlying the figures of the report. 

Besides these data files, the *aggregator* script also creates aggregated versions
of these files, indicated by the suffix "_agg". These files have the same 
structure as explained above, with sectors/regions aggregated.

## Running the code

The code is setup to run with [version 3.8.2 of EXIOBASE](https://doi.org/10.5281/zenodo.5589597).

A full run of the code consists of the following steps:

- [building the environment (first time)](#setting-up-the-environment) - only necessary the first time 
- activating the environment:

        conda activate cafean

- [Downloading and extracting the EXIOBASE raw data](#preparing-the-exiobase-raw-data). This step is the most time/resource intensive. It is not necessary to run it when using EXIOBASE 3.8.2 as the extracted output is included in the code.

        python exio_prepper.py

- [Linking the Norwegian IO data with EXIOBASE results](#linking-the-norwegian-io-data-with-exiobase)

        python nor_exio_snac.py

- Aggregating the results

        python aggregator.py

- Converting the results to an excel format (optional)

        python excel_converter.py

- [Visualizing the results](#visualizing-the-results)

        streamlit run vis_app.py

It is **not necessary to run the full sequence.** For example, if you only want to change the aggregation of the results, you can skip the first two steps and start with the *aggregator* script. See [section use cases below](#use-cases) for some standard work flows you might encounter.

In particular, the [visualization app](#visualizing-the-results) can be run without running any of the other scripts.

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


### Visualizing the results

The results can be visualized with the *vis_app* script. This is an interactive visualization app based on [streamlit](https://streamlit.io).
To run the app, simply activate the conda environment and run the script:

```bash
streamlit run vis_app.py
```

The app will guide you through the different visualizations.

### Use cases

Here we give some example use cases and how to run the code for them. We go from simple cases to more complex ones.

**WARNING:** Running the code overwrite previous results! Make sure to backup the results before running the code again. 
Alternatively, all scripts define the output folder in the first lines of the script.

#### Changing the aggregation of the results

For the report we aggregated the results to 17 sectors and 4 regions.
The aggregations are defined in ./data/region_agg_spec.tsv and ./data/sector_agg_spec.tsv.

To change the aggregation, simply change the aggregation specification in these files and rerun:

- the *aggregator* script to aggregate the results
- the *excel_converter* script to convert the results to an excel format (optional)

Then the new results can be visualized with the *vis_app* script.

It is NOT necessary to rerun the *exio_prepper* or the *nor_exio_snac* script.
The aggregation works on the results of the *nor_exio_snac* script.

#### Changing the exchange rates

Exchange rates are given in the file ./data/exchange_rates.csv. 
These are average Eurostat exchange rates. In case other data sources for these should be used, just change them in this file.
Then rerun: 

- the *nor_exio_snac* script to rerun the model with the new exchange rates
- the *aggregator* script to aggregate the results
- the *excel_converter* script to convert the results to an excel format (optional)

Then the new results can be visualized with the *vis_app* script.

It is NOT necessary to rerun the *exio_prepper* script, as the exchange rates are not used in the pre-processing of the EXIOBASE data.

#### Updating the Norwegian Emission or IO data

As long as the updated data keeps the same format/excel structure this can be done by simply replacing the old data with the new one in the ./data folder. 
Then rerun: 

- the *nor_exio_snac* script to link the new data with the EXIOBASE data
- the *aggregator* script to aggregate the results
- the *excel_converter* script to convert the results to an excel format (optional)

Then the new results can be visualized with the *vis_app* script.

Again, it is not necessary to rerun the *exio_prepper* script, as the Norwegian data is not used in the pre-processing of the EXIOBASE data.

In case the structure of the data changed, the *nor_exio_snac* script needs to be adapted accordingly.

#### Updating the GHG characterisation factors 

The characterisation factor from GHG emissions to GHG totals are based on the GWP100 factors of the IPCC Fifth Assessment Report.
These are given in the file ./data/emis_conv_char.xlsx. The files contains a sheet "meta" given some information on the structure.
In short, the characterisation factors are given in the sheet "nor_char" and can be updated there.

Then rerun: 

- the *exio_prepper* script to update the characterisation for the EXIOBASE data
- the *nor_exio_snac* script to characterise the Norwegian GHG emissions with the new factors and link to the new EXIOBASE data
- the *aggregator* script to aggregate the results
- the *excel_converter* script to convert the results to an excel format (optional)

Running the exio_prepper script takes some time and a decent amount of RAM memory (around 20 GB). See [Preparing the EXIOBASE raw data](#preparing-the-exiobase-raw-data) for more details.

Finally, the new results can be visualized with the *vis_app* script.

#### Running for a new EXIOBASE version

The code is setup to run with [version 3.8.2 of EXIOBASE](https://doi.org/10.5281/zenodo.5589597).
When new EXIOBASE data gets available, the script can be adapted to run with the new version.
For that, change the doi number in the *exio_prepper* script to the new version and rerun the script.

- the *exio_prepper* script to download the new EXIOBASE data and extract the Norwegian monetary and GHG imports
- the *nor_exio_snac* for linking to the new data
- the *aggregator* script to aggregate the results
- the *excel_converter* script to convert the results to an excel format (optional)

Running the exio_prepper script takes some time and a decent amount of RAM memory (around 20 GB). See [Preparing the EXIOBASE raw data](#preparing-the-exiobase-raw-data) for more details.

Finally, the new results can be visualized with the *vis_app* script.


### Preparing the EXIOBASE raw data

**This step is optional!**

We include the prepared GHG and monetary imports for Norway in the repository, so you can skip this step if you want to run the model with the prepared data. 

In case you want to rerun the pre-processing of the EXIOBASE data, the *exio_prepper* script downloads the EXIOBASE raw data from the Zenodo repository and extracts the Norwegian monetary and GHG imports. This can take a while and also requires a lot of RAM memory (around 20 GB). 

### Linking the Norwegian IO data with EXIOBASE

This is done with the *nor_exio_snac* script. It requires the EXIOBASE data to be pre-processed (see above). 


## Development guide

The code is provided as Jupyter Notebooks (.ipynb) and Python scripts (.py). 
These two are equivalent and can be used interchangeably.
They are linked via [jupytext](https://jupytext.readthedocs.io/en/latest/), thus any changes in one file will be reflected in the other.

The .py version is formatted in [percent
format](https://jupytext.readthedocs.io/en/latest/formats-scripts.html) and the
linking is defined in ./jupytext.toml (so no further manual pairing necessary).


