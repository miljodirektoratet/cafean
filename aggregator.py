# -*- coding: utf-8 -*-
# %% [markdown]
#  # Small aggregator utility script
#
# This can be used to aggregate sector and/or regions of the results produced by the `nor_exio_snac` script.
#
# The aggregation keys are read from
#
# ./data/sector_agg_spec.tsv
# ./data/region_agg_spec.tsv
#
# Changing the aggregation keys should happen in these files, then rerun the script.
#
# Requires pandas
#

# %% [markdown]
# Import of required external packages.

# %%
from pathlib import Path
import pandas as pd

# %% [markdown]
# Folder where the results are stored.

# %%
# Set the work path to the directory where this script is located
# and if this is not available to the current working directory
try:
    work_path = Path(__file__).parent.absolute()  # when running as script
except NameError:
    work_path = Path.cwd()

results_folder = work_path / "results"

# %% [markdown]
# File(s) to aggregate.
# The script can also be used to just aggregate a specific file,
# just pass the path to that file within a list.

# %%
src_files = list(results_folder.glob("*.tsv"))

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
    pd.read_csv(work_path / "data" / "region_agg_spec.tsv", sep="\t")
    .apply(lambda x: x.str.strip())
    .set_index("src")
    .to_dict(orient="dict")
    .get("agg")
)

# sec_agg = None
sec_agg = pd.read_csv(work_path / "data" / "sector_agg_spec.tsv", sep="\t")

sec_agg = (
    pd.read_csv(work_path / "data" / "sector_agg_spec.tsv", sep="\t")
    .apply(lambda x: x.str.strip())
    .set_index("src")
    .to_dict(orient="dict")
    .get("agg")
)

# %% [markdown]
# The actual aggregation function


# %%
def aggregate(df, agg_dict, level_to_agg):
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
