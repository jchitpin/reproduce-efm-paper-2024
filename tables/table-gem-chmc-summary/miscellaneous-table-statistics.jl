## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the aggregated carbon and nitrogen CHMC statistics.
# (2) Computes some additional statistics reported in text of the manuscript.
# ------------------------------------------------------------------------------

# Export filenames
im_df    = "tables/table-gem-chmc-summary/table-gem-chmc-summary.csv"

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using MarkovWeightedEFMs, CSV, Tables, JLD2, DataFrames, Statistics
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

df = CSV.read(im_df, DataFrame)

ratioC1 = df[:,:CHMC_Transitions_Carbon] ./ df[:,:CHMC_State_Space_Carbon]
ratioN1 = df[:,:CHMC_Transitions_Nitrogen] ./ df[:,:CHMC_State_Space_Nitrogen]

mean(ratioC1) # 1.43
std(ratioC1)  # 0.14

mean(ratioN1) # 1.53
std(ratioN1)  # 0.16
# ------------------------------------------------------------------------------

