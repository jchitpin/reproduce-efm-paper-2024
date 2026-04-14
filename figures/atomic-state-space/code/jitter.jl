## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the aggregated carbon and nitrogen CHMC output files and
#     precomputed dictionaries.
# (2) Compute the total number of MC states for each metabolite/atom.
#     
# (3) Export the statistics to dataframe for plotting.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
datasets = ["e_coli_core", "iAB_RBC_283", "iIT341", "iSB619", "hepg2"]

# Import and export directories
load1 = "data/"
load2 = "gems/"
load3 = "/processed/"
save1 = "figures/atomic-state-space/subpanels/"

# Import filenames
im_carbon_chmcs    = "/gems-aggregated/gems-carbon.jld2"
im_nitrogen_chmcs  = "/gems-aggregated/gems-nitrogen.jld2"

im_carbon_dicts = "dictionary-atom-tracing-carbon.csv"
im_nitrogen_dicts = "dictionary-atom-tracing-nitrogen.csv"

# Export filenames
ex_carbon_mc_states   = "jitter-carbon.dat"
ex_nitrogen_mc_states = "jitter-nitrogen.dat"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using Statistics
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))

# Setting seed for datapoint jittering!
Random.seed!(0)
# ------------------------------------------------------------------------------

## LOAD AGGREGATED DATASETS ----------------------------------------------------
@info("Loading aggregated datasets.")
C = load_aggregated_chmcs(load1 * im_carbon_chmcs, datasets)
N = load_aggregated_chmcs(load1 * im_nitrogen_chmcs, datasets)

# Load atom tracing dictionary
D_C = [#
    import_atom_tracing_dictionary(load1 * load2 * d * load3 * im_carbon_dicts)
    for d in datasets
]
D_N = [#
    import_atom_tracing_dictionary(load1 * load2 * d * load3 * im_nitrogen_dicts)
    for d in datasets
]
# ------------------------------------------------------------------------------

## COMPUTE SUMMARY STATISTICS FOR EACH GEM -------------------------------------
# Carbon and nitrogen atomic state counts
C1 = get_atomic_mc_state_space(C)
N1 = get_atomic_mc_state_space(N)

C_num = [reduce(vcat, c) for c in C1]
N_num = [reduce(vcat, n) for n in N1]

# Total carbon and nitrogen atomic states
C_denom = length.(D_C) # 1127, 3444, 6273, 9896, 8194 total carbon atom states
N_denom = length.(D_N) # 265, 944, 1677, 2533, 1908 total nitrogen atom states

C_frac = C_num ./ C_denom
N_frac = N_num ./ N_denom

mean.(C_frac) # 0.013, 0.0025, 0.0035, 0.0021, 0.018 percent
mean.(N_frac) # 0.018, 0.012, 0.020, 0.014, 0.034 percent

mean(reduce(vcat, C_frac)) # 0.0041
mean(reduce(vcat, N_frac)) # 0.017

# ------------------------------------------------------------------------------

## RESHAPE DATA INTO DATAFRAMES ------------------------------------------------
# Set jitter (variance of a normal distribution)
jit = 0.1

# Assign labels based on source metabolite super class defined by HMDB
# x = dataset index, y = log10 statistic, z = jittered dataset index
dfc = get_atomic_state_count_summary_dataframe(C_frac, C_num, jit)
dfn = get_atomic_state_count_summary_dataframe(N_frac, N_num, jit)
# ------------------------------------------------------------------------------

## EXPORT SUPER CLASS MATRICES TO DATA FILES -----------------------------------
CSV.write(save1 * ex_carbon_mc_states,   dfc, header = true, delim = " ")
CSV.write(save1 * ex_nitrogen_mc_states, dfn, header = true, delim = " ")
# ------------------------------------------------------------------------------

