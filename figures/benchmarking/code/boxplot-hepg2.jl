## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Computes the boxplot statistics for the top 33 HepG2 carbon CHMC run
# times across serial and multithreads.
# (2) Computes the slope and intercept of the top 33 HepG2 carbon CHMC serial
# run times.
# ------------------------------------------------------------------------------

## DIRECTORIES -----------------------------------------------------------------
# Dataset directories
datasets = ["e_coli_core", "iAB_RBC_283", "iIT341", "iSB619", "hepg2"]
dataset_names = ["E. coli core", "iAB RBC 283", "iIT341", "iSB619", "HepG2"]

# Import and export directories
load1 = "data/"
load2 = "/gems/"

# Import filenames (without suffix and .jld2)
im_carbon_chmcs    = "/gems-aggregated/gems-carbon.jld2"
im_y_carbon        = "/efm-enumeration-times-carbon" 
im_y_nitrogen      = "/efm-enumeration-times-nitrogen" 
im_suffix          = ["", "-T2", "-T4", "-T8", "-T16", "-T32"]

# Export filenames (without suffix and .tex)
ex_fig = "figures/benchmarking/subpanel-d/figure-boxoplot-hepg2-julia.tex"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## CUSTOM TIKZ PREAMBLE --------------------------------------------------------
# Define custom colours
pgfplotsx_preamble()
# ------------------------------------------------------------------------------

## LOAD DATASETS ---------------------------------------------------------------
# Load atomic CHMC run times
C = [#
    load_atomic_chmc_times(#
        load1 * load2, datasets, im_y_carbon * s * ".jld2", "ctimes"
    )
    for s in im_suffix
]
N = [#
    load_atomic_chmc_times(#
        load1 * load2, datasets, im_y_nitrogen * s * ".jld2", "ntimes"
    )
    for s in im_suffix
]
# ------------------------------------------------------------------------------
#
## COMPUTE SPEEDUP AND BOXPLOT STATISTICS --------------------------------------
# Find all source metabolite/atoms with serial times ≥ 5 for HepG2 GEM
min_threshold = 5
src_met_atoms = Vector{Tuple{Int64, Int64}}()
for j in 1:length(C[1][5])
    if !isempty(C[1][5][j]) && any(C[1][5][j] .>= min_threshold)
        for l in findall(C[1][5][j] .>= min_threshold)
            push!(src_met_atoms, (j, l))
        end
    end
end

hepg2_serial_thread_times = Vector{Vector{Float64}}(undef, length(C))
for i in eachindex(C)
    tmp_i = Vector{Float64}()
    for j in src_met_atoms
        push!(tmp_i, C[i][5][j[1]][j[2]])
    end
    hepg2_serial_thread_times[i] = tmp_i
end

# Compute speedup
hepg2_speedup = [hepg2_serial_thread_times[1] ./ h for h in hepg2_serial_thread_times]

# Max speedup across serial/parallel runs and the 33 CHMCS:
maximum(reduce(vcat, hepg2_speedup))

# Statistics for the biggest ACHMC (HepG2 dataset)
Carbon_ACHMCs = load_aggregated_chmcs(load1 * im_carbon_chmcs, datasets)
idx = src_met_atoms[findfirst(==(maximum(hepg2_speedup[6])), hepg2_speedup[6])]
length(Carbon_ACHMCs[5][idx[1]][idx[2]].dchmc) # 1,120,290 CHMC states
length(Carbon_ACHMCs[5][idx[1]][idx[2]].dmc) # 1,158 MC states
length(Carbon_ACHMCs[5][idx[1]][idx[2]].e) # 168,793 EFMs
hepg2_speedup[6][findfirst(==(maximum(hepg2_speedup[6])), hepg2_speedup[6])] # 9.26 times speedup

# Boxplot statistics
hepg2_median = [median(h) for h in hepg2_speedup]
hepg2_lq = [percentile(h, 25) for h in hepg2_speedup]
hepg2_uq = [percentile(h, 75) for h in hepg2_speedup]
hepg2_whiskers = [lower_upper_whisker(h) for h in hepg2_speedup[2:end]]
hepg2_lw = first.(hepg2_whiskers)
hepg2_uw = last.(hepg2_whiskers)
pushfirst!(hepg2_lw, 1.0)
pushfirst!(hepg2_uw, 1.0)
# ------------------------------------------------------------------------------

## Calculating slope of the largest set of HepG2 scatterplot points ------------
xy_c = [get_efm_total_vs_time_benchmark(Carbon_ACHMCs, c) for c in C]
y = hepg2_serial_thread_times[1] # times

xy_c[1][1][5] # number of carbon EFMs for HepG2
xy_c[1][2][5] # serial times for HepG2

ids = [findfirst(xy_c[1][2][5] .== y[i]) for i in eachindex(y)]
x = xy_c[1][1][5][ids]

using LinearRegression

lr = linregress(log10.(x), log10.(y))
lr.coeffs[1] # slope:     +1.90
lr.coeffs[2] # intercept: -6.54
# log₁₀(run time (s)) = log₁₀(# AEFMs)*1.9 - 6.54
# ------------------------------------------------------------------------------


