## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the aggregated carbon and nitrogen CHMC output files.
# (2) Calculates the length of each atomic EFM in each GEM.
# (3) Plots the histogram curves and exports the LaTeX (.tex) figure.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
datasets = ["e_coli_core", "iAB_RBC_283", "iIT341", "iSB619", "hepg2"]
dataset_names = ["E. coli core", "iAB RBC 283", "iIT341", "iSB619", "HepG2"]

# Import and export directories
load1 = "data/"
load2 = "/gems/"

# Import filenames
im_carbon   = "/gems-aggregated/gems-carbon.jld2"
im_nitrogen = "/gems-aggregated/gems-nitrogen.jld2"

# Export filenames
ex_fig_c = "figures/structural/subpanel-c-d/figure-efm-length-curves-carbon.tex"
ex_fig_n = "figures/structural/subpanel-c-d/figure-efm-length-curves-nitrogen.tex"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
#using MarkovWeightedEFMs, JLD2, PGFPlotsX, StatsBase, LinearAlgebra
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD AGGREGATED DATASETS ----------------------------------------------------
@info("Loading aggregated datasets.")
C = load_aggregated_chmcs(load1 * im_carbon, datasets)
N = load_aggregated_chmcs(load1 * im_nitrogen, datasets)
# ------------------------------------------------------------------------------

## COLLECT EFM LENGTHS ---------------------------------------------------------
efm_c = get_efm_length(C)
efm_n = get_efm_length(N)
# ------------------------------------------------------------------------------

## CUSTOM TIKZ PREAMBLE --------------------------------------------------------
# Define custom colours
pgfplotsx_preamble()
# ------------------------------------------------------------------------------

## PLOT EFM LENGTH CURVES ------------------------------------------------------
fit_mle(Normal, efm_c[5]) # μ = 61.71; σ² = 15.46
fit_mle(Normal, efm_n[5]) # μ = 12.03; σ² = 8.17

histc = Vector{Vector{Vector{Real}}}(undef, length(datasets))
maximum(maximum.(reduce(vcat, efm_c)))
xmax = 120
nbinsc = [20, 20, 20, 20, 20]
for i in 1:length(datasets)
    hst = fit(Histogram, efm_c[i], closed = :left, nbins = nbinsc[i])
    hst = LinearAlgebra.normalize(hst, mode = :pdf)
    x = (collect(hst.edges[1])[2:end] .+ collect(hst.edges[1])[1:end-1]) / 2
    histc[i] = [[0; x; xmax], [0; hst.weights; 0]]
end

pc = @pgf Axis(#
    {#
        height = "7cm",
        width  = "7cm",
        no_markers,
        "xlabel = Carbon AEFM length",
        "ylabel = PDF",
        "scaled y ticks = base 10:2",
        "legend cell align = {left}",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        xmax = xmax,
        xtick = [0, 30, 60, 90, 120],
        ymax = 0.15,
        ytick = [0, 0.05, 0.1, 0.15]
    },
    Plot(#
        {"cred"},
        Table(histc[1][1], histc[1][2])
    ),
    LegendEntry(dataset_names[1]),
    Plot(#
        {"cgreen"},
        Table(histc[2][1], histc[2][2])
    ),
    LegendEntry(dataset_names[2]),
    Plot(#
        {"cpurple"},
        Table(histc[3][1], histc[3][2])
    ),
    LegendEntry(dataset_names[3]),
    Plot(#
        {"corange"},
        Table(histc[4][1], histc[4][2])
    ),
    LegendEntry(dataset_names[4]),
    Plot(#
        {"cblue"},
        Table(histc[5][1], histc[5][2])
    ),
    LegendEntry(dataset_names[5]),
)

# Nitrogen
histn = Vector{Vector{Vector{Real}}}(undef, length(datasets))
maximum(maximum.(reduce(vcat, efm_n)))
xmax = 40
nbinsn = [5, 15, 10, 10, 10]
for i in 1:length(datasets)
    hst = fit(Histogram, efm_n[i], closed = :left, nbins = nbinsn[i])
    hst = LinearAlgebra.normalize(hst, mode = :pdf)
    x = (collect(hst.edges[1])[2:end] .+ collect(hst.edges[1])[1:end-1]) / 2
    histn[i] = [[0; x; xmax], [0; hst.weights; 0]]
end

pn = @pgf Axis(#
    {#
        height = "7cm",
        width  = "7cm",
        no_markers,
        "xlabel = Nitrogen AEFM length",
        "ylabel = PDF",
        "scaled y ticks = base 10:1",
        "legend cell align = {left}",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        xmax = xmax,
        xtick = [0, 10, 20, 30, 40],
        ymax = 0.6,
        "extra y ticks={-0.06}",
        "extra y tick labels={10}",
        "extra y tick style={yticklabel style={color=white}}",
        ytick = [0, 0.2, 0.4, 0.6],
        yticklabels = [0, 2, 4, 6],
    },
    Plot(#
        {"cred"},
        Table(histn[1][1], histn[1][2])
    ),
    LegendEntry(dataset_names[1]),
    Plot(#
        {"cgreen"},
        Table(histn[2][1], histn[2][2])
    ),
    LegendEntry(dataset_names[2]),
    Plot(#
        {"cpurple"},
        Table(histn[3][1], histn[3][2])
    ),
    LegendEntry(dataset_names[3]),
    Plot(#
        {"corange"},
        Table(histn[4][1], histn[4][2])
    ),
    LegendEntry(dataset_names[4]),
    Plot(#
        {"cblue"},
        Table(histn[5][1], histn[5][2])
    ),
    LegendEntry(dataset_names[5]),
)
# ------------------------------------------------------------------------------

## SAVE EFM LENGTH CURVES ------------------------------------------------------
pgfsave(ex_fig_c, pc) # carbon
pgfsave(ex_fig_n, pn) # nitrogen
# ------------------------------------------------------------------------------



