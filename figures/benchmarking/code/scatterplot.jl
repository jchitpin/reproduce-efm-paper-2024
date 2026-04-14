## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the aggregated carbon and nitrogen CHMC output files.
# (2) Aggregates the run times of the atomic CHMCs.
# (3) Plots run times as a function of number of atomic (carbon/nitrogen) efms.
# ------------------------------------------------------------------------------

## DIRECTORIES -----------------------------------------------------------------
# Dataset directories
datasets = ["e_coli_core", "iAB_RBC_283", "iIT341", "iSB619", "hepg2"]
dataset_names = ["E. coli core", "iAB RBC 283", "iIT341", "iSB619", "HepG2"]

# Import and export directories
load1 = "data/"
load2 = "/gems/"

# Import filenames (without suffix and .jld2)
im_x_carbon   = "gems-aggregated/gems-carbon.jld2"
im_x_nitrogen = "gems-aggregated/gems-nitrogen.jld2"
im_y_carbon   = "/efm-enumeration-times-carbon" 
im_y_nitrogen = "/efm-enumeration-times-nitrogen" 
im_suffix     = ["", "-T2", "-T4", "-T8", "-T16", "-T32"]

# Export filenames (without suffix and .tex)
ex_fig = "figures/benchmarking/subpanel-c//figure-scatterplot-julia-benchmark"
ex_suffix     = ["", "-T2", "-T4", "-T8", "-T16", "-T32"]
ex_fig_combined_log_log = "figures/benchmarking/subpanel-c/figure-scatterplot-julia-benchmark-combined-log-log.tex"
ex_fig_combined_log_lin = "figures/benchmarking/subpanel-c/figure-scatterplot-julia-benchmark-combined-log-lin.tex"
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
# Load atomic CHMCs
C = load_aggregated_chmcs(load1 * im_x_carbon, datasets)
N = load_aggregated_chmcs(load1 * im_x_nitrogen, datasets)

# Load atomic CHMC run times
Cy = [#
    load_atomic_chmc_times(#
        load1 * load2, datasets, im_y_carbon * s * ".jld2", "ctimes"
    )
    for s in im_suffix
]
Ny = [#
    load_atomic_chmc_times(#
        load1 * load2, datasets, im_y_nitrogen * s * ".jld2", "ntimes"
    )
    for s in im_suffix
]

# Get x/y coordinates of number of atomic EFMs and run time
xy_c = [get_efm_total_vs_time_benchmark(C, c) for c in Cy]
xy_n = [get_efm_total_vs_time_benchmark(N, n) for n in Ny]
# ------------------------------------------------------------------------------

## PLOT BENCHMARK TIMES --------------------------------------------------------
# xy_c[A][B][C] # A for threads; B = 1 for x-coord, 2 for y-coord; C'th dataset
# xy_n[A][B][C] # A for threads; B = 1 for x-coord, 2 for y-coord; C'th dataset

# Total number of carbon EFMs (for manual subpanel-a)
sum(xy_c[1][1][1]) # 6,519 E. coli carbon EFMs
sum(xy_c[1][1][2]) # 7,958 iAB_RBC_283 carbon EFMs
sum(xy_c[1][1][3]) # 30,633 iIT341 carbon EFMs
sum(xy_c[1][1][4]) # 59,444 iSB619 carbon EFMs
sum(xy_c[1][1][5]) # 2,180,613 HepG2 carbon EFMs

# Total number of nitrogen EFMs (for manual subpanel-a)
sum(xy_n[1][1][1]) # 204 E. coli carbon EFMs
sum(xy_n[1][1][2]) # 2,091 iAB_RBC_283 carbon EFMs
sum(xy_n[1][1][3]) # 10,440 iIT341 carbon EFMs
sum(xy_n[1][1][4]) # 14,214 iSB619 carbon EFMs
sum(xy_n[1][1][5]) # 5,105 HepG2 carbon EFMs

# Carbon Julia serial time
sum(xy_c[1][2][1]) # E. coli carbon time
sum(xy_c[1][2][2]) # iAB_RBC_283 carbon time
sum(xy_c[1][2][3]) # iIT341 carbon time
sum(xy_c[1][2][4]) # iSB619 carbon time
sum(xy_c[1][2][5]) # HepG2  carbon time

# Nitrogen Julia serial time
sum(xy_n[1][2][1]) # E. coli nitrogen time
sum(xy_n[1][2][2]) # iAB_RBC_283 nitrogen time
sum(xy_n[1][2][3]) # iIT341 nitrogen time
sum(xy_n[1][2][4]) # iSB619 nitrogen time
sum(xy_n[1][2][5]) # HepG2  nitrogen time

# Total Julia serial time
sum(xy_c[1][2][1]) + sum(xy_n[1][2][1]) # E. coli carbon+nitrogen time
sum(xy_c[1][2][2]) + sum(xy_n[1][2][2]) # iAB_RBC_283 carbon+nitrogen time
sum(xy_c[1][2][3]) + sum(xy_n[1][2][3]) # iIT341 carbon+nitrogen time
sum(xy_c[1][2][4]) + sum(xy_n[1][2][4]) # iSB619 carbon+nitrogen time
sum(xy_c[1][2][5]) + sum(xy_n[1][2][5]) # HepG2 carbon+nitrogen time

# JULIA LOG-LOG SCALING PLOTS
for i in eachindex(ex_suffix)
    p = @pgf Axis(#
        {#
            height = "7cm",
            width  = "7cm",
            "xmode = log",
            "ymode = log",
            "xlabel = AEFMs (\\#)",
            "ylabel = Serial run time (s)",
            "legend cell align = {left}",
            xmajorgrids = false,
            ymajorgrids = false,
            "xtick pos = bottom",
            "ytick pos = left",
            "legend style={font=\\tiny}",
            "legend entries={Carbon, Nitrogen, E. coli core, iAB RBC 283, iIT341, iSB619, HepG2}",
            "legend pos=north west",
            "legend image post style={scale=0.5}",
            xmax = 1000000,
            xtick = [1, 100, 10000, 1000000],
            ymax = 1000000,
            ytick = [0.0001, 1, 10000],
        },
        "\\addlegendimage{only marks,mark=o}",
        "\\addlegendimage{only marks,mark=square}",
        "\\addlegendimage{no markers, cred}",
        "\\addlegendimage{no markers, cblue}",
        "\\addlegendimage{no markers, cgreen}",
        "\\addlegendimage{no markers, cpurple}",
        "\\addlegendimage{no markers, corange}",
        "\\addlegendimage{no markers, cyellow}",
        Plot(#
            {#
                "only marks",
                "black",
                "mark=*",
                "mark options={fill=cred, fill opacity=1.0}"
            },
            Table([xy_c[i][1][1], xy_c[i][2][1]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=square*",
                "mark options={fill=cred, fill opacity=1.0}"
            },
            Table([xy_n[i][1][1], xy_n[i][2][1]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=*",
                "mark options={fill=cgreen, fill opacity=1.0}"
            },
            Table([xy_c[i][1][2], xy_c[i][2][2]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=square*",
                "mark options={fill=cgreen, fill opacity=1.0}"
            },
            Table([xy_n[i][1][2], xy_n[i][2][2]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=*",
                "mark options={fill=cpurple, fill opacity=1.0}"
            },
            Table([xy_c[i][1][3], xy_c[i][2][3]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=square*",
                "mark options={fill=cpurple, fill opacity=1.0}"
            },
            Table([xy_n[i][1][3], xy_n[i][2][3]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=*",
                "mark options={fill=corange, fill opacity=1.0}"
            },
            Table([xy_c[i][1][4], xy_c[i][2][4]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=square*",
                "mark options={fill=corange, fill opacity=1.0}"
            },
            Table([xy_n[i][1][4], xy_n[i][2][4]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=*",
                "mark options={fill=cblue, fill opacity=1.0}"
            },
            Table([xy_c[i][1][5], xy_c[i][2][5]]),
        ),
        Plot(#
            {#
                "only marks",
                "black",
                "mark=square*",
                "mark options={fill=cblue, fill opacity=1.0}"
            },
            Table([xy_n[i][1][5], xy_n[i][2][5]]),
        ),
    )
    pgfsave(ex_fig * ex_suffix[i] * ".tex", p)
end

# COMBINED JULIA LOG-LINEAR/LOG-LOG SCALING PLOTS
xcs = [reduce(vcat, xy_c[i][1]) for i in eachindex(xy_c)]
ycs = [reduce(vcat, xy_c[i][2]) for i in eachindex(xy_c)]
xns = [reduce(vcat, xy_n[i][1]) for i in eachindex(xy_n)]
yns = [reduce(vcat, xy_n[i][2]) for i in eachindex(xy_n)]

xs = [[xcs[i]; xns[i]] for i in eachindex(xcs)]
ys = [[ycs[i]; yns[i]] for i in eachindex(ycs)]

p = @pgf Axis(#
    {#
        height = "7cm",
        width  = "7cm",
        "xmode = log",
        "ymode = log",
        "xlabel = AEFMs (\\#)",
        "ylabel = Run time (s)",
        "legend cell align = {left}",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        "legend style={font=\\tiny}",
        "legend pos=north west",
        "legend image post style={scale=1.0}",
        xmax = 1000000,
        xtick = [1, 100, 10000, 1000000],
        ymax = 10000,
        ytick = [0.0001, 0.01, 1, 100, 10000],
    },
    raw"""\addlegendimage{empty legend}""",

    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cred, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[1], ys[1]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cgreen, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[2], ys[2]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cpurple, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[3], ys[3]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=corange, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[4], ys[4]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cblue, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[5], ys[5]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cyellow, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[6], ys[6]),
    ),
    "\\addlegendentry{\\hspace{-.23cm}\\textbf{Threads (\\#)}}",
    "\\addlegendentry{1}",
    "\\addlegendentry{2}",
    "\\addlegendentry{4}",
    "\\addlegendentry{8}",
    "\\addlegendentry{16}",
    "\\addlegendentry{32}",
)
pgfsave(ex_fig_combined_log_log, p)

p = @pgf Axis(#
    {#
        height = "7cm",
        width  = "7cm",
        "xmode = log",
        "xlabel = AEFMs (\\#)",
        "ylabel = Run time (s)",
        "legend cell align = {left}",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        "legend style={font=\\tiny}",
        "legend pos=north west",
        "legend image post style={scale=1.0}",
        xmax = 1000000,
        xtick = [1, 100, 10000, 1000000],
        ymin = -250,
        ymax = 2500,
        ytick = [0, 500, 1000, 1500, 2000, 2500]
    },
    raw"""\addlegendimage{empty legend}""",

    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cred, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[1], ys[1]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cgreen, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[2], ys[2]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cpurple, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[3], ys[3]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=corange, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[4], ys[4]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cblue, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[5], ys[5]),
    ),
    Plot(#
        {#
            "only marks",
            "black",
            "mark=*",
            "mark options={fill=cyellow, fill opacity=1.0, draw opacity=1.0}"
        },
        Table(xs[6], ys[6]),
    ),
    "\\addlegendentry{\\hspace{-.23cm}\\textbf{Threads (\\#)}}",
    "\\addlegendentry{1}",
    "\\addlegendentry{2}",
    "\\addlegendentry{4}",
    "\\addlegendentry{8}",
    "\\addlegendentry{16}",
    "\\addlegendentry{32}",
)
pgfsave(ex_fig_combined_log_lin, p)

