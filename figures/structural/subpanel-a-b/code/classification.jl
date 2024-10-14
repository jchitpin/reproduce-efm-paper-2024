## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the aggregated carbon and nitrogen CHMC output files.
# (2) Calculates the fraction of pathed-versus-looped carbon/nitrogen EFMs.
# (3) Computes fraction of weight in pathed-versus-looped carbon/nitrogen EFMs.
# (3) Plots the barplot and exports the LaTeX (.tex) figure.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
datasets = ["e_coli_core", "iAB_RBC_283", "iIT341", "iSB619", "hepg2"]
dataset_names = ["E. coli core", "iAB RBC 283", "iIT341", "iSB619", "HepG2"]

# Import and export directories
load1 = "data/"
load2 = "/gems/"
save1 = "figures/structural/subpanel-a-b/"

# Import filenames
im_carbon         = "/gems-aggregated/gems-carbon.jld2"
im_nitrogen       = "/gems-aggregated/gems-nitrogen.jld2"
im_carbon_hepg2   = "gems/hepg2/atomic-efm-weights-carbon"
im_nitrogen_hepg2 = "gems/hepg2/atomic-efm-weights-nitrogen"

# Export filenames
ex_fig_c1    = save1 * "/figure-barplot-tree-classification-carbon.tex"
ex_fig_n1    = save1 * "/figure-barplot-tree-classification-nitrogen.tex"
ex_fig_c2    = save1 * "/figure-barplot-tree-classification-carbon-scaled.tex"
ex_fig_n2    = save1 * "/figure-barplot-tree-classification-nitrogen-scaled.tex"
ex_fig_hepg2 = save1 * "/figure-barplot-pathed-versus-looped-efm-flux-hepg2.tex"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD AGGREGATED DATASETS ----------------------------------------------------
@info("Loading aggregated datasets.")
C = load_aggregated_chmcs(load1 * im_carbon, datasets)
N = load_aggregated_chmcs(load1 * im_nitrogen, datasets)

Ch = load_hepg2_chmcs(load1 * im_carbon_hepg2)
Nh = load_hepg2_chmcs(load1 * im_nitrogen_hepg2)
# ------------------------------------------------------------------------------

## CLASSIFY ATOMIC EFMS --------------------------------------------------------
class_C = classify_atomic_efms(C)
class_N = classify_atomic_efms(N)

# Scale classifications into percentages
format_C = [[t.n_pathed_no_rec, t.n_pathed_rec, t.n_looped_no_rec, t.n_looped_rec] for t in class_C]
format_N = [[t.n_pathed_no_rec, t.n_pathed_rec, t.n_looped_no_rec, t.n_looped_rec] for t in class_N]
percent_C = [100 .* (t ./ sum(t)) for t in format_C]
percent_N = [100 .* (t ./ sum(t)) for t in format_N]
# ------------------------------------------------------------------------------

## HEPG2 FLUX THROUGH S-T AND LOOPED ATOMIC EFMS -------------------------------
pathed_c, looped_c = sum_pathed_and_looped_efm_fluxes(Ch)
pathed_n, looped_n = sum_pathed_and_looped_efm_fluxes(Nh)
# ------------------------------------------------------------------------------

## CUSTOM TIKZ PREAMBLE --------------------------------------------------------
# Define custom colours
pgfplotsx_preamble()
# ------------------------------------------------------------------------------

## BARPLOTS --------------------------------------------------------------------
# Carbon
coords_c1_1 = [(i, format_C[i][1]) for i in eachindex(format_C)]
coords_c1_2 = [(i, format_C[i][2]) for i in eachindex(format_C)]
coords_c1_3 = [(i, format_C[i][3]) for i in eachindex(format_C)]
coords_c1_4 = [(i, format_C[i][4]) for i in eachindex(format_C)]

pc1 = @pgf SemiLogYAxis(
    {
        height = "7cm",
        width  = "7cm",
        "ybar = 0pt",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        "bar width = 6pt",
        "enlarge x limits = 0.15",
        "enlarge y limits = false",
        "legend image code/.code={\\draw [#1] (0cm,-0.1cm) rectangle (0.2cm,0.25cm); }",
        legend_style = {
            legend_columns = 2,
            "font=\\tiny",
        },
        xmax = length(datasets),
        xtick = [1, 2, 3, 4, 5],
        xticklabels = dataset_names,
        "xticklabel style = {font=\\tiny, align=center,rotate=45}",
        ylabel = "Carbon AEFM classes (\\#)",
        ymin = 0,
    },
        Plot(#
        {#
            "black",
            "fill=julia1"
        },
        Coordinates(coords_c1_1)
    ),
    Plot(#
        {#
            "black",
            "fill=julia2",
        },
        Coordinates(coords_c1_2)
    ),
    Plot(#
        {#
            "black",
            "fill=julia3"
        },
        Coordinates(coords_c1_3)
    ),
    Plot(#
        {#
            "black",
            "fill=julia4",
        },
        Coordinates(coords_c1_4)
    ),
    #Legend([#
        #"Src-to-sink (no revisit)", "Src-to-sink (revisit)",
        #"Loop (no revisit)", "Looped (revisit)",
    #]),
)

# Nitrogen
coords_n1_1 = [(i, format_N[i][1]) for i in eachindex(format_N)]
coords_n1_2 = [(i, format_N[i][2]) for i in eachindex(format_N)]
coords_n1_3 = [(i, format_N[i][3]) for i in eachindex(format_N)]
coords_n1_4 = [(i, format_N[i][4]) for i in eachindex(format_N)]

pn1 = @pgf SemiLogYAxis(
    {
        height = "7cm",
        width  = "7cm",
        "ybar = 0pt",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        "bar width = 6pt",
        "enlarge x limits = 0.15",
        "enlarge y limits = false",
        "legend image code/.code={\\draw [#1] (0cm,-0.1cm) rectangle (0.2cm,0.25cm); }",
        legend_style = {
            legend_columns = 2,
            "font=\\tiny",
        },
        xmax = length(datasets),
        xtick = [1, 2, 3, 4, 5],
        xticklabels = dataset_names,
        "xticklabel style = {font=\\tiny, align=center,rotate=45}",
        ylabel = "Nitrogen AEFM classes (\\#)",
        ymin = 0,
    },
        Plot(#
        {#
            "black",
            "fill=julia1"
        },
        Coordinates(coords_n1_1)
    ),
    Plot(#
        {#
            "black",
            "fill=julia2",
        },
        Coordinates(coords_n1_2)
    ),
    Plot(#
        {#
            "black",
            "fill=julia3"
        },
        Coordinates(coords_n1_3)
    ),
    Plot(#
        {#
            "black",
            "fill=julia4",
        },
        Coordinates(coords_n1_4)
    ),
    Legend([#
        "Src-to-sink (no revisit)", "Src-to-sink (revisit)",
        "Loop (no revisit)", "Looped (revisit)",
    ]),
)

# Carbon (scaled)
coords_c2_1 = [(i, percent_C[i][1]) for i in eachindex(percent_C)]
coords_c2_2 = [(i, percent_C[i][2]) for i in eachindex(percent_C)]
coords_c2_3 = [(i, percent_C[i][3]) for i in eachindex(percent_C)]
coords_c2_4 = [(i, percent_C[i][4]) for i in eachindex(percent_C)]

pc2 = @pgf Axis(
    {
        height = "7cm",
        width  = "7cm",
        "ybar = 0pt",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        "bar width = 6pt",
        "enlarge x limits = 0.15",
        "enlarge y limits = false",
        "legend image code/.code={\\draw [#1] (0cm,-0.1cm) rectangle (0.2cm,0.25cm); }",
        legend_style = {
            legend_columns = 2,
            "font=\\tiny",
        },
        ymax = 100,
        ymin = 0,
        xmax = length(datasets),
        xtick = [1, 2, 3, 4, 5],
        xticklabels = dataset_names,
        "xticklabel style = {font=\\tiny, align=center,rotate=45}",
        ytick = [0, 25, 50, 75, 100],
        ylabel = "Carbon AEFM classes (\\%)"
    },
        Plot(#
        {#
            "black",
            "fill=julia1"
        },
        Coordinates(coords_c2_1)
    ),
    Plot(#
        {#
            "black",
            "fill=julia2",
        },
        Coordinates(coords_c2_2)
    ),
    Plot(#
        {#
            "black",
            "fill=julia3"
        },
        Coordinates(coords_c2_3)
    ),
    Plot(#
        {#
            "black",
            "fill=julia4",
        },
        Coordinates(coords_c2_4)
    ),
    #Legend([#
        #"Src-to-sink (no revisit)", "Src-to-sink (revisit)",
        #"Loop (no revisit)", "Looped (revisit)",
    #]),
)

# Nitrogen (scaled)
coords_n2_1 = [(i, percent_N[i][1]) for i in eachindex(percent_N)]
coords_n2_2 = [(i, percent_N[i][2]) for i in eachindex(percent_N)]
coords_n2_3 = [(i, percent_N[i][3]) for i in eachindex(percent_N)]
coords_n2_4 = [(i, percent_N[i][4]) for i in eachindex(percent_N)]

pn2 = @pgf Axis(
    {
        height = "7cm",
        width  = "7cm",
        "ybar = 0pt",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        "bar width = 6pt",
        "enlarge x limits = 0.15",
        "enlarge y limits = false",
        "legend image code/.code={\\draw [#1] (0cm,-0.1cm) rectangle (0.2cm,0.25cm); }",
        legend_style = {
            legend_columns = 2,
            "font=\\tiny",
        },
        ymax = 100,
        ymin = 0,
        xmax = length(datasets),
        xtick = [1, 2, 3, 4, 5],
        xticklabels = dataset_names,
        "xticklabel style = {font=\\tiny, align=center,rotate=45}",
        ytick = [0, 25, 50, 75, 100],
        ylabel = "Nitrogen AEFM classes (\\%)"
    },
        Plot(#
        {#
            "black",
            "fill=julia1"
        },
        Coordinates(coords_n2_1)
    ),
    Plot(#
        {#
            "black",
            "fill=julia2",
        },
        Coordinates(coords_n2_2)
    ),
    Plot(#
        {#
            "black",
            "fill=julia3"
        },
        Coordinates(coords_n2_3)
    ),
    Plot(#
        {#
            "black",
            "fill=julia4",
        },
        Coordinates(coords_n2_4)
    ),
    Legend([#
        "Src-to-sink (no revisit)", "Src-to-sink (revisit)",
        "Loop (no revisit)", "Looped (revisit)",
    ]),
)
# ------------------------------------------------------------------------------

## BARPLOTS (HEPG2 FLUX) -------------------------------------------------------
flux_coords_c = [(1, pathed_c), (2, pathed_n)]
flux_coords_n = [(1, looped_c), (2, looped_n)]
frac_flux = [#
    round(pathed_c / (pathed_c + looped_c) * 100, digits = 3),
    round(pathed_n / (pathed_n + looped_n) * 100, digits = 3)
]
ph = @pgf Axis(
    {
        height = "7cm",
        width  = "5cm",
        "ybar = 0pt",
        "ymode = log",
        xmajorgrids = false,
        ymajorgrids = false,
        "xtick pos = bottom",
        "ytick pos = left",
        "bar width = 12pt",
        enlargelimits = 0.10,
        legend_style =
        {
            at = Coordinate(0.5, -0.25),
            anchor = "north",
            legend_columns = 1,
            "legend cell align={left}",
        },
        xmin = 0.5,
        xmax = 2.5,
        xtick = [1, 2],
        xticklabels = ["Carbon", "Nitrogen"],
        xlabel = "HepG2",
        ylabel = "Total atomic flux (\\%)",
    },
        Plot(#
        {#
            "black",
            "fill=cgrey"
        },
        Coordinates(flux_coords_c)
    ),
    Plot(#
        {#
            "black",
            "fill=none"
        },
        Coordinates(flux_coords_n)
    ),
    Legend(["Src-to-sink AEFMs", "Loop AEFMs"]),
    "\\node [above,font=\\tiny] at (axis cs: 0.9, 360000) {$(frac_flux[1])\\%};",
    "\\node [above,font=\\tiny] at (axis cs: 1.9, 85000) {$(frac_flux[2])\\%};",
)
# ------------------------------------------------------------------------------

## SAVE ALL PLOTS TO TEX -------------------------------------------------------
pgfsave(ex_fig_c1, pc1) # carbon
pgfsave(ex_fig_n1, pn1) # nitrogen
pgfsave(ex_fig_c2, pc2) # carbon (scaled)
pgfsave(ex_fig_n2, pn2) # nitrogen (scaled)
pgfsave(ex_fig_hepg2, ph) # HepG2 weights
# ------------------------------------------------------------------------------

## MISCELLANEOUS ANALYSES ------------------------------------------------------
met_names = load_gem_metabolite_names(load1 * load2, datasets)
resC = get_looped_revisitation_efm_indices(C)
resN = get_looped_revisitation_efm_indices(N)

loop_types_C, loop_counts_C = group_looped_revisitations(resC, met_names, C)
loop_types_N, loop_counts_N = group_looped_revisitations(resN, met_names, N)

findfirst(==(maximum(loop_counts_C[5])), loop_counts_C[5])
loop_types_C[5][17]

loop_types_C[1][22]

loop_types_N[4][2]
# ------------------------------------------------------------------------------

