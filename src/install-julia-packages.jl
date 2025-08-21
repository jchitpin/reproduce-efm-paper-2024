## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Installs all Julia packages necessary to run the scripts.
#     Note: GLMakie may return an error if OpenGL is not installed.
#     See https://github.com/JuliaPlots/GLMakie.jl for troubleshooting.
#     It is only required for one optional plotting function in
#     MarkovWeightedEFMs.jl that is not used in this repository.
# (2) An additional step is required to install MarkovWeightedEFMs.jl
#     to enumerate/compute atomic EFM weights.
#     The requirement is installing RXNMapper.
#
## RXNMapper installation guide (pacman is the package manager for Arch Linux)
# * Tested with Python version 3.10
# $ pacman -S rust # Rust seemed to be required but your mileage may vary
# $ pip install virtualenv
# $ virtualenv --python="/usr/bin/python310" "/path/to/new/virtualenv1/"
# Also tested that python3.11 works
# $ source virtualenv1/bin/activate
# $ pip install rxnmapper
# $ pip install rdkit
# $ pip install requests
# $ pip install tdqm
# $ pip install bs4
# $ pip install CTSgetPy
# (virtualenv1) $ julia
# julia> using Pkg, PyCall
# julia> ENV["PYTHON"] = joinpath(ENV["VIRTUAL_ENV"], "bin", "python")
# julia> Pkg.build("PyCall") 
# ------------------------------------------------------------------------------

## INSTALLING MarkovWeightedEFMs.jl -------------------------------------------
using Pkg
#Pkg.add(url="https://github.com/jchitpin/MarkovWeightedEFMs.jl.git")
# ------------------------------------------------------------------------------

## INSTALLING PACKAGE DEPENDENCIES ---------------------------------------------
Pkg.add([#
    "BenchmarkTools",   # benchmarking run times
    "CairoMakie",       # plotting backend for heatmaps
    "Clustering",       # clustering for heatmaps
    "Distances",        # distance metrics for heatmaps
    "Distributions",    # boxplot jittering
    "GLM",              # linear regression
    "GLMakie",          # plotting backend for MarkovWeightedEFMs
    "GR",               # plotting backend for density shade plots
    "GRUtils",          # plotting backend for density shade plots
    "JLD2",             # file compression
    "LaTeXStrings",     # LaTeX support for Julia plots
    "LinearRegression", # linear regression
    "NaturalSort",      # natural sorting
    "PGFPlotsX",        # plotting backend for non-heatmaps/Sankey diagrams
    "PlotlyJS",         # plotting backend for Sankey diagrams
    "Plots",            # plotting package
    "Random",           # seed for jittering
    "SBML",             # loading SBML-specific files
    "StatsBase",        # common statistical calculations
    "CSV",              # CSV utilities
    "Catalyst",         # for constructing and exporting graphs for Gephi
    "ChunkSplitters",   
    "DataFrames",       # dataframe support
    "Dates",            # date utilities
    "ExtendableSparse",
    "GeometryBasics",
    "GraphMakie",
    "Graphs",
    "LinearAlgebra",    # linear algebra utilities
    "LinearSolve",
    "Luxor",
    "Makie",
    "MolecularGraph",   # checking molecular formulas from SMILES strings
    "NetworkLayout",
    "Printf",
    "ProgressMeter",
    "PubChemCrawler",   # molecular formula utilities
    "RDKitMinimalLib",
    "Reexport",
    "SparseArrays",
    "Statistics",
    "Tables",           # CSV helper
    "FilePaths",
    "FileIO"
])
# ------------------------------------------------------------------------------

Pkg.instantiate()
