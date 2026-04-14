## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Enumerates and saves the atomic EFMs in the iIT341 network.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set parent directory containing this repo
verbose = false # verbose mode while enumerating EFMs
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILENAMES ---------------------------------------------------
# Set parent directory
set = "iIT341"

# Import and export subdirectories
load       = "data/gems/$set/processed/"
saveC      = "data/gems/$set/atomic-efms-carbon/"
saveN      = "data/gems/$set/atomic-efms-nitrogen/"

# File names
stoich_loc = load * "stoichiometry-matrix-processed.csv"
smiles_loc = load * "smiles-isomeric-processed.csv"
mets_loc   = load * "metabolites-processed.csv"
D_C_loc    = load * "dictionary-atom-tracing-carbon.csv"
D_N_loc    = load * "dictionary-atom-tracing-nitrogen.csv"
ms_loc     = load * "mapped-reaction-smiles-strings-processed.csv"
rxns_loc   = load * "reactions-processed.csv"
time_C_loc = "data/gems/$set/efm-enumeration-times-carbon-T2.jld2"
time_N_loc = "data/gems/$set/efm-enumeration-times-nitrogen-T2.jld2"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using CSV, Tables, MarkovWeightedEFMs, JLD2
# ------------------------------------------------------------------------------

## LOAD FINAL DATA -------------------------------------------------------------
# Load stoichiometry matrix
S = Int16.(CSV.read(stoich_loc, Tables.matrix, header = false))

# Load SMILES strings matching the stoichiometry rows
smiles = vec(CSV.read(smiles_loc, Tables.matrix, header = false))

# Load metabolites
mets = vec(CSV.read(mets_loc, Tables.matrix, header = false))

# Load atom tracing dictionary
D_C = import_atom_tracing_dictionary(D_C_loc)
D_N = import_atom_tracing_dictionary(D_N_loc)

# Load mapped reaction smiles strings
ms = vec(CSV.read(#
    ms_loc, Tables.matrix, delim = ';', ignoreemptyrows = false, header = false
))
g(x) = ismissing(x) ? "" : x
ms = g.(ms)

# Load reactions
rxns = String.(vec(CSV.read(rxns_loc, Tables.matrix, header = false)))
# ------------------------------------------------------------------------------

## ENUMERATE ATOMIC EFMS -------------------------------------------------------
# Identify indices of all source metabolites and number of carbons/nitrogens
srcs = get_source_metabolites(Int16.(S))
amax_C = get_max_atoms(smiles, :C)
amax_N = get_max_atoms(smiles, :N)

# Carbon EFMs (with precompilation)
enumerate_atomic_efms(S, ms, (srcs[1], 1, :C), D_C, verbose = verbose)
ctimes = Vector{Vector{Float64}}()
for k in eachindex(srcs)
    res_dir = saveC * "$k-" * mets[srcs[k]]
    isdir(res_dir) || mkdir(res_dir)
    ctmes = Vector{Float64}()
    for l in 1:amax_C[srcs[k]]
        I = (srcs[k], l, :C)
        stime = time()
        res = enumerate_atomic_efms(S, ms, I, D_C, verbose = verbose)
        etime = time()
        push!(ctmes, etime - stime)
    end
    push!(ctimes, ctmes)
end
jldsave(time_C_loc; ctimes)

# Nitrogen EFMs
ntimes = Vector{Vector{Float64}}()
for k in eachindex(srcs)
    res_dir = saveN * "$k-" * mets[srcs[k]]
    isdir(res_dir) || mkdir(res_dir)
    ntmes = Vector{Float64}()
    for l in 1:amax_N[srcs[k]]
        I = (srcs[k], l, :N)
        stime = time()
        res = enumerate_atomic_efms(S, ms, I, D_N, verbose = verbose)
        etime = time()
        push!(ntmes, etime - stime)
    end
    push!(ntimes, ntmes)
end
jldsave(time_N_loc; ntimes)
# ------------------------------------------------------------------------------

