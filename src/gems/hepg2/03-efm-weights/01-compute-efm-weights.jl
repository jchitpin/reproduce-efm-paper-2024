## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Computes carbon/nitrogen atomic EFM weights in the HepG2 network.
# (2) Saves the atomic EFMs as JLD2 files.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
verbose = true # verbose mode while enumerating EFMs and computing their weights
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILENAMES ---------------------------------------------------
# Set parent directory
set = "hepg2"

# Import subdirectories
load       = "data/gems/$set/processed/"
saveC      = "data/gems/$set/atomic-efm-weights-carbon/"
saveN      = "data/gems/$set/atomic-efm-weights-nitrogen/"

# File names
stoich_loc = load * "stoichiometry-matrix-processed.csv"
v_loc      = load * "flux-vector-processed.csv"
smiles_loc = load * "smiles-isomeric-processed.csv"
mets_loc   = load * "metabolites-processed.csv"
D_C_loc    = load * "dictionary-atom-tracing-carbon.csv"
D_N_loc    = load * "dictionary-atom-tracing-nitrogen.csv"
ms_loc     = load * "mapped-reaction-smiles-strings-processed.csv"
rxns_loc   = load * "reactions-processed.csv"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using CSV, Tables, MarkovWeightedEFMs, JLD2, BenchmarkTools, Dates
# ------------------------------------------------------------------------------

## LOAD FINAL DATA -------------------------------------------------------------
# Load stoichiometry matrix
S = Int16.(CSV.read(stoich_loc, Tables.matrix, header = false))

# Load fluxes
v = vec(CSV.read(v_loc, Tables.matrix, header = false))

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
met_remove = "3-beta-D-glucuronosyl-3-beta-D-galactosyl-4-beta-D-galactosyl-O-beta-D-xylosylprotein[g]"
idx = findfirst(==(met_remove), mets[srcs])
popat!(srcs, idx) # remove 3-beta-D-glucuronsyl because it undergoes no internal reaction
amax_C = get_max_atoms(smiles, :C)
amax_N = get_max_atoms(smiles, :N)

# Carbon EFMs
stime = Dates.now()
for k in eachindex(srcs)
    res_dir = saveC * "$k-" * mets[srcs[k]]
    isdir(res_dir) || mkdir(res_dir)
    for l in 1:amax_C[srcs[k]]
        @info("$k. $l")
        I = (srcs[k], l, :C)
        res = steady_state_efm_distribution(S, v, ms, I, D_C; verbose = verbose)
        jldsave(saveC * "$k-" * mets[srcs[k]] * "/atomic-efm-weights-$l.jld2"; res)
    end
end
etime = Dates.now()
cefms = etime - stime

# Nitrogen EFMs
stime = Dates.now()
for k in eachindex(srcs)
    res_dir = saveN * "$k-" * mets[srcs[k]]
    isdir(res_dir) || mkdir(res_dir)
    for l in 1:amax_N[srcs[k]]
        @info("$k. $l")
        I = (srcs[k], l, :N)
        res = steady_state_efm_distribution(S, v, ms, I, D_N; verbose = verbose)
        jldsave(saveN * "$k-" * mets[srcs[k]] * "/atomic-efm-weights-$l.jld2"; res)
    end
end
etime = Dates.now()
nefms = etime - stime
# ------------------------------------------------------------------------------

@info "It took $cefms to compute all carbon EFM weights for the hepg2 GEM."
@info "It took $nefms to compute all nitrogen EFM weights for the hepg2 GEM."

