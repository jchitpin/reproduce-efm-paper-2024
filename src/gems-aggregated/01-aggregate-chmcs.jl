## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads carbon and nitrogen CHMCs from all genome-scale models and
#     aggregates them into a single carbon and nitrogen JLD2 file.
# NOTE: Each JLD2 file is an array of an array of an array of structs.
# 1st array: Dataset (5 total)
# 2nd array: Source metabolites
# 3rd array: Atom index
# ------------------------------------------------------------------------------

## DIRECTORIES -----------------------------------------------------------------
# Dataset directories
datasets = ["e_coli_core", "hepg2", "iAB_RBC_283", "iIT341", "iSB619"]

# Import directories
load  = "data/gems/"
loadp = "/processed/"
loadc = "/atomic-efms-carbon/"
loadn = "/atomic-efms-nitrogen/"

# Import filenames
im_stoich = loadp * "stoichiometry-matrix-processed.csv"
im_smiles = loadp * "smiles-isomeric-processed.csv"
im_metabs = loadp * "metabolites-processed.csv"

# Export filenames
ex_carbon   = "data/gems-aggregated/gems-carbon.jld2";
ex_nitrogen = "data/gems-aggregated/gems-nitrogen.jld2";
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using CSV, Tables, MarkovWeightedEFMs, JLD2
dir_functions = pwd() * "/functions/gems-aggregated/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## AGGREGATE CHMCS -------------------------------------------------------------
C = Vector{Vector{Vector{CHMCAtomicSummary}}}()
N = Vector{Vector{Vector{CHMCAtomicSummary}}}()
for set in datasets
    @info set

    # Load stoichiometry matrix
    S = Int16.(CSV.read(load * set * im_stoich, Tables.matrix, header = false))

    # Load SMILES strings matching the stoichiometry rows
    smiles = vec(CSV.read(load * set * im_smiles, Tables.matrix, header = false))

    # Load metabolites
    mets = vec(CSV.read(load * set * im_metabs, Tables.matrix, header = false))

    # Source/sink metabolite indices
    srcs = get_source_metabolites(Int16.(S))
    if set == "hepg2"
        met_remove = "3-beta-D-glucuronosyl-3-beta-D-galactosyl-4-beta-D-galactosyl-O-beta-D-xylosylprotein[g]"
        idx = findfirst(==(met_remove), mets[srcs])
        popat!(srcs, idx) # remove 3-beta-D-glucuronsyl because it undergoes no internal reaction
    end

    # Total number of carbons/nitrogens in each metabolite
    amaxC = get_max_atoms(smiles, :C)
    amaxN = get_max_atoms(smiles, :N)

    # Load atomics CHMCs
    push!(C, load_all_atomic_efms(load * set * loadc, mets, srcs, amaxC))
    push!(N, load_all_atomic_efms(load * set * loadn, mets, srcs, amaxN))
end
# ------------------------------------------------------------------------------

## SAVE ATOMIC CHMC ------------------------------------------------------------
jldsave(#
    ex_carbon;
    e_coli_core = C[1],
    hepg2       = C[2],
    iAB_RBC_283 = C[3],
    iIT341      = C[4],
    iSB619      = C[5],
)

jldsave(#
    ex_nitrogen;
    e_coli_core = N[1],
    hepg2       = N[2],
    iAB_RBC_283 = N[3],
    iIT341      = N[4],
    iSB619      = N[5],
)
# ------------------------------------------------------------------------------


