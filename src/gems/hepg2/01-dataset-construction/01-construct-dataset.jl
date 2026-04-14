## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Load original/raw PNAS HepG2 dataset. Stoichiometry matrix, flux vector,
#     and names of the reactions and metabolites.
# (2) Uses helper functions from MarkovWeightedEFMs to pre-process dataset.
# (3) Outputs the processed metabolic model data as text files.
# ------------------------------------------------------------------------------

## DIRECTORIES -----------------------------------------------------------------
set = "hepg2"

# Import and export directories
load = "data/gems/$set/raw/"
save = "data/gems/$set/processed/"

# Import file names
im_stoich = load * "pnas-fig3-stoich.csv"
im_flxs   = load * "pnas-fig3-fluxes.csv"
im_mets   = load * "pnas-fig3-metabolites.csv"
im_rxns   = load * "pnas-fig3-reactions.csv"
im_smiles = load * "hmr-pnas-dataset-isomeric.csv"

# Export file names
ex_stoich = save * "stoichiometry-matrix-processed.csv"
ex_flxs   = save * "flux-vector-processed.csv"
ex_mets   = save * "metabolites-processed.csv"
ex_rxns   = save * "reactions-processed.csv"
ex_smiles = save * "smiles-isomeric-processed.csv"
ex_reaction_smiles = save * "reaction-smiles-processed.csv"
ex_mapped_reaction_smiles = save * "mapped-reaction-smiles-strings-processed.csv"
ex_dict_carbon = save * "dictionary-atom-tracing-carbon.csv"
ex_dict_nitrogen = save * "dictionary-atom-tracing-nitrogen.csv"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using FilePaths, FileIO
using CSV, Tables, MarkovWeightedEFMs, BenchmarkTools, DataFrames
# ------------------------------------------------------------------------------

## LOAD METABOLIC MODEL --------------------------------------------------------
# These data were exported from MATLAB file figure3_metaboliteBalances.m from
# https://github.com/SysBioChalmers/LiverCellMetabolismSimulation
# Note only reactions/stoichiometry columns with non-zero fluxes are kept

# Stoichiometry matrix
S = CSV.read(im_stoich, Tables.matrix, header = false)

# Steady state fluxes
v = vec(CSV.read(im_flxs, Tables.matrix, header = true))

# Metabolites
mets = vec(CSV.read(im_mets, Tables.matrix, header = true))

# Reactions
rxns = String.(vec(CSV.read(im_rxns, Tables.matrix, header = true)))
# ------------------------------------------------------------------------------

## PRE-PROCESS DATASET ---------------------------------------------------------
# (1) Identify problems with S/v inputs
errors = find_atomic_chmc_input_errors(S, v)
print(errors) # summary of errors associated with S/v

# (2) Clean S/v inputs
S2, v2, mets2, rxns2 = correct_atomic_chmc_input_errors(errors, S, v, mets, rxns)
print(find_atomic_chmc_input_errors(S2, v2)) # confirm errors have been fixed

# (3) Construct vector of smiles corresponding to the remaining metabolites in S
df = CSV.read(im_smiles, DataFrame, header = true)
smiles3 = String.(strip.(df[!, :Smiles]))

# (4) Remove pseudometabolites and reactions exceeding RXNMapper character limit
S4, v4, mets4, rxns4, smiles4, i4 = correct_atomic_chmc_input_smiles(#
  S2, v2, mets2, rxns2, smiles3
)
i4.dropped_rows_pseudometabolites # 33 pseudometabolite rows removed from S2
i4.dropped_cols_pseudometabolites # 46 pseudometabolite reactions removed from S2
i4.dropped_cols_rxnmapper_limit # 3 reactions in S2 removed bc of RXNMapper limit
print(find_atomic_chmc_input_errors(S4, v4)) # confirm no errors

# (5) Construct the reaction strings and map atoms via RXNMAPPER
smiles5 = canonicalize_smiles(smiles4) # smiles strings must be canonicalized!
rs5, ms5 = map_reaction_strings(S4, smiles5, rxns, false)

# Manually fix incorrect atom mappings involving transaminases
include("./manually-corrected-reaction-mappings.jl")
ms5 = manual_fix(rxns4, ms5)

# (6) Precompute atom tracing dictionary (for carbons)
amax = get_max_atoms(smiles5, :C)
D_C = precompute_atom_tracing_dictionary(S4, ms5, amax, :C)

# (7) Precompute atom tracing dictionary (for nitrogens)
amax = get_max_atoms(smiles5, :N)
D_N = precompute_atom_tracing_dictionary(S4, ms5, amax, :N)
# ------------------------------------------------------------------------------

## SAVE FINAL STOICHIOMETRY MATRIX, METABOLITES, REACTIONS, FLUX VECTOR, ETC. --
CSV.write(ex_stoich, Tables.table(S4), header = false)
CSV.write(ex_flxs, Tables.table(v4), header = false)
CSV.write(ex_mets, Tables.table(mets4), header = false)
CSV.write(ex_rxns, Tables.table(rxns4), header = false)
CSV.write(ex_smiles, Tables.table(smiles5), header = false)
CSV.write(ex_reaction_smiles, Tables.table(rs5), header = false)
CSV.write(ex_mapped_reaction_smiles, Tables.table(ms5), header = false)
CSV.write(ex_dict_carbon, D_C, header = false)
CSV.write(ex_dict_nitrogen, D_N, header = false)
# ------------------------------------------------------------------------------

@info("Dropped rows pseudomets: $(length(i4.dropped_rows_pseudometabolites))") # 33
@info("Dropped cols pseudomets: $(length(i4.dropped_cols_pseudometabolites))") # 46
@info("Dropped cols rxnmapper:  $(length(i4.dropped_cols_rxnmapper_limit))")   # 3
# ------------------------------------------------------------------------------

## Inspecting some mapped reaction SMILES strings
#using ElectronDisplay

# HMR_4964 (citrate[m] to citrate[c])
#id1 = findfirst(==("HMR_4964"), rxns4)
#plot_mapped_reaction(ms5[id1], view = true)

# HMR_4454 (citrate[c] to isocitrate[c] results in the central carbon popped off as CO2)
# This mapping needs to be corrected...
# 6 to 11 is incorrect (citrate to isocitrate)
# 5 to 9 is incorrect (citrate to isocitrate)
#id2 = findfirst(==("HMR_4454"), rxns4)
#plot_mapped_reaction(ms5[id2], view = true)

# HMR_0710 (isocitrate[c] to AKG[c])
#id3 = findfirst(==("HMR_0710"), rxns4)
#plot_mapped_reaction(ms5[id3], view = true)

# HMR_3829 (AKG[c] to glutamate[c])
#id4 = findfirst(==("HMR_3829"), rxns4)
#plot_mapped_reaction(ms5[id4], view = true)

# HMR_3825 (glutamate[c] to glutamate[m])
#id5 = findfirst(==("HMR_3825"), rxns4)
#plot_mapped_reaction(ms5[id5], view = true)

# HMR_3827 (glutamate[m] to AKG[m])
#id6 = findfirst(==("HMR_3827"), rxns4)
#plot_mapped_reaction(ms5[id6], view = true)

# HMR_5297 (AKG[m] to succinyl-coA[m] involves CO2 production)
#id7 = findfirst(==("HMR_5297"), rxns4)
#plot_mapped_reaction(ms5[id7], view = true)

# HMR_4152 (succinyl-coA[m] to succinate[m])
#id8 = findfirst(==("HMR_4152"), rxns4)
#plot_mapped_reaction(ms5[id8], view = true)

# HMR_4652 (succinate[m] to fumarate[m])
#id9 = findfirst(==("HMR_4652"), rxns4)
#plot_mapped_reaction(ms5[id9], view = true)

# HMR_4410 (fumarate[m] to malate[m])
#id10 = findfirst(==("HMR_4410"), rxns4)
#plot_mapped_reaction(ms5[id10], view = true)

# HMR_4141 (malate[m] to oaa[m])
#id11 = findfirst(==("HMR_4141"), rxns4)
#plot_mapped_reaction(ms5[id11], view = true)

# HMR_4145 (oaa[m] to citrate[m])
#id12 = findfirst(==("HMR_4145"), rxns4)
#plot_mapped_reaction(ms5[id12], view = true)

