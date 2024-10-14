## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Load original iIT341 dataset. Stoichiometry matrix and
#     names of the reactions and metabolites.
# (2) Uses helper functions from MarkovWeightedEFMs to pre-process dataset.
# (3) Outputs the processed metabolic model data as text files.
# ------------------------------------------------------------------------------

## DIRECTORIES -----------------------------------------------------------------
set = "iIT341"

# Import and export directories
load1 = "data/gems/$set/raw/"
save = "data/gems/$set/processed/"

# Import file names
im_stoich = load1 * "stoich.csv"
im_mets   = load1 * "metabolites-compartments.csv"
im_rxns   = load1 * "reactions.csv"
im_smiles = load1 * "smiles-isomeric.csv"

# Export file names
ex_stoich = save * "stoichiometry-matrix-processed.csv"
ex_mets   = save * "metabolites-processed.csv"
ex_rxns   = save * "reactions-processed.csv"
ex_smiles = save * "smiles-isomeric-processed.csv"
ex_reaction_smiles        = save * "reaction-smiles-processed.csv"
ex_mapped_reaction_smiles = save * "mapped-reaction-smiles-strings-processed.csv"
ex_dict_carbon   = save * "dictionary-atom-tracing-carbon.csv"
ex_dict_nitrogen = save * "dictionary-atom-tracing-nitrogen.csv"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using CSV, Tables, MarkovWeightedEFMs, BenchmarkTools
# ------------------------------------------------------------------------------

## LOAD METABOLIC MODEL --------------------------------------------------------
# Stoichiometry matrix
S = CSV.read(im_stoich, Tables.matrix, header = false)

# Metabolites
mets = vec(CSV.read(im_mets, Tables.matrix, header = false))

# Reactions
rxns = String.(vec(CSV.read(im_rxns, Tables.matrix, header = false)))
# ------------------------------------------------------------------------------

## PRE-PROCESS DATASET ---------------------------------------------------------
# (1) Identify problems with S/v inputs
errors = find_atomic_chmc_input_errors(S)
print(errors) # summary of errors associated with S

# (2) Clean S inputs
S2, mets2, rxns2 = correct_atomic_chmc_input_errors(errors, S, mets, rxns)
print(find_atomic_chmc_input_errors(S2)) # confirm errors have been fixed

# (3) Construct vector of smiles corresponding to the remaining metabolites in S
# The SMILES strings for pseudometabolites with no defined chemical structure
# are given an arbitrary SMILES of 'R' (or character that does not represent
# a periodic table element)
# SMILES strings matching S2
smiles3 = vec(CSV.read(im_smiles, Tables.matrix, header = false))

# (4) Remove pseudometabolites and reactions exceeding RXNMapper character limit
S4, mets4, rxns4, smiles4, i4 = correct_atomic_chmc_input_smiles(#
  S2, mets2, rxns2, smiles3
)
i4.dropped_rows_pseudometabolites # pseudometabolite rows removed from S2
i4.dropped_cols_pseudometabolites # pseudometabolite reactions removed from S2
i4.dropped_cols_rxnmapper_limit # reactions in S2 removed bc of RXNMapper limit
print(find_atomic_chmc_input_errors(S4)) # confirm no errors

# (5) Construct the reaction strings and map atoms via RXNMAPPER
smiles5 = canonicalize_smiles(smiles4) # smiles strings must be canonicalized!
rs5, ms5 = map_reaction_strings(S4, smiles5, rxns4, false)

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

## SAVE FINAL STOICHIOMETRY MATRIX, METABOLITES, REACTIONS, ETC. ---------------
CSV.write(ex_stoich, Tables.table(S4), header = false)
CSV.write(ex_mets, Tables.table(mets4), header = false)
CSV.write(ex_rxns, Tables.table(rxns4), header = false)
CSV.write(ex_smiles, Tables.table(smiles5), header = false)
CSV.write(ex_reaction_smiles, Tables.table(rs5), header = false)
CSV.write(ex_mapped_reaction_smiles, Tables.table(ms5), header = false)
CSV.write(ex_dict_carbon, D_C, header = false)
CSV.write(ex_dict_nitrogen, D_N, header = false)
# ------------------------------------------------------------------------------

@info("Dropped rows pseudomets: $(length(i4.dropped_rows_pseudometabolites))") # 52
@info("Dropped cols pseudomets: $(length(i4.dropped_cols_pseudometabolites))") # 93
@info("Dropped cols rxnmapper:  $(length(i4.dropped_cols_rxnmapper_limit))")   #  4

