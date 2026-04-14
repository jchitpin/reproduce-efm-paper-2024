## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the SBML file for the e_coli_core SBML file.
# (2) Exports the stoichiometry matrix, metabolite names, reaction names, and
#     additional meta-data.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
set = "e_coli_core"

# Import and export directories
load = "data/gems/$set/raw/"

# Import file names
im_sbml = load * "$set.xml"

# Export file names
ex_stoich                   = load * "stoich.csv"
ex_metabolites              = load * "metabolites.csv"
ex_metabolites_compartments = load * "metabolites-compartments.csv"
ex_reactions                = load * "reactions.csv"
ex_formulas                 = load * "metabolite-formulas.csv"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using SBML, CSV, Tables
# ------------------------------------------------------------------------------

## LOAD SBML FILE AND EXTRACT STOICHIOMETRY MATRIX, METABOLITES, REACTIONS -----
mdl = readSBML(im_sbml)
metabolites, reactions, S = stoichiometry_matrix(mdl)

# Dense stoichiometry matrix
S = Array(S)

# Metabolite names
mets = [mdl.species[m].name for m in metabolites]
mets = replace.(mets, " " => "_")

# Metabolite names concatenated with compartment
metsc = [#
    join([mdl.species[m].name, mdl.species[m].compartment], "_")
    for m in metabolites
]
metsc = replace.(metsc, " " => "_")

# Metabolite formulas
formulas = [mdl.species[m].formula for m in metabolites]
formulas[isnothing.(formulas)] .= ""

# Reaction names
rxns = [mdl.reactions[r].name for r in reactions]
# ------------------------------------------------------------------------------

## EXPORT STOICHIOMETRY MATRIX, METABOLITES, REACTIONS -------------------------
CSV.write(ex_stoich, Tables.table(S), header = false)
CSV.write(ex_metabolites, Tables.table(mets), header = false)
CSV.write(ex_metabolites_compartments, Tables.table(metsc), header = false)
CSV.write(ex_reactions, Tables.table(rxns), header = false)
CSV.write(ex_formulas, Tables.table(formulas), header = false)
# ------------------------------------------------------------------------------

