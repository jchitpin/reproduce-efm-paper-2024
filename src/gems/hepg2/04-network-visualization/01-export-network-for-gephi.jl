## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads and modifies the metabolic network for Catalyst graph construction.
# (2) Exports a metabolic network in a Graphviz dot file.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILENAMES ---------------------------------------------------
# Set parent directory
set = "hepg2"

# Import subdirectories
load         = "data/gems/$set/processed/"
save         = "data/gems/$set/network-visualization/"

# Import file names
stoich_loc   = load * "stoichiometry-matrix-processed.csv"
mets_loc     = load * "metabolites-processed.csv"

# Export file names
catalyst_loc = pwd() * "/" * save * "catalyst-network-$set.jl"
graphviz_loc = pwd() * "/" * save * "catalyst-network-$set.dot"
node_loc = pwd() * "/" * save * "catalyst-network-$set-nodes.csv"
edge_loc = pwd() * "/" * save * "catalyst-network-$set-edges.csv"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD FINAL DATA -------------------------------------------------------------
# Load stoichiometry matrix
S = Int16.(CSV.read(stoich_loc, Tables.matrix, header = false))

# Load metabolites
mets = vec(CSV.read(mets_loc, Tables.matrix, header = false))
# ------------------------------------------------------------------------------

## CONSTRUCT CATALYST GRAPH AND EXPORT AS A GRAPHVIZ DOT FILE ------------------
# Generates a Catalyst network to graphviz/node and edge table script
construct_catalyst_network(S, mets, catalyst_loc, graphviz_loc, node_loc, edge_loc)

# Runs the script
include(catalyst_loc)
# ------------------------------------------------------------------------------

## Gephi notes: (for the node/edge table)
# 1. Load from spreadsheet the node table then append the edge table.
# 2. Partition nodes by "cat".
# 3. Set node colourscheme: #00A2E8 and #7F7F7F.
# 4. Set node size proportional to degree: minimum of 10 and maximum of 30.
# 5. Set edge colourscheme: #00A3E8.
# 6. Yifan-Hu proportional network layout with default parameters.

