## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the carbon CHMC output files from the HepG2 dataset.
# (2) Converts the top 5 EFMs for a given carbon into a linked list for
#     plotting in R. 
# (3) Exports the linked lists as a text file.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
# Import and export directories
load1 = "data/gems/hepg2/"
save1 = "figures/hepg2-weights/subpanels-c-f/"

# Import filenames
im_carbon     = "atomic-efm-weights-carbon"
im_subsystems = "raw/HMRdatabase2_00.xlsx"
im_mets       = "processed/metabolites-processed.csv"

# Export filenames
ex_glutamine_1     = save1 * "glutamine/carbon-1/01-for-R/"
ex_glutamine_2_3_4 = save1 * "glutamine/carbon-2-3-4/01-for-R/"
ex_glutamine_5     = save1 * "glutamine/carbon-5/01-for-R/"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using PlotlyJS
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD HEPG2 DATASETS ---------------------------------------------------------
Ch = load_hepg2_chmcs(load1 * im_carbon)

# Check that no atomic EFM weights are negative.
for i in eachindex(Ch)
    for j in 1:length(Ch[i])
        if minimum(Ch[i][j].w) < 0
            @info("$i. $j.")
        end
    end
end

# Load metabolites
mets = vec(CSV.read(load1 * im_mets, Tables.matrix, header=false))
# ------------------------------------------------------------------------------

## BINNING ATOMIC EFMS IN SORTED ORDER -----------------------------------------
# Glutamine carbon 1
m = 14    # glutamine source metabolite
a = 1     # glutamine carbon index (corresponding to greatest number of EFMs)
s = 5     # number of atomic EFMs in each bin (bin size)
b = 1     # number of bins
t = s * b # number of EFMs to pull out

# Extract EFM sequences with the t'th greatest weight in sorted order
seq_14 = efm_state_to_metabolite_seq(Ch[m][a])
top_ids = sortperm(Ch[m][a].w, rev = true)[1:t]
subseq_14 = seq_14[top_ids]
subwgt_14 = Ch[m][a].w[top_ids]

# Reshape EFMs to a node and flow dataframe
x = [1:5]
sseq_14 = [subseq_14[x[i]] for i in 1:b]
wseq_14 = [subwgt_14[x[i]] for i in 1:b]
nodes_14_bin = efms_to_r_node_dataframe.(sseq_14, Ref(mets))
flows_14_bin = efms_to_r_flow_dataframe.(sseq_14, wseq_14)

# Export to text file for R Sankey diagrams
CSV.write(ex_glutamine_1 * "nodes.csv", nodes_14_bin[1])
CSV.write(ex_glutamine_1 * "flows.csv", flows_14_bin[1])

# Extract EFM sequences for carbon 2, 3, 4
a = 2     # glutamine carbon index
seq_14 = efm_state_to_metabolite_seq(Ch[m][a])
top_ids = sortperm(Ch[m][a].w, rev = true)[1:t]
subseq_14 = seq_14[top_ids]
subwgt_14 = Ch[m][a].w[top_ids]

# Reshape EFMs to a node and flow dataframe
x = [1:5]
sseq_14 = [subseq_14[x[i]] for i in 1:b]
wseq_14 = [subwgt_14[x[i]] for i in 1:b]
nodes_14_bin = efms_to_r_node_dataframe.(sseq_14, Ref(mets))
flows_14_bin = efms_to_r_flow_dataframe.(sseq_14, wseq_14)

# Export to text file for R Sankey diagrams
CSV.write(ex_glutamine_2_3_4 * "nodes.csv", nodes_14_bin[1])
CSV.write(ex_glutamine_2_3_4 * "flows.csv", flows_14_bin[1])

# Extract EFM sequences for carbon 5
a = 5     # glutamine carbon index
seq_14 = efm_state_to_metabolite_seq(Ch[m][a])
top_ids = sortperm(Ch[m][a].w, rev = true)[1:t]
subseq_14 = seq_14[top_ids]
subwgt_14 = Ch[m][a].w[top_ids]

# Reshape EFMs to a node and flow dataframe
x = [1:5]
sseq_14 = [subseq_14[x[i]] for i in 1:b]
wseq_14 = [subwgt_14[x[i]] for i in 1:b]
nodes_14_bin = efms_to_r_node_dataframe.(sseq_14, Ref(mets))
flows_14_bin = efms_to_r_flow_dataframe.(sseq_14, wseq_14)

# Export to text file for R Sankey diagrams
CSV.write(ex_glutamine_5 * "nodes.csv", nodes_14_bin[1])
CSV.write(ex_glutamine_5 * "flows.csv", flows_14_bin[1])
# ------------------------------------------------------------------------------

