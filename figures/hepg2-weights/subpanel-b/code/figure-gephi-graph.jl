## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the carbon and nitrogen CHMC output files from the HepG2 dataset.
# (2) Subsets for the top X atomic EFMs for specific source metabolites.
# (3) Exports the edge lists for plotting in Gephi (or Cytoscape)
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
dataset = "hepg2"

# Import and export directories
load1 = "data/gems/hepg2/"
save1 = "figures/hepg2-weights/subpanel-b/gephi/"

# Import filenames
im_mets = "processed/metabolites-processed.csv"
im_rxns = "processed/reactions-processed.csv"
im_dictionary_c = "processed/linked-list-reactions-dict-carbon.jld2"
im_dictionary_n = "processed/linked-list-reactions-dict-nitrogen.jld2"
im_carbon_hepg2   = "atomic-efm-weights-carbon"
im_nitrogen_hepg2 = "atomic-efm-weights-nitrogen"

# Export filenames
ex_glutamine_0 = save1 * "gephi-glutamine-top-5.csv"
ex_glutamine_1 = save1 * "gephi-glutamine-top-10.csv"
ex_glutamine_2 = save1 * "gephi-glutamine-top-50.csv"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD HEPG2 DATASETS ---------------------------------------------------------
# Load metabolites
mets = vec(CSV.read(load1 * im_mets, Tables.matrix, header = false))

# Load reactions
rxns = String.(vec(CSV.read(load1 * im_rxns, Tables.matrix, header = false)))

# Load dictionaries
dict_c = jldopen(load1 * im_dictionary_c)["dict"]
dict_n = jldopen(load1 * im_dictionary_n)["dict"]

# Load atomic CHMCs
Ch = load_hepg2_chmcs(load1 * im_carbon_hepg2)
#Nh = load_hepg2_chmcs(load1 * im_nitrogen_hepg2)
# ------------------------------------------------------------------------------

## EXPORT GRAPH OF TOP 5, 10 AND 50 CARBON AEFMS FROM GLUTAMINE ----------------
y = 14
df14_0 = export_top_X_efms_as_edge_list(Ch[y], 5, dict_c[y], mets, rxns)  # top 5  each
df14_1 = export_top_X_efms_as_edge_list(Ch[y], 10, dict_c[y], mets, rxns) # top 10 each
df14_2 = export_top_X_efms_as_edge_list(Ch[y], 50, dict_c[y], mets, rxns) # top 50 each

mets_1 = unique([df14_1.Source; df14_1.Target])
filter!(!=("Sink"), mets_1)
filter!(!=("Source"), mets_1)

p14 = percentage_explained_source_flux_by_efms(Ch[y])

p14_denom = sum(sum.([p14[i].efm_fluxes for i in eachindex(p14)]))
p14_0_num = sum(sum.([p14[i].efm_fluxes[1:5] for i in 1:length(Ch[y])]))
p14_1_num = sum(sum.([p14[i].efm_fluxes[1:10] for i in 1:length(Ch[y])]))
p14_2_num = sum(sum.([p14[i].efm_fluxes[1:50] for i in 1:length(Ch[y])]))

unique([df14_0.Source; df14_0.Target]) # 22 metabolites + Source/Sink
unique([df14_1.Source; df14_1.Target]) # 32 metabolites + Source/Sink
unique([df14_2.Source; df14_2.Target]) # 102 metabolites + Source/Sink

# These numbers should roughly equal that from code in subpanel-a
p14_0 = p14_0_num / p14_denom # 44.66% explained atomic fluxes (across all 5 carbons); 24 metabolites; 31 edges (UNIQUE EDGES COLLAPSING METABOLITE-ATOMS ONTO METABOLITES).
p14_1 = p14_1_num / p14_denom # 57.25% explained atomic fluxes (across all 5 carbons); 54 metabolites; 46 edges (UNIQUE EDGES COLLAPSING METABOLITE-ATOMS ONTO METABOLITES).
p14_2 = p14_2_num / p14_denom # 86.94% explained atomic fluxes (across all 5 carbons); 124 metabolites. 146 edges (UNIQUE EDGES COLLAPSING METABOLITE-ATOMS ONTO METABOLITES.

df14_0[:, :Source] .= uppercasefirst.(df14_0[:, :Source])
df14_0[:, :Target] .= uppercasefirst.(df14_0[:, :Target])

CSV.write(ex_glutamine_0, df14_0)
CSV.write(ex_glutamine_1, df14_1)
CSV.write(ex_glutamine_2, df14_2)
# ------------------------------------------------------------------------------

## GET LENGTHS OF TOP AEFMS EXCLUDING SOURCE/SINK PSEUDONODES ------------------

gln_efm_1 = efm_state_to_metabolite_atom_seq(Ch[14][1])
gln_efm_2 = efm_state_to_metabolite_atom_seq(Ch[14][2])
gln_efm_3 = efm_state_to_metabolite_atom_seq(Ch[14][3])
gln_efm_4 = efm_state_to_metabolite_atom_seq(Ch[14][4])
gln_efm_5 = efm_state_to_metabolite_atom_seq(Ch[14][5])

gln_efm_1 = efm_state_to_metabolite_seq(Ch[14][1])
gln_efm_2 = efm_state_to_metabolite_seq(Ch[14][2])
gln_efm_3 = efm_state_to_metabolite_seq(Ch[14][3])
gln_efm_4 = efm_state_to_metabolite_seq(Ch[14][4])
gln_efm_5 = efm_state_to_metabolite_seq(Ch[14][5])

ids_1 = sortperm(Ch[14][1].w, rev=true)[1:5]
ids_2 = sortperm(Ch[14][2].w, rev=true)[1:5]
ids_3 = sortperm(Ch[14][3].w, rev=true)[1:5]
ids_4 = sortperm(Ch[14][4].w, rev=true)[1:5]
ids_5 = sortperm(Ch[14][5].w, rev=true)[1:5]

top_5_glutamine_efm_lengths = [#
    length.(gln_efm_1[ids_1]);
    length.(gln_efm_2[ids_2]);
    length.(gln_efm_3[ids_3]);
    length.(gln_efm_4[ids_4]);
    length.(gln_efm_5[ids_5]);
]
mean(top_5_glutamine_efm_lengths) # 6.0 length

top_5_glutamine_efm_lengths[1:5]   # 7, 7, 5,  3, 15
top_5_glutamine_efm_lengths[6:10]  # 7, 7, 5,  7,  3
top_5_glutamine_efm_lengths[11:15] # 7, 7, 5,  7,  3
top_5_glutamine_efm_lengths[16:20] # 7, 7, 5,  7,  3
top_5_glutamine_efm_lengths[21:25] # 5, 8, 3,  2,  8

efm_state_to_metabolite_seq(Ch[14][1])[ids_1]
efm_state_to_metabolite_atom_seq(Ch[14][1])[ids_1]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][1])[ids_1][1])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][1])[ids_1][2])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][1])[ids_1][3])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][1])[ids_1][4])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][1])[ids_1][5])[1:15]]

efm_state_to_metabolite_seq(Ch[14][2])[ids_2]
efm_state_to_metabolite_atom_seq(Ch[14][2])[ids_2]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][2])[ids_2][1])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][2])[ids_2][2])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][2])[ids_2][3])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][2])[ids_2][4])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][2])[ids_2][5])]

efm_state_to_metabolite_seq(Ch[14][3])[ids_3]
efm_state_to_metabolite_atom_seq(Ch[14][3])[ids_3]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][3])[ids_3][1])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][3])[ids_3][2])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][3])[ids_3][3])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][3])[ids_3][4])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][3])[ids_3][5])]

efm_state_to_metabolite_seq(Ch[14][4])[ids_4]
efm_state_to_metabolite_atom_seq(Ch[14][4])[ids_4]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][4])[ids_4][1])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][4])[ids_4][2])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][4])[ids_4][3])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][4])[ids_4][4])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][4])[ids_4][5])]

efm_state_to_metabolite_seq(Ch[14][5])[ids_5]
efm_state_to_metabolite_atom_seq(Ch[14][5])[ids_5]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][5])[ids_5][1])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][5])[ids_5][2])[1:8]]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][5])[ids_5][3])]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][5])[ids_5][4])[1:2]]
mets[first.(efm_state_to_metabolite_atom_seq(Ch[14][5])[ids_5][5])[1:8]]
# ------------------------------------------------------------------------------

