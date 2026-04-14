## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the carbon and nitrogen CHMC output files from the HepG2 dataset.
# (2) Constructs a Sankey diagram in the form of weighted linked lists
#     for atomic EFM bins.
# (3) Exports the linked lists and compressed linked lists.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
# Import and export directories
load1 = "data/gems/hepg2/"
save1 = "figures/hepg2-weights/subpanels-c-f/"

# Import filenames
im_carbon     = "atomic-efm-weights-carbon"
im_nitrogen   = "atomic-efm-weights-nitrogen"
im_subsystems = "raw/HMRdatabase2_00.xlsx"
im_mets       = "processed/metabolites-processed.csv"
im_rxns       = "processed/reactions-processed.csv"
im_dict_c     = "processed/linked-list-reactions-dict-carbon.jld2"

# Export filenames
ex_glutamine  = save1 * "glutamine/interactive-sankey-diagrams/"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using PlotlyJS
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD HEPG2 DATASETS ---------------------------------------------------------
Ch = load_hepg2_chmcs(load1 * im_carbon)

# Check that no atomic EFM weights are negative.
for i in 1:length(Ch)
    for j in 1:length(Ch[i])
        if minimum(Ch[i][j].w) < 0
            @info("$i. $j.")
        end
    end
end
# ------------------------------------------------------------------------------

## LOAD KEGG ANNOTATIONS FOR EACH REACTION -------------------------------------
# Load metabolites
mets = vec(CSV.read(load1 * im_mets, Tables.matrix, header=false))

# Load reactions
rxns = String.(vec(CSV.read(load1 * im_rxns, Tables.matrix, header=false)))

# Subsystems for each reaction
subsystems = kegg_subsystems_from_hmr(load1 * im_subsystems, rxns)

# Linked list dictionaries
dcr = jldopen(load1 * im_dict_c)["dict"]
# ------------------------------------------------------------------------------

## VISUALIZING TOP 5 AEFMS FOR EACH GLUTAMINE CARBON AS A SANKEY DIAGRAM -------
m = 14 # source metabolite index
a = 1  # carbon atom index
linked_lists_14, roots_14 = kegg_linked_list(Ch[m][[a]], dcr[m][[a]], subsystems, 1:5)
p = sankey_atomic_chmc(#
    linked_lists_14,
    roots_14,
    mets,
    subsystems,
    cutoff = 1e-200000,
    atom_type = :Carbon
);
savefig(p, ex_glutamine * "sankey-diagram-carbon-$a.html")

m = 14 # source metabolite index
b = 2  # carbon atom index
linked_lists_14, roots_14 = kegg_linked_list(Ch[m][[b]], dcr[m][[b]], subsystems, 1:5)
p = sankey_atomic_chmc(#
    linked_lists_14,
    roots_14,
    mets,
    subsystems,
    cutoff = 1e-200000,
    atom_type = :Carbon
);
savefig(p, ex_glutamine * "sankey-diagram-carbon-$b.html")

m = 14 # source metabolite index
c = 3  # carbon atom index
linked_lists_14, roots_14 = kegg_linked_list(Ch[m][[c]], dcr[m][[c]], subsystems, 1:5)
p = sankey_atomic_chmc(#
    linked_lists_14,
    roots_14,
    mets,
    subsystems,
    cutoff = 1e-200000,
    atom_type = :Carbon
);
savefig(p, ex_glutamine * "sankey-diagram-carbon-$c.html")

m = 14 # source metabolite index
d = 4  # carbon atom index
linked_lists_14, roots_14 = kegg_linked_list(Ch[m][[d]], dcr[m][[d]], subsystems, 1:5)
p = sankey_atomic_chmc(#
    linked_lists_14,
    roots_14,
    mets,
    subsystems,
    cutoff = 1e-200000,
    atom_type = :Carbon
);
savefig(p, ex_glutamine * "sankey-diagram-carbon-$d.html")

m = 14 # source metabolite index
e = 5  # carbon atom index
linked_lists_14, roots_14 = kegg_linked_list(Ch[m][[e]], dcr[m][[e]], subsystems, 1:5)
p = sankey_atomic_chmc(#
    linked_lists_14,
    roots_14,
    mets,
    subsystems,
    cutoff = 1e-200000,
    atom_type = :Carbon
);
savefig(p, ex_glutamine * "sankey-diagram-carbon-$e.html")
# ------------------------------------------------------------------------------

