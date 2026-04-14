## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the carbon and nitrogen CHMC output files from the HepG2 dataset.
# (2) Constructs a Sankey diagram in the form of weighted linked lists.
# (3) Exports the linked lists and compressed linked lists.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
dataset = "hepg2"

# Import and export directories
load1 = "data/gems/hepg2/"

# Import filenames
im_carbon     = "atomic-efm-weights-carbon"
im_nitrogen   = "atomic-efm-weights-nitrogen"
im_subsystems = "raw/HMRdatabase2_00.xlsx"
im_mets       = "processed/metabolites-processed.csv"
im_rxns       = "processed/reactions-processed.csv"
im_dict_c     = "processed/linked-list-reactions-dict-carbon.jld2"
im_dict_n     = "processed/linked-list-reactions-dict-nitrogen.jld2"

# Export filenames
ex_linked_c     = "processed/sankey-diagram-linked-lists-carbon.jld2"
ex_linked_n     = "processed/sankey-diagram-linked-lists-nitrogen.jld2"
ex_compressed_c = "processed/sankey-diagram-linked-lists-compressed-carbon.jld2"
ex_compressed_n = "processed/sankey-diagram-linked-lists-compressed-nitrogen.jld2"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using PlotlyJS
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD HEPG2 DATASETS ---------------------------------------------------------
C = load_hepg2_chmcs(load1 * im_carbon)
N = load_hepg2_chmcs(load1 * im_nitrogen)

# Check that no atomic EFM weights are negative.
for i in 1:length(C)
    for j in 1:length(C[i])
        if minimum(C[i][j].w) < 0
            @info("$i. $j.")
        end
    end
end
for i in 1:length(N)
    for j in 1:length(N[i])
        if minimum(N[i][j].w) < 0
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
dnr = jldopen(load1 * im_dict_n)["dict"]
# ------------------------------------------------------------------------------

## EXPORT FULL LINKED LIST FOR SANKEY DIAGRAMS ---------------------------------
# Linked lists for PlotlyJS Sankey diagram
linked_lists_C = Vector{Vector{KEGGArc}}(undef, length(C))
linked_lists_N = Vector{Vector{KEGGArc}}(undef, length(N))

# Linked lists with linear intermediates compressed/removed
compressed_lists_C = Vector{Vector{KEGGArc}}(undef, length(C))
compressed_lists_N = Vector{Vector{KEGGArc}}(undef, length(N))

# Indices corresponding to linear intermediates that have been removed
assembly_lists_C = Vector{Vector{Vector{Int16}}}(undef, length(C))
assembly_lists_N = Vector{Vector{Vector{Int16}}}(undef, length(N))

# Root index of each atomic CHMC
roots_C = Vector{Int16}(undef, length(C))
roots_N = Vector{Int16}(undef, length(N))

# Fill carbon linked lists
for i in eachindex(C)
    @info i
    if !isempty(C[i])
        linked_lists_C[i], roots_C[i] = kegg_linked_list(C[i], dcr[i], subsystems)
        compressed_lists_C[i], assembly_lists_C[i] = compress(linked_lists_C[i], roots_C[i])
    end
end

# Fill nitrogen linked lists
for i in eachindex(N)
    @info i
    if !isempty(N[i])
        linked_lists_N[i], roots_N[i] = kegg_linked_list(N[i], dnr[i], subsystems)
        compressed_lists_N[i], assembly_lists_N[i] = compress(linked_lists_N[i], roots_N[i])
    end
end
# ------------------------------------------------------------------------------

## SAVE LINKED LISTS -----------------------------------------------------------
# Carbon
jldopen(load1 * ex_linked_c, "w") do file
    for i in eachindex(C)
        @info i
        fname = "source_met_$i"
        if !isempty(C[i])
            file[fname] = (linked_lists_C[i], roots_C[i])
        end
    end
end
jldopen(load1 * ex_compressed_c, "w") do file
    for i in eachindex(C)
        fname = "source_met_$i"
        if !isempty(C[i])
            file[fname] = (compressed_lists_C[i], assembly_lists_C[i])
        end
    end
end

# Nitrogen
jldopen(load1 * ex_linked_n, "w") do file
    for i in eachindex(N)
        fname = "source_met_$i"
        if !isempty(N[i])
            file[fname] = (linked_lists_N[i], roots_N[i])
        end
    end
end
jldopen(load1 * ex_compressed_n, "w") do file
    for i in eachindex(N)
        fname = "source_met_$i"
        if !isempty(N[i])
            file[fname] = (compressed_lists_N[i], assembly_lists_N[i])
        end
    end
end
# ------------------------------------------------------------------------------

