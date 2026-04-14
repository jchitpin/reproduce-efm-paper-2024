## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the carbon and nitrogen CHMC output files from the HepG2 dataset.
# (2) Precompute dictionaries of MC state pairs and reaction index for each
#     atomic CHMC.
# (3) Export dictionaries to JLD2.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
dataset = "hepg2"

# Import and export directories
load1 = "data/gems/hepg2/"
save1 = "data/gems/hepg2/processed/"

# Import filenames
im_carbon     = "atomic-efm-weights-carbon"
im_nitrogen   = "atomic-efm-weights-nitrogen"

# Export filenames
ex_dictionary_c = "linked-list-reactions-dict-carbon.jld2"
ex_dictionary_n = "linked-list-reactions-dict-nitrogen.jld2"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD HEPG2 DATASETS ---------------------------------------------------------
C = load_hepg2_chmcs(load1 * im_carbon)
N = load_hepg2_chmcs(load1 * im_nitrogen)
# ------------------------------------------------------------------------------

## PRECOMPUTE DICTIONARY RELATING MC STATE PAIRS AND THEIR REACTION INDICES ----
linked_list_reaction_dict_c = Vector{Vector{Dict{Tuple{Int64, Int64}, Int64}}}(undef, length(C))
linked_list_reaction_dict_n = Vector{Vector{Dict{Tuple{Int64, Int64}, Int64}}}(undef, length(N))
Threads.@threads for i in eachindex(C)
    @info "$i/$(length(C))"
    linked_list_reaction_dict_c[i] = precompute_linked_list_reactions(C[i])
end

Threads.@threads for i in eachindex(N)
    @info "$i/$(length(N))"
    linked_list_reaction_dict_n[i] = precompute_linked_list_reactions(N[i])
end
# ------------------------------------------------------------------------------

## EXPORT DICTIONARIES ---------------------------------------------------------
jldsave(save1 * ex_dictionary_c; dict = linked_list_reaction_dict_c)
jldsave(save1 * ex_dictionary_n; dict = linked_list_reaction_dict_n)
# ------------------------------------------------------------------------------

