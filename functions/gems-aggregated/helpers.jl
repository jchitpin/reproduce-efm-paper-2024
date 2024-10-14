# Load all atomic CHMC JLD2 files
function load_all_atomic_efms(#
    dir::String,
    mets::Vector{String},
    src_mets::Vector{Int64},
    amax::Vector{Int64}
)
    C = Vector{Vector{CHMCAtomicSummary}}()
    for k in eachindex(src_mets)
        tmp = Vector{CHMCAtomicSummary}()
        for l in 1:amax[src_mets[k]]
            fname = dir * "$k-" * mets[src_mets[k]] * "/atomic-efms-$l.jld2"
            push!(tmp, jldopen(fname, "r")["res"])
        end
        push!(C, tmp)
    end

  return C
end




