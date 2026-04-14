struct Bowtie
    source::Tuple{String,Int64}
    sink::Tuple{String,Int64}
    moi::Tuple{String,Int64}
    source_to_moi::Vector{Int64}
    moi_to_sink::Vector{Int64}
    weight::Vector{Float64}
    atom_idx::Vector{Int64}
    efm_idx::Vector{Int64}
end

function compute_moi_bowties(#
    moi::Vector{Tuple{String,Int64}},
    Ch::Vector{Vector{CHMCAtomicSummary}}
)

    moi_bowties = Vector{Vector{Vector{Bowtie}}}(undef, length(moi))
    g(x) = last(x) == (0, 0) ? x[end-1] : last(x)
    h(x) = (0, 0) in x ? true : false

    # Metabolite of interest
    for i in 1:length(moi)
        moi_outer = Vector{Vector{Bowtie}}(undef, last(moi[i]))
        k = findfirst(==(first(moi[i])), mets)
        # Carbon of interest
        for j in 1:last(moi[i])
            @info("$i. $j")
            moi_inner = Vector{Bowtie}()
            for m in 1:length(Ch)
                for n in 1:length(Ch[m])
                    # (k, j) is the internal metabolite/atom state of interest
                    idx = findfirst(==((k, j)), collect(values(Ch[m][n].dmc)))
                    if !isnothing(idx)
                        efms = efm_state_to_metabolite_atom_seq(Ch[m][n])
                        for e in findall(h.(efms))
                            val = findfirst(x -> x == (k, j), efms[e])
                            if !isnothing(val)
                                src = (#
                                    mets[first(efms[e])[1]],
                                    Int64(first(efms[e])[2])
                                )
                                snk = (#
                                    mets[first(g(efms[e]))],
                                    Int64(last(g(efms[e])))
                                )
                                soi = (mets[k], j)
                                d1 = val - 1
                                d2 = length(efms[e]) - val
                                push!(#
                                    moi_inner,
                                    Bowtie(#
                                        src, snk, soi, [d1], [d2],
                                        [Ch[m][n].w[e]], [j], [e]
                                    )
                                )

                            end
                        end
                    end
                end
            end
            moi_outer[j] = moi_inner
        end
        moi_bowties[i] = moi_outer
    end
    return moi_bowties
end

function aggregate_bowties(D::Vector{Bowtie})
    E = Vector{Bowtie}()
    uniqs = Vector{Tuple{Tuple{String,Int64}, Tuple{String,Int64}}}()
    D = deepcopy(D)
    while !isempty(D)
        if (D[1].source, D[1].sink) ∉ uniqs
            push!(uniqs, (D[1].source, D[1].sink))
            push!(E, popfirst!(D))
        else
            idx = findfirst(==((D[1].source, D[1].sink)), uniqs)
            E[idx] = Bowtie(#
                E[idx].source,
                E[idx].sink,
                E[idx].moi,
                [E[idx].source_to_moi; D[1].source_to_moi],
                [E[idx].moi_to_sink; D[1].moi_to_sink],
                [E[idx].weight; D[1].weight],
                [E[idx].atom_idx; D[1].atom_idx],
                [E[idx].efm_idx; D[1].atom_idx]
            )
            popfirst!(D)
        end
    end
    return E
end

function top_x_bowtie(B::Vector{Bowtie}, x::Int64)
    sids = sortperm([b.weight for b in B], rev = true)
    return B[sids[1:x]]
end

function kegg_pathway_id_to_common_name()
    d = Dict(# 
        "map00010" => "Glycolysis / Gluconeogenesis",
        "map00020" => "Citrate cycle (TCA cycle)",
        "map00630" => "Glyoxylate and dicarboxylate metabolism",
        "map00250" => "Alanine, aspartate and glutamate metabolism",
        "map00620" => "Pyruvate metabolism",
        "map00030" => "Pentose phosphate pathway",
        "map00785" => "Lipoic acid metabolism",
        "map00260" => "Glycine, serine and threonine metabolism",
        "map00220" => "Arginine biosynthesis",
        "map00650" => "Butanoate metabolism",
        "map00524" => "Neomycin, kanamycin and gentamicin biosynthesis",
        "map00051" => "Fructose and mannose metabolism",
        "map00640" => "Propanoate metabolism",
        "map00562" => "Inositol phosphate metabolism",
        "map00561" => "Glycerolipid metabolism",
        "map00910" => "Nitrogen metabolism",
        "map00500" => "Starch and sucrose metabolism",
        "map00040" => "Pentose and glucuronate interconversions",
        "map00410" => "beta-Alanine metabolism",
        "map00270" => "Cysteine and methionine metabolism",
        "map00760" => "Nicotinate and nicotinamide metabolism",
        "map00330" => "Arginine and proline metabolism",
        "map00340" => "Histidine metabolism",
        "map00061" => "Fatty acid biosynthesis",
        "map00480" => "Glutathione metabolism",
        "map00310" => "Lysine degradation",
        "map00564" => "Glycerophospholipid metabolism",
        "map00780" => "Biotin metabolism",
        "map00062" => "Fatty acid elongation",
        "map00280" => "Valine, leucine and isoleucine degradation",
        "map00380" => "Tryptophan metabolism",
        "map00350" => "Tyrosine metabolism",
        "map00470" => "D-Amino acid metabolism",
        "map00900" => "Terpenoid backbone biosynthesis",
        "map00450" => "Selenocompound metabolism",
        "map00770" => "Pantothenate and CoA biosynthesis",
        "map00052" => "Galactose metabolism",
        "map00860" => "Porphyrin metabolism",
        "map00071" => "Fatty acid degradation",
        "map00520" => "Amino sugar and nucleotide sugar metabolism",
        "map00290" => "Valine, leucine and isoleucine biosynthesis",
        "map00053" => "Ascorbate and aldarate metabolism",
        "map00230" => "Purine metabolism",
        "map00790" => "Folate biosynthesis",
        "map01040" => "Biosynthesis of unsaturated fatty acids",
        "map00240" => "Pyrimidine metabolism",
        "map00515" => "Mannose type O-glycan biosynthesis",
        "map00120" => "Primary bile acid biosynthesis",
        "map00670" => "One carbon pool by folate",
        "map00400" => "Phenylalanine, tyrosine and tryptophan biosynthesis",
        "map00130" => "Ubiquinone and other terpenoid-quinone biosynthesis",
        "map00600" => "Sphingolipid metabolism",
        "map00563" => "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
        "map00360" => "Phenylalanine metabolism",
        "map00730" => "Thiamine metabolism",
        "map00430" => "Taurine and hypotaurine metabolism",
        "map00920" => "Sulfur metabolism",
        "map00100" => "Steroid biosynthesis",
        "map00140" => "Steroid hormone biosynthesis"
    )
    return d
end

function hierarchical_cluster(#
    M::Matrix,
    dist::Distances.UnionMetric,
    linkage::Symbol,
    met_names::Vector{String}
)
    dismat_rows = pairwise(dist, M, dims = 1)
    hcl_row = hclust(dismat_rows, linkage = linkage)
    return M[hcl_row.order, :], met_names[hcl_row.order]
end

function hierarchical_cluster_2(#
    M::Matrix,
    dist::Distances.UnionMetric,
    linkage::Symbol,
    row_names::Vector{String},
    col_names::Vector{String}
)
    dismat_rows = pairwise(dist, M, dims = 1)
    dismat_cols = pairwise(dist, M, dims = 2)
    hcl_row = hclust(dismat_rows, linkage = linkage)
    hcl_col = hclust(dismat_cols, linkage = linkage)
    return M[hcl_row.order, hcl_col.order], row_names[hcl_row.order], col_names[hcl_col.order]
end

function write_pgfplots_heatmap(#
    fdata::String,
    fname::String;
    height::Real,
    width::Real,
    xlabel::String,
    ylabel::String,
    xticklabels::Vector{String},
    yticklabels::Vector{String},
    meta_min::Real,
    meta_max::Real,
    nan_val::Real,
    colourbar_title::String
)
    xticklabels = replace.(xticklabels, "_" => " ")
    yticklabels = replace.(yticklabels, "'" => "")
    yticklabels = replace.(yticklabels, "_" => " ")

    xt = "{" * join(1:length(xticklabels), ", ") * "}"
    yt = "{" * join(1:length(yticklabels), ", ") * "}"
    xtl = "{" * join(xticklabels, ", ") * "}"
    ytl = "{" * join(yticklabels, ", ") * "}"

    io = open(fname, "w");
    # Preamble
    write(io, "\\RequirePackage{luatex85}\n")
    write(io, "\\documentclass[tikz]{standalone}\n")
    write(io, "% Default preamble\n")
    write(io, "\\usepackage{pgfplots}\n")
    write(io, "\\pgfplotsset{compat=newest}\n")
    write(io, "\\usepgfplotslibrary{groupplots}\n")
    write(io, "\\usepgfplotslibrary{polar}\n")
    write(io, "\\usepgfplotslibrary{smithchart}\n")
    write(io, "\\usepgfplotslibrary{statistics}\n")
    write(io, "\\usepgfplotslibrary{dateplot}\n")
    write(io, "\\usepgfplotslibrary{ternary}\n")
    write(io, "\\usepackage{xcolor}\n")

    # Custom colourbar
    write(io, "\\definecolor{viridis1}{HTML}{FDE725}\n")
    write(io, "\\definecolor{viridis2}{HTML}{21918C}\n")
    write(io, "\\definecolor{viridis3}{HTML}{440154}\n")

    write(io, "\\pgfplotsset{compat=newest,\n")
    write(io, "colormap={mycolourmap}{")
    write(io, join(["color(", nan_val, ")=(black)\n"]))
    write(io, join(["color(", nan_val+0.001, ")=(viridis1)\n"]))
    write(io, join(["color(", meta_min, ")=(viridis1)\n"]))
    write(io, join(["color(", (meta_max-meta_min)/2, ")=(viridis2)\n"]))
    write(io, join(["color(", meta_max, ")=(viridis3)\n"]))
    write(io, "}\n")
    write(io, "}\n")

    # Begin document with plot
    write(io, "\\begin{document}\n")
    write(io, "\\begin{tikzpicture}\n")
    write(io, "\\begin{axis}[\n")

    write(io, "view={0}{90}, colorbar,\n")
    write(io, "colormap name=mycolourmap,\n")
    write(io, "no markers, xmajorgrids={false}, ymajorgrids={false},\n")
    write(io, "enlargelimits=false,\n")
    write(io, "x tick label style={font=\\tiny},\n")
    write(io, "y tick label style={font=\\tiny},\n")
    write(io, join(["xtick=", xt, ",\n"]))
    write(io, join(["ytick=", yt, ",\n"]))
    write(io, join(["xticklabels=", xtl, ",\n"]))
    write(io, join(["yticklabels=", ytl, ",\n"]))
    write(io, join(["point meta min=", string(nan_val), ",\n"]))
    write(io, join(["point meta max=", string(meta_max), ",\n"]))
    write(io, "colorbar style={\n")
    write(io, join(["title=", colourbar_title, ",\n"]))
    write(io, "title style={align=center},\n")
    write(io, join(["ymin=", meta_min, "\n"]))
    write(io, "},\n")
    write(io, join(["height=\"", string(height), "cm\",\n"]))
    write(io, join(["width=\"", string(width), "cm\",\n"]))
    write(io, join(["xlabel={", xlabel, "},\n"]))
    write(io, join(["ylabel={", ylabel, "},\n"]))
    write(io, "]\n")

    # Plot data
    write(io, join(["\\addplot[matrix plot*,point meta=explicit] file {", fdata, "};\n"]))
    write(io, "\\end{axis}\n")
    write(io, "\\end{tikzpicture}\n")
    write(io, "\\end{document}\n")
    close(io)
end

# Not used
function histogram_pdf_curves(#
    data::Vector{Vector{Float64}},
    xmax::Int64,
    bins::Vector{Int64}
)
    hist = Vector{Vector{Vector{Real}}}(undef, length(data))
    for i in eachindex(data)
        hst = fit(Histogram, data[i], closed = :left, nbins = bins[i])
        hst = LinearAlgebra.normalize(hst, mode = :pdf)
        x = (collect(hst.edges[1])[2:end] .+ collect(hst.edges[1])[1:end - 1]) / 2
        hist[i] = [#
            [0; x; x[end] + hst.edges[1][2]; xmax],
            [0; hst.weights; 0; 0]
        ]
    end
    return hist
end

# Not used
function get_dict_gem_met_ids_to_combined_met_ids(source_mets::Vector{Vector{String}})
    # Combined/unique metabolite names
    source_mets_combined = unique(Iterators.flatten(source_mets))

    # Dictionary mapping GEM metabolite index to combined metabolite name index
    d = Vector{Dict{Int64,Int64}}(undef, length(source_mets))
    for i in eachindex(source_mets)
        dict = Dict{Int64,Int64}()
        for j in eachindex(source_mets_combined)
            k = findfirst(==(source_mets_combined[j]), source_mets[i])
            if !isnothing(k)
                dict[k] = j
            end
        end
        d[i] = dict
    end

    return d
end

# Not used
function subset_matrix_by_whitespace_rows(#
    M::Matrix,
    N::Vector{String},
    r::AbstractRange
)
    K = Vector{Matrix}(undef, length(r))
    O = Vector{Vector{String}}(undef, length(r))
    c = 0
    for k in r
        c += 1
        rids = findall(#
            [sum(M[i, :] .== -1) >= k ? true : false for i in 1:size(M, 1)]
        )
        O[c] = N[setdiff(1:length(N), rids)]
        K[c] = M[setdiff(1:size(M, 1), rids), :]
    end
    return K, O
end

# Not used
function matrix_to_xyz(M::Matrix)
    xyz = Matrix{Float64}(undef, length(M), 3)
    c = 0
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            c += 1
            xyz[c,1] = j
            xyz[c,2] = i
            xyz[c,3] = M[i,j]
        end
    end
    return xyz
end

# Not used
function export_xyz_to_pgfplots_matrix_plot(xyz::Matrix, fname::String)
    x_max = maximum(xyz[:,1])
    open(fname, "w") do f
        for i in 1:size(xyz, 1)
            write(#
                f,
                join([xyz[i, 1], " ", xyz[i, 2], " ", xyz[i, 3], "\n"])
            )
            if xyz[i, 1] == x_max
                write(f, "\n")
            end
        end
    end 
end

# Not used
function count_metabolite_occurrences_over_efms(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    function aggregate_metabolites(D::Vector{Vector{CHMCAtomicSummary}})
        res = Vector{Vector{Int16}}()
        for j in 1:length(D)
            for k in 1:length(D[j])
                e = efm_state_to_metabolite_seq(D[j][k])
                for l in 1:length(e)
                    if e[l][1] == e[l][end] # remove repeated metabolites for internally looped EFMs
                        pop!(e[l])
                    end
                end
                append!(res, e)
            end
        end
        return countmap(reduce(vcat, res))
    end

    D = Vector{Dict{Int16,Int64}}(undef, length(C))
    for i in 1:length(C)
        D[i] = aggregate_metabolites(C[i])
    end
    return D
end

# Not used
function count_metabolite_occurrences_over_efms(#
    C::Vector{CHMCAtomicSummary},
    x::Int64,
    mets::Vector{String},
    kegg::Dict{String, Union{Missing, String}},
    kegg_pathway_ids::Vector{String},
    kegg_pathway_sets::Vector{Vector{String}}
)

    kegg_mets = unique(reduce(vcat, kegg_pathway_sets))
    M = zeros(Int64, length(kegg_mets), x * length(C))

    sids = [sortperm(C[i].w, rev = true)[1:x] for i in 1:length(C)]
    efms = [efm_state_to_metabolite_seq(C[i])[sids[i]] for i in 1:length(C)]
    efms = reduce(vcat, efms)

    for j in 1:size(M, 2)
        # Convert EFM metabolite sequence to KEGG IDs
        efm_kegg_ids = [kegg[m] for m in mets[efms[j]]]
        filter!(!ismissing, efm_kegg_ids)
        if efm_kegg_ids[1] == efm_kegg_ids[end]
            pop!(efm_kegg_ids)
        end

        # Count +1 if EFM metabolite matches a KEGG pathway metabolite
        for e in efm_kegg_ids
            i = findfirst(==(e), kegg_mets)
            if !isnothing(i)
                M[i,j] += 1
            end
        end
    end

    g(x) = x[1] == x[end] ? length(x) - 1 : length(x)
    efm_lengths = g.(efms)

    return copy(transpose(transpose(M) ./ efm_lengths))
    #return M
end

# Not used
function remove_empty_rows(M::Matrix, n::Vector{String})
    idx = Vector{Int64}()
    for i in 1:size(M, 1)
        if all(M[i,:] .== -1)
            push!(idx, i)
        end
    end
    return M[setdiff(1:size(M, 1), idx), :], n[setdiff(1:length(n), idx)]
end

# Not used
function countloop(itr::Vector{Int16})
    d = Dict{Int64, Int64}()
    for val in itr
        d[val] = get(d, val, -1) + 1
    end
    return collect(filter(((k,v),) -> v > 0, d))
end

# Not used
function metabolite_loops(res::CHMCAtomicSummary)
    seq = efm_state_to_metabolite_seq(res)
    loops = Vector{Vector{Pair{Int64,Int64}}}()
    for i in eachindex(seq)
        push!(loops, countloop(seq[i]))
    end
    return loops
end

# Not used
function loops_per_metabolite(loops::Vector{Pair{Int64,Int64}})
    d = Dict{Int64,Int64}()
    for (k, v) in loops
        if !haskey(d, k)
            d[k] = v
        else
            d[k] += v
        end
    end
    return d
end

# Not used
function get_reaction(E::CHMCAtomicSummary, c::Tuple{Int64, Int64})
    m(x) = x[2]
    ii = findall(==(c[1]), first.(E.R))
    jj = findall(==(c[2]), m.(E.R))
    return E.R[ii[findfirst(x -> x in jj, ii)]].k
end

# Not used
function fraction_metabolite_efm_over_kegg(#
    C::Vector{CHMCAtomicSummary},
    x::Int64,
    mets::Vector{String},
    kegg::Dict{String, Union{Missing, String}},
    kegg_pathway_ids::Vector{String},
    kegg_pathway_sets::Vector{Vector{String}}
)
    M = zeros(Int64, length(kegg_pathway_ids), x * length(C))

    sids = [sortperm(C[i].w, rev = true)[1:x] for i in 1:length(C)]
    efms = [efm_state_to_metabolite_seq(C[i])[sids[i]] for i in 1:length(C)]
    efms = reduce(vcat, efms)

    for j in 1:size(M, 2)
        # Convert EFM metabolite sequence to KEGG IDs
        efm_kegg_ids = [kegg[m] for m in mets[efms[j]]]
        filter!(!ismissing, efm_kegg_ids)

        # Count +1 if EFM metabolite matches a KEGG pathway
        for e in efm_kegg_ids
            xs = findall(Ref(e) .∈ kegg_pathway_sets)
            M[xs,j] .+= 1
        end

    end

    return M ./ length.(kegg_pathway_sets)
end

# Not used
function is_pathed_or_looped_efms(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    pathed_or_looped_ids = Vector{Vector{Vector{Vector{Bool}}}}(undef, length(C))
    for i in 1:length(C)
        tmp1 = Vector{Vector{Vector{Bool}}}(undef, length(C[i]))
        for j in 1:length(C[i])
            tmp2 = Vector{Vector{Bool}}(undef, length(C[i][j]))
            for k in 1:length(C[i][j])
                efms = first.(C[i][j][k].e)
                tmp3 = Vector{Bool}(undef, length(efms))
                for l in 1:length(efms)
                    if efms[l][2] == 1
                        tmp3[l] = true
                    else
                        tmp3[l] = false
                    end
                end
                tmp2[k] = tmp3
            end
            tmp1[j] = tmp2
        end
        pathed_or_looped_ids[i] = tmp1
    end
    @info "1 corresponds to pathed EFMs. 0 corresponds to looped EFMs."
    return pathed_or_looped_ids
end

