# Get substrates, products, stoichiometries, and parameters from stoichiometry
function parse_catalyst_reactions(S::Matrix{Int16}, met_names::Vector{String})
    @assert(size(S,1) == length(met_names))

    catalyst = Vector{#
        NamedTuple{#
            (:kinetic, :substrates, :products, :substoich, :prodstoich),
            Tuple{#
                String,
                Vector{String},
                Vector{String},
                Vector{Float64},
                Vector{Float64}
            }
        }
    }(undef, size(S,2))

    k = 1
    for j in eachindex(catalyst)
        subs_idx = findall(<(0), @view S[:, j])
        prods_idx = findall(>(0), @view S[:, j])

        if isempty(subs_idx)
            catalyst[j] = (#
                kinetic = "k$k",
                substrates = ["nothing"],
                products = met_names[prods_idx],
                substoich = (S[subs_idx, j]),
                prodstoich = (S[prods_idx, j])
            )
        elseif isempty(prods_idx)
            catalyst[j] = (#
                kinetic = "k$k",
                substrates = met_names[subs_idx],
                products = ["nothing"],
                substoich = (S[subs_idx, j]),
                prodstoich = (S[prods_idx, j])
            )
        else
            catalyst[j] = (#
                kinetic = "k$k",
                substrates = met_names[subs_idx],
                products = met_names[prods_idx],
                substoich = (S[subs_idx, j]),
                prodstoich = (S[prods_idx, j])
            )
        end

        k += 1
    end
    return catalyst
end

function construct_catalyst_network(#
    S::Matrix{Int16},
    met_names::Vector{String},
    jlname::String,
    dotname::String,
    nname::String,
    ename::String
)
    @assert(size(S,1) == length(met_names))

    # @species macro doesn't work with special characters
    idx = findall(.!isnothing.(findfirst.(r"^[0-9]", met_names)))
    for i in eachindex(idx)
        met_names[idx[i]] = "_" * met_names[idx][i]
    end
    #met_names = replace.(met_names, r"\[|\]|-|\(|\)|\+|,|[0-9]|\'|:" => x -> "z")
    met_names = replace.(met_names, r"\["    => x -> "a")
    met_names = replace.(met_names, r"\]"    => x -> "b")
    met_names = replace.(met_names, r"-"     => x -> "c")
    met_names = replace.(met_names, r"\("    => x -> "d")
    met_names = replace.(met_names, r"\)"    => x -> "e")
    met_names = replace.(met_names, r"\+"    => x -> "f")
    met_names = replace.(met_names, r","     => x -> "g")
    met_names = replace.(met_names, r"[0-9]" => x -> "h")
    met_names = replace.(met_names, r"\'"    => x -> "i")
    met_names = replace.(met_names, r":"     => x -> "j")
    met_names = replace.(met_names, " "      => "_")
    R = parse_catalyst_reactions(S, met_names)
    met_names[end] = join([met_names[end], "(t)"])
    kinetics = first.(R)

    open(jlname, "w") do io
        write(io, join(["# CATALYST reaction system generator", "\n"]))
        write(io, join(["# This file was generated on: ", now(), "\n"]))
        write(io, "\n")
        write(io, "using Catalyst\n")
        write(io, "\n")
        write(io, "# Model parameters\n")
        write(io, join(["@parameters ", join(kinetics, " "), "\n"]))
        write(io, "\n")
        write(io, "# Model variables\n")
        write(io, "@variables t\n")
        write(io, "\n")
        write(io, "# Model species\n")
        write(io, join(["@species ", join(met_names, "(t) "), "\n\n"]))
        write(io, "rxs = [\n")
        nl = ",\n"
        for i in eachindex(R)
            if i == length(R)
                nl = "\n"
            end
            if R[i].substrates == ["nothing"]
                subs = "nothing"
                prods = join(["[", join(R[i].products, ", "), "]"])
                sub_stoich = "nothing"
                prod_stoich = join(["[", join(R[i].prodstoich, ", "), "]"])
            elseif R[i].products == ["nothing"]
                subs = join(["[", join(R[i].substrates, ", "), "]"])
                prods = "nothing"
                sub_stoich = join(["[", join(R[i].substoich, ", "), "]"])
                prod_stoich = "nothing"
            else
                subs = join(["[", join(R[i].substrates, ", "), "]"])
                prods = join(["[", join(R[i].products, ", "), "]"])
                sub_stoich = join(["[", join(R[i].substoich, ", "), "]"])
                prod_stoich = join(["[", join(R[i].prodstoich, ", "), "]"])
            end
            write(#
                io,
                join([#
                    "\t",
                    "Catalyst.Reaction(",
                    R[i].kinetic,
                    ", ",
                    #join(["[", join(R[i].substrates, ", "), "]"]),
                    subs,
                    ", ",
                    #join(["[", join(R[i].products, ", "), "]"]),
                    prods,
                    ", ",
                    #join(["[", join(R[i].substoich, ", "), "]"]),
                    sub_stoich,
                    ", ",
                    #join(["[", join(R[i].prodstoich, ", "), "]"]),
                    prod_stoich,
                    ")",
                    nl
                ])
            )
        end
        write(io, "];\n")
        write(io, "@named rs = ReactionSystem(rxs, t)\n")
        write(io, "g = Graph(rs)\n")
        write(io, "catalyst_to_graphviz(g, \"$dotname\")\n")
        write(io, "catalyst_to_node_edge_table(g, \"$nname\", \"$ename\")\n")
    end
end

function catalyst_to_graphviz(g::Graph, fname::String)
    # Fill edge and node lists in GraphViz DOT format
    edge_list = Vector{String}()
    node_list = Vector{String}()
    for i in eachindex(g.stmts)
        line = g.stmts[i]
        if string(typeof(line)) == "Catalyst.Node"
            push!(#
                node_list,
                join([#
                    line.name, " [label=\"", line.name, "\"];"
                ])
            )
        end
        if string(typeof(line)) == "Catalyst.Edge"
            push!(#
                edge_list,
                join([#
                    line.path[1].name,
                    " -> ",
                    line.path[2].name,
                    " [ label = \"",
                    line.attrs[:label],
                    "\" ];"
                ])
            )
        end
    end

    open(fname, "w") do io
        write(io, "digraph figure {\n")
        for e in edge_list
            write(io, join([e, "\n"]))
        end
        for n in node_list
            write(io, join([n, "\n"]))
        end
        write(io, "}\n")
    end
end

function catalyst_to_node_edge_table(g::Graph, nname::String, ename::String)

    node_ids = Vector{String}()
    node_cats = Vector{String}()
    edge_src = Vector{String}()
    edge_snc = Vector{String}()
    edge_typ = Vector{String}()
    for i in eachindex(g.stmts)
        if string(typeof(g.stmts[i])) == "Catalyst.Node"
            push!(node_ids, g.stmts[i].name)
            if contains(g.stmts[i].name, "rx_")
                push!(node_cats, "Reaction")
            else
                push!(node_cats, "Species")
            end
        end
        if string(typeof(g.stmts[i])) == "Catalyst.Edge"
            push!(edge_src, g.stmts[i].path[1].name)
            push!(edge_snc, g.stmts[i].path[2].name)
            push!(edge_typ, "Directed")
        end
    end
    ntable = Matrix{String}(undef, length(node_ids), 3)
    ntable[:,1] = node_ids
    ntable[:,2] = node_ids
    ntable[:,3] = node_cats

    etable = Matrix{String}(undef, length(edge_src), 3)
    etable[:,1] = edge_src
    etable[:,2] = edge_snc
    etable[:,3] = edge_typ

    CSV.write(nname, Tables.table(ntable), header = [:Id, :Label, :Cat])
    CSV.write(ename, Tables.table(etable), header = [:Source, :Target, :Type])
end
