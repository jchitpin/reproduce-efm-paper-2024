function manual_fix(rxns::Vector{String}, ms::Vector{String})
    ex_path = "figures/incorrect-reaction-atom-mappings/iAB_RBC_283/"

    problem_ids = findall(x -> occursin("ransaminase", x), rxns)
    ms[19]

    plot_mapped_reaction(ms[19],  ex_path * "19-"  * replace(rxns[19],  " " => "-") * ".png")

    ms[19] = "[O:6]=[C:5]([OH:7])[CH2:4][CH2:3][C:14](=[O:15])[C:12](=[O:11])[OH:13].[NH2:1][C@@H:2]([CH2:16][c:17]1[cH:18][cH:19][cH:20][cH:21][cH:22]1)[C:8](=[O:9])[OH:10]>>[NH3+:1][CH:14]([CH2:3][CH2:4][C:5](=[O:6])[OH:7])[C:12](=[O:11])[OH:13].[O:9]=[C:8]([OH:10])[C:2](=[O:15])[CH2:16][c:17]1[cH:18][cH:19][cH:20][cH:21][cH:22]1"
    plot_mapped_reaction(ms[19],  ex_path * "19-"  * replace(rxns[19],  " " => "-") * "-manually-corrected.png")

    return ms
end

