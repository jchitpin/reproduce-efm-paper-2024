function manual_fix(rxns::Vector{String}, ms::Vector{String})
    ex_path = "figures/incorrect-reaction-atom-mappings/iIT341/"

    problem_ids = findall(x -> occursin("ransaminase", x), rxns)
    ms[32] # C12, C14 and C2, C8 should be swapped on the product side
    ms[35] # C5, C7 and C11, C17 should be swapped on the product side
    ms[58] # C11, C17 and C5, C7 should be swapped on the product side
    ms[125] # C2, C8 and C12, C14 should be swapped on the product side
    ms[156] # C6, C8 and C12, C19 should be swapped on the product side
    ms[186] # C4, C6 and C10, C16 should be swapped on the product side
    ms[348] # completely wrong.
    ms[374] # completely wrong.

    plot_mapped_reaction(ms[32],  ex_path * "32-"  * replace(rxns[32],  " " => "-") * ".png")
    plot_mapped_reaction(ms[35],  ex_path * "35-"  * replace(rxns[35],  " " => "-") * ".png")
    plot_mapped_reaction(ms[58],  ex_path * "58-"  * replace(rxns[58],  " " => "-") * ".png")
    plot_mapped_reaction(ms[125], ex_path * "125-" * replace(rxns[125], " " => "-") * ".png")
    plot_mapped_reaction(ms[156], ex_path * "156-" * replace(rxns[156], " " => "-") * ".png")
    plot_mapped_reaction(ms[186], ex_path * "186-" * replace(rxns[186], " " => "-") * ".png")
    plot_mapped_reaction(ms[348], ex_path * "348-" * replace(rxns[348], " " => "-") * ".png")
    plot_mapped_reaction(ms[374], ex_path * "374-" * replace(rxns[374], " " => "-") * ".png")

    ms[32] = "[O:6]=[C:5]([OH:7])[CH2:4][CH2:3][C:14](=[O:15])[C:12](=[O:11])[OH:13].[NH2:1][C@@H:2]([CH2:16][c:17]1[cH:18][cH:19][cH:20][cH:21][cH:22]1)[C:8](=[O:9])[OH:10]>>[NH3+:1][CH:14]([CH2:3][CH2:4][C:5](=[O:6])[OH:7])[C:12](=[O:11])[OH:13].[O:9]=[C:8]([OH:10])[C:2](=[O:15])[CH2:16][c:17]1[cH:18][cH:19][cH:20][cH:21][cH:22]1"
    plot_mapped_reaction(ms[32],  ex_path * "32-"  * replace(rxns[32],  " " => "-") * "-manually-corrected.png")

    ms[35] = "[O:15]=[C:14]([OH:16])[CH2:13][CH2:12][C:5](=[O:6])[C:7](=[O:8])[OH:9].[CH3:1][CH2:2][C@H:3]([CH3:4])[C@H:11]([NH2:10])[C:17](=[O:18])[OH:19]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[C:11](=[O:6])[C:17](=[O:18])[OH:19].[NH3+:10][CH:5]([CH2:12][CH2:13][C:14](=[O:15])[OH:16])[C:7](=[O:8])[OH:9]"
    plot_mapped_reaction(ms[35],  ex_path * "35-"  * replace(rxns[35],  " " => "-") * "-manually-corrected.png")

    ms[58] = "[CH3:1][CH:2]([CH3:3])[CH2:4][C@H:11]([NH2:10])[C:17](=[O:18])[OH:19].[O:15]=[C:14]([OH:16])[CH2:13][CH2:12][C:5](=[O:6])[C:7](=[O:8])[OH:9]>>[CH3:1][CH:2]([CH3:3])[CH2:4][C:11](=[O:6])[C:17](=[O:18])[OH:19].[NH3+:10][CH:5]([CH2:12][CH2:13][C:14](=[O:15])[OH:16])[C:7](=[O:8])[OH:9]"
    plot_mapped_reaction(ms[58],  ex_path * "58-"  * replace(rxns[58],  " " => "-") * "-manually-corrected.png")

    ms[125] = "[NH2:1][C@@H:2]([CH2:16][c:17]1[cH:18][cH:19][c:20]([OH:21])[cH:22][cH:23]1)[C:8](=[O:9])[OH:10].[O:6]=[C:5]([OH:7])[CH2:4][CH2:3][C:14](=[O:15])[C:12](=[O:11])[OH:13]>>[NH3+:1][CH:14]([CH2:3][CH2:4][C:5](=[O:6])[OH:7])[C:12](=[O:11])[OH:13].[O:9]=[C:8]([OH:10])[C:2](=[O:15])[CH2:16][c:17]1[cH:18][cH:19][c:20]([OH:21])[cH:22][cH:23]1"
    plot_mapped_reaction(ms[125], ex_path * "125-" * replace(rxns[125], " " => "-") * "-manually-corrected.png")

    ms[156] = "[NH3+:11][CH:12]([CH2:5][CH2:4][C:2](=[O:1])[OH:3])[C:19](=[O:20])[OH:21].[O:9]=[C:8]([OH:10])[C:6](=[O:7])[CH2:13][O:14][P:15](=[O:16])([OH:17])[OH:18]>>[O:1]=[C:2]([OH:3])[CH2:4][CH2:5][C:12](=[O:7])[C:19](=[O:20])[OH:21].[NH3+:11][C@@H:6]([CH2:13][O:14][P:15](=[O:16])([OH:17])[OH:18])[C:8](=[O:9])[OH:10]"
    plot_mapped_reaction(ms[156], ex_path * "156-" * replace(rxns[156], " " => "-") * "-manually-corrected.png")

    ms[186] = "[O:14]=[C:13]([OH:15])[CH2:12][CH2:11][C:4](=[O:5])[C:6](=[O:7])[OH:8].[CH3:1][CH:2]([CH3:3])[C@H:10]([NH2:9])[C:16](=[O:17])[OH:18]>>[CH3:1][CH:2]([CH3:3])[C:10](=[O:5])[C:16](=[O:17])[OH:18].[NH3+:9][CH:4]([CH2:11][CH2:12][C:13](=[O:14])[OH:15])[C:6](=[O:7])[OH:8]"
    plot_mapped_reaction(ms[186], ex_path * "186-" * replace(rxns[186], " " => "-") * "-manually-corrected.png")

    ms[348] = "[O:11]=[C:12]([OH:13])[CH2:3][CH2:14][C:15](=[O:16])[C:17](=[O:18])[OH:19].[NH2:1][C@@H:2]([CH2:4][C:5](=[O:6])[OH:7])[C:8](=[O:9])[OH:10]>>[NH3+:1][CH:15]([CH2:14][CH2:3][C:12](=[O:11])[OH:13])[C:17](=[O:18])[OH:19].[O:6]=[C:5]([OH:7])[CH2:4][C:2](=[O:16])[C:8](=[O:9])[OH:10]"
    plot_mapped_reaction(ms[348], ex_path * "348-" * replace(rxns[348], " " => "-") * "-manually-corrected.png")

    ms[374] = "[NH2:1][C@@H:2]([CH2:3][CH2:21][CH2:20][C@H:19]([NH:18][C:16](=[O:17])[CH2:15][CH2:14][C:12](=[O:11])[OH:13])[C:28](=[O:29])[OH:30])[C:8](=[O:9])[OH:10].[O:6]=[C:5]([OH:7])[CH2:4][CH2:22][C:23](=[O:24])[C:25](=[O:26])[OH:27]>>[NH3+:1][CH:23]([CH2:22][CH2:4][C:5](=[O:6])[OH:7])[C:25](=[O:26])[OH:27].[O:11]=[C:12]([OH:13])[CH2:14][CH2:15][C:16](=[O:17])[NH:18][C@@H:19]([CH2:20][CH2:21][CH2:3][C:2](=[O:24])[C:8](=[O:9])[OH:10])[C:28](=[O:29])[OH:30]"
    plot_mapped_reaction(ms[374], ex_path * "374-" * replace(rxns[374], " " => "-") * "-manually-corrected.png")

    return ms
end

