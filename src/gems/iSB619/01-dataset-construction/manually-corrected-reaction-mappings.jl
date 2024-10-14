function manual_fix(rxns::Vector{String}, ms::Vector{String})
    ex_path = "figures/incorrect-reaction-atom-mappings/iSB619/"

    problem_ids = findall(x -> occursin("ransaminase", x), rxns)
    ms[6]
    ms[23]
    ms[171]
    ms[284]
    ms[285]
    ms[290]
    ms[340]
    ms[346]
    ms[370]
    ms[507]
    ms[561]

    plot_mapped_reaction(ms[6],   ex_path * "6-"   * replace(rxns[6],   " " => "-") * ".png")
    plot_mapped_reaction(ms[23],  ex_path * "23-"  * replace(rxns[23],  " " => "-") * ".png")
    plot_mapped_reaction(ms[171], ex_path * "171-" * replace(rxns[171], " " => "-") * ".png")
    plot_mapped_reaction(ms[284], ex_path * "284-" * replace(rxns[284], " " => "-") * ".png")
    plot_mapped_reaction(ms[285], ex_path * "285-" * replace(rxns[285], " " => "-") * ".png")
    plot_mapped_reaction(ms[290], ex_path * "290-" * replace(rxns[290], " " => "-") * ".png")
    plot_mapped_reaction(ms[340], ex_path * "340-" * replace(rxns[340], " " => "-") * ".png")
    plot_mapped_reaction(ms[346], ex_path * "346-" * replace(rxns[346], " " => "-") * ".png")
    plot_mapped_reaction(ms[370], ex_path * "370-" * replace(rxns[370], " " => "-") * ".png")
    plot_mapped_reaction(ms[507], ex_path * "507-" * replace(rxns[507], " " => "-") * ".png")
    plot_mapped_reaction(ms[561], ex_path * "561-" * replace(rxns[561], " " => "-") * ".png")

    ms[6] = "[NH2:13][C@@H:14]([CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1)[C:20](=[O:21])[OH:22].[O:18]=[C:17]([OH:19])[CH2:16][CH2:15][C:4](=[O:5])[C:2](=[O:1])[OH:3]>>[O:21]=[C:20]([OH:22])[C:14](=[O:5])[CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1.[NH3+:13][CH:4]([CH2:15][CH2:16][C:17](=[O:18])[OH:19])[C:2](=[O:1])[OH:3]"
    plot_mapped_reaction(ms[6],   ex_path * "6-"   * replace(rxns[6],   " " => "-") * "-manually-corrected.png")

    ms[23] = "[CH3:1][CH:2]([CH3:3])[CH2:4][C@H:11]([NH2:10])[C:17](=[O:18])[OH:19].[O:15]=[C:14]([OH:16])[CH2:13][CH2:12][C:5](=[O:6])[C:7](=[O:8])[OH:9]>>[CH3:1][CH:2]([CH3:3])[CH2:4][C:11](=[O:6])[C:17](=[O:18])[OH:19].[NH3+:10][CH:5]([CH2:12][CH2:13][C:14](=[O:15])[OH:16])[C:7](=[O:8])[OH:9]"
    plot_mapped_reaction(ms[23],  ex_path * "23-"  * replace(rxns[23],  " " => "-") * "-manually-corrected.png")

    ms[171] = "[CH3:1][CH2:2][C@H:3]([CH3:4])[C@H:11]([NH2:10])[C:17](=[O:18])[OH:19].[O:15]=[C:14]([OH:16])[CH2:13][CH2:12][C:5](=[O:6])[C:7](=[O:8])[OH:9]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[C:11](=[O:6])[C:17](=[O:18])[OH:19].[NH3+:10][CH:5]([CH2:12][CH2:13][C:14](=[O:15])[OH:16])[C:7](=[O:8])[OH:9]"
    plot_mapped_reaction(ms[171], ex_path * "171-" * replace(rxns[171], " " => "-") * "-manually-corrected.png")

    ms[284] = "[O:13]=[C:12]([OH:14])[CH2:11][CH2:10][C:2](=[O:1])[C:15](=[O:16])[OH:17].[NH3+:8][CH2:9][CH2:3][CH2:4][C:5](=[O:6])[OH:7]>>[O:1]=[CH:9][CH2:3][CH2:4][C:5](=[O:6])[OH:7].[NH3+:8][CH:2]([CH2:10][CH2:11][C:12](=[O:13])[OH:14])[C:15](=[O:16])[OH:17]"
    plot_mapped_reaction(ms[284], ex_path * "284-" * replace(rxns[284], " " => "-") * "-manually-corrected.png")

    ms[285] = "[NH2:21][C@@H:22]([CH2:23][CH2:11][CH2:10][C@H:9]([NH:8][C:6](=[O:7])[CH2:5][CH2:4][C:2](=[O:1])[OH:3])[C:18](=[O:19])[OH:20])[C:28](=[O:29])[OH:30].[O:26]=[C:25]([OH:27])[CH2:24][CH2:12][C:13](=[O:14])[C:15](=[O:16])[OH:17]>>[O:1]=[C:2]([OH:3])[CH2:4][CH2:5][C:6](=[O:7])[NH:8][C@@H:9]([CH2:10][CH2:11][CH2:23][C:22](=[O:14])[C:28](=[O:29])[OH:30])[C:18](=[O:19])[OH:20].[NH3+:21][CH:13]([CH2:12][CH2:24][C:25](=[O:26])[OH:27])[C:15](=[O:16])[OH:17]"
    plot_mapped_reaction(ms[285], ex_path * "285-" * replace(rxns[285], " " => "-") * "-manually-corrected.png")

    ms[290] = "[O:21]=[C:20]([CH2:3][O:4][P:5](=[O:6])([OH:7])[OH:8])[CH2:9][c:10]1[cH:11][nH:12][cH:13][n:14]1.[NH3+:1][CH:2]([CH2:19][CH2:18][C:16](=[O:15])[OH:17])[C:22](=[O:23])[OH:24]>>[NH3+:1][C@H:20]([CH2:3][O:4][P:5](=[O:6])([OH:7])[OH:8])[CH2:9][c:10]1[cH:11][nH:12][cH:13][n:14]1.[O:15]=[C:16]([OH:17])[CH2:18][CH2:19][C:2](=[O:21])[C:22](=[O:23])[OH:24]"
    plot_mapped_reaction(ms[290], ex_path * "290-" * replace(rxns[290], " " => "-") * "-manually-corrected.png")

    ms[340] = "[CH3:11][C@@H:2]([NH2:1])[C:8](=[O:9])[OH:10].[O:6]=[C:5]([OH:7])[CH2:4][CH2:3][C:12](=[O:13])[C:14](=[O:15])[OH:16]>>[NH2:1][C@H:12]([CH2:3][CH2:4][C:5](=[O:6])[OH:7])[C:14](=[O:15])[OH:16].[CH3:11][C:2](=[O:13])[C:8](=[O:9])[OH:10]"
    plot_mapped_reaction(ms[340], ex_path * "340-" * replace(rxns[340], " " => "-") * "-manually-corrected.png")

    ms[346] = "[NH2:14][C@@H:15]([CH2:6][c:7]1[cH:8][cH:9][c:10]([OH:11])[cH:12][cH:13]1)[C:21](=[O:22])[OH:23].[O:19]=[C:18]([OH:20])[CH2:17][CH2:16][C:4](=[O:5])[C:2](=[O:1])[OH:3]>>[O:22]=[C:21]([OH:23])[C:15](=[O:5])[CH2:6][c:7]1[cH:8][cH:9][c:10]([OH:11])[cH:12][cH:13]1.[NH3+:14][CH:4]([CH2:16][CH2:17][C:18](=[O:19])[OH:20])[C:2](=[O:1])[OH:3]"
    plot_mapped_reaction(ms[346], ex_path * "346-" * replace(rxns[346], " " => "-") * "-manually-corrected.png")

    ms[370] = "[CH3:1][CH:2]([CH3:3])[C@H:10]([NH2:9])[C:16](=[O:17])[OH:18].[O:14]=[C:13]([OH:15])[CH2:12][CH2:11][C:4](=[O:5])[C:6](=[O:7])[OH:8]>>[CH3:1][CH:2]([CH3:3])[C:10](=[O:5])[C:16](=[O:17])[OH:18].[NH3+:9][CH:4]([CH2:11][CH2:12][C:13](=[O:14])[OH:15])[C:6](=[O:7])[OH:8]"
    plot_mapped_reaction(ms[370], ex_path * "370-" * replace(rxns[370], " " => "-") * "-manually-corrected.png")

    ms[507] = "[O:9]=[C:8]([OH:10])[C:6](=[O:7])[CH2:13][O:14][P:15](=[O:16])([OH:17])[OH:18].[NH3+:11][CH:12]([CH2:5][CH2:4][C:2](=[O:1])[OH:3])[C:19](=[O:20])[OH:21]>>[O:1]=[C:2]([OH:3])[CH2:4][CH2:5][C:12](=[O:7])[C:19](=[O:20])[OH:21].[NH3+:11][C@@H:6]([CH2:13][O:14][P:15](=[O:16])([OH:17])[OH:18])[C:8](=[O:9])[OH:10]"
    plot_mapped_reaction(ms[507], ex_path * "507-" * replace(rxns[507], " " => "-") * "-manually-corrected.png")

    ms[561] = "[NH2:10][C@@H:11]([CH2:4][C:2](=[O:1])[OH:3])[C:17](=[O:18])[OH:19].[O:15]=[C:14]([OH:16])[CH2:13][CH2:12][C:5](=[O:6])[C:7](=[O:8])[OH:9]>>[O:1]=[C:2]([OH:3])[CH2:4][C:11](=[O:6])[C:17](=[O:18])[OH:19].[NH3+:10][CH:5]([CH2:12][CH2:13][C:14](=[O:15])[OH:16])[C:7](=[O:8])[OH:9]"
    plot_mapped_reaction(ms[561], ex_path * "561-" * replace(rxns[561], " " => "-") * "-manually-corrected.png")

    return ms
end

