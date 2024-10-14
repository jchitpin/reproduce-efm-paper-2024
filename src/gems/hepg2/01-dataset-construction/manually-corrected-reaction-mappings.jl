function manual_fix(rxns::Vector{String}, ms::Vector{String})
    ex_path = "figures/incorrect-reaction-atom-mappings/hepg2/"

    # Have to manually convert HMR reaction IDs to names to determine aminotransferase enzymes that are likely problematic
    r1  = findfirst(==("HMR_3827"), rxns)
    r2  = findfirst(==("HMR_3829"), rxns)
    r3  = findfirst(==("HMR_3903"), rxns) # correctly mapped
    r4  = findfirst(==("HMR_4109"), rxns)
    r5  = findfirst(==("HMR_3807"), rxns)
    r6  = findfirst(==("HMR_3841"), rxns) # correctly mapped
    r7  = findfirst(==("HMR_4788"), rxns)
    r8  = findfirst(==("HMR_4597"), rxns)
    r9  = findfirst(==("HMR_3747"), rxns)
    r10 = findfirst(==("HMR_3777"), rxns)
    r11 = findfirst(==("HMR_6923"), rxns)
    r12 = findfirst(==("HMR_6725"), rxns)
    r13 = findfirst(==("HMR_6768"), rxns)
    r14 = findfirst(==("HMR_4300"), rxns) # correctly mapped
    
    #r15 = findfirst(==("HMR_4437"), rxns) # histidine reaction maps correctly

    plot_mapped_reaction(ms[r1],  ex_path * "$r1-"  * replace(rxns[r1],  " " => "-") * ".png")
    plot_mapped_reaction(ms[r2],  ex_path * "$r2-"  * replace(rxns[r2],  " " => "-") * ".png")
    plot_mapped_reaction(ms[r4],  ex_path * "$r4-"  * replace(rxns[r4],  " " => "-") * ".png")
    plot_mapped_reaction(ms[r5],  ex_path * "$r5-"  * replace(rxns[r5],  " " => "-") * ".png")
    plot_mapped_reaction(ms[r7],  ex_path * "$r7-"  * replace(rxns[r7],  " " => "-") * ".png")
    plot_mapped_reaction(ms[r8],  ex_path * "$r8-"  * replace(rxns[r8],  " " => "-") * ".png")
    plot_mapped_reaction(ms[r9],  ex_path * "$r9-"  * replace(rxns[r9],  " " => "-") * ".png")
    plot_mapped_reaction(ms[r10], ex_path * "$r10-" * replace(rxns[r10], " " => "-") * ".png")
    plot_mapped_reaction(ms[r11], ex_path * "$r11-" * replace(rxns[r11], " " => "-") * ".png")
    plot_mapped_reaction(ms[r12], ex_path * "$r12-" * replace(rxns[r12], " " => "-") * ".png")
    plot_mapped_reaction(ms[r13], ex_path * "$r13-" * replace(rxns[r13], " " => "-") * ".png")

    ms[r1] = "[NH3+:11][CH:12]([CH2:5][CH2:4][C:2](=[O:1])[O-:3])[C:17](=[O:18])[O-:19].[O:15]=[C:14]([OH:16])[CH2:13][C:6](=[O:7])[C:8](=[O:9])[OH:10]>>[O:1]=[C:2]([OH:3])[CH2:4][CH2:5][C:12](=[O:7])[C:17](=[O:18])[OH:19].[NH2:11][C@@H:6]([CH2:13][C:14](=[O:15])[OH:16])[C:8](=[O:9])[OH:10]"
    plot_mapped_reaction(ms[r1],  ex_path * "$r1-"  * replace(rxns[r1],  " " => "-") * "-manually-corrected.png")

    ms[r2] = "[O:1]=[C:2]([OH:3])[CH2:4][CH2:5][C:6](=[O:7])[C:8](=[O:9])[OH:10].[NH2:11][C@@H:12]([CH2:13][C:14](=[O:15])[OH:16])[C:17](=[O:18])[OH:19]>>[NH3+:11][CH:6]([CH2:5][CH2:4][C:2](=[O:1])[O-:3])[C:8](=[O:9])[O-:10].[O:15]=[C:14]([OH:16])[CH2:13][C:12](=[O:7])[C:17](=[O:18])[OH:19]"
    plot_mapped_reaction(ms[r2],  ex_path * "$r2-"  * replace(rxns[r2],  " " => "-") * "-manually-corrected.png")

    ms[r5] = "[O:6]=[C:5]([OH:7])[CH2:4][CH2:3][C:15](=[O:16])[C:8](=[O:9])[OH:10].[NH2:1][CH2:2][CH2:14][CH2:13][C@H:12]([NH2:11])[C:17](=[O:18])[OH:19]>>[NH3+:1][CH:15]([CH2:3][CH2:4][C:5](=[O:6])[O-:7])[C:8](=[O:9])[O-:10].[NH3+:11][C@@H:12]([CH2:13][CH2:14][CH:2]=[O:16])[C:17](=[O:18])[OH:19]"
    plot_mapped_reaction(ms[r5],  ex_path * "$r5-"  * replace(rxns[r5],  " " => "-") * "-manually-corrected.png")

    ms[r7] = "[CH3:6][C@H:2]([NH2:1])[C:3](=[O:4])[OH:5].[O:8]=[CH:7][C:9](=[O:10])[OH:11]>>[NH2:1][CH2:7][C:9](=[O:10])[OH:11].[CH3:6][C:2](=[O:8])[C:3](=[O:4])[OH:5]"
    plot_mapped_reaction(ms[r7],  ex_path * "$r7-"  * replace(rxns[r7],  " " => "-") * "-manually-corrected.png")

    ms[r8] = "[O:1]=[C:2]([OH:3])[CH2:4][CH2:6][C:7](=[O:8])[C:9](=[O:10])[OH:11].[NH2:12][C@@H:13]([CH2:14][CH2:5][CH2:15][C:16](=[O:17])[OH:18])[C:19](=[O:20])[OH:21]>>[O:17]=[C:16]([OH:18])[CH2:15][CH2:5][CH2:14][C:13](=[O:8])[C:19](=[O:20])[OH:21].[NH3+:12][CH:7]([CH2:6][CH2:4][C:2](=[O:1])[O-:3])[C:9](=[O:10])[O-:11]"
    plot_mapped_reaction(ms[r8],  ex_path * "$r8-"  * replace(rxns[r8],  " " => "-") * "-manually-corrected.png")

    ms[r9] = "[O:14]=[C:13]([OH:15])[CH2:12][CH2:11][C:4](=[O:5])[C:6](=[O:7])[OH:8].[CH3:1][CH:2]([CH3:3])[C@H:10]([NH2:9])[C:16](=[O:17])[OH:18]>>[CH3:1][CH:2]([CH3:3])[C:10](=[O:5])[C:16](=[O:17])[OH:18].[NH3+:9][CH:4]([CH2:11][CH2:12][C:13](=[O:14])[O-:15])[C:6](=[O:7])[O-:8]"
    plot_mapped_reaction(ms[r9],  ex_path * "$r9-"  * replace(rxns[r9],  " " => "-") * "-manually-corrected.png")

    ms[r10] = "[O:15]=[C:14]([OH:16])[CH2:13][CH2:12][C:5](=[O:6])[C:7](=[O:8])[OH:9].[CH3:1][CH2:2][C@H:3]([CH3:4])[C@H:11]([NH2:10])[C:17](=[O:18])[OH:19]>>[CH3:1][CH2:2][CH:3]([CH3:4])[C:11](=[O:6])[C:17](=[O:18])[OH:19].[NH3+:10][CH:5]([CH2:12][CH2:13][C:14](=[O:15])[O-:16])[C:7](=[O:8])[O-:9]"
    plot_mapped_reaction(ms[r10], ex_path * "$r10-" * replace(rxns[r10], " " => "-") * "-manually-corrected.png")

    ms[r11] = "[O:15]=[C:14]([OH:16])[CH2:13][CH2:12][C:5](=[O:6])[C:7](=[O:8])[OH:9].[CH3:1][CH:2]([CH3:3])[CH2:4][C@H:11]([NH2:10])[C:17](=[O:18])[OH:19]>>[CH3:1][CH:2]([CH3:3])[CH2:4][C:11](=[O:6])[C:17](=[O:18])[OH:19].[NH3+:10][CH:5]([CH2:12][CH2:13][C:14](=[O:15])[O-:16])[C:7](=[O:8])[O-:9]"
    plot_mapped_reaction(ms[r11], ex_path * "$r11-" * replace(rxns[r11], " " => "-") * "-manually-corrected.png")

    ms[r12] = "[O:6]=[C:5]([OH:7])[CH2:4][CH2:3][C:14](=[O:15])[C:12](=[O:11])[OH:13].[NH2:1][C@@H:2]([CH2:16][c:17]1[cH:18][cH:19][cH:20][cH:21][cH:22]1)[C:8](=[O:9])[OH:10]>>[NH3+:1][CH:14]([CH2:3][CH2:4][C:5](=[O:6])[O-:7])[C:12](=[O:11])[O-:13].[O:9]=[C:8]([OH:10])[C:2](=[O:15])[CH2:16][c:17]1[cH:18][cH:19][cH:20][cH:21][cH:22]1"
    plot_mapped_reaction(ms[r12], ex_path * "$r12-" * replace(rxns[r12], " " => "-") * "-manually-corrected.png")

    ms[r13] = "[O:19]=[C:18]([OH:20])[CH2:17][CH2:16][C:4](=[O:5])[C:2](=[O:1])[OH:3].[NH2:14][C@@H:15]([CH2:6][c:7]1[cH:8][cH:9][c:10]([OH:11])[cH:12][cH:13]1)[C:21](=[O:22])[OH:23]>>[O:22]=[C:21]([OH:23])[C:15](=[O:5])[CH2:6][c:7]1[cH:8][cH:9][c:10]([OH:11])[cH:12][cH:13]1.[NH3+:14][CH:4]([CH2:16][CH2:17][C:18](=[O:19])[O-:20])[C:2](=[O:1])[O-:3]"
    plot_mapped_reaction(ms[r13], ex_path * "$r13-" * replace(rxns[r13], " " => "-") * "-manually-corrected.png")

    return ms
end
