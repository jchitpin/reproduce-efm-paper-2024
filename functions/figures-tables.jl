using MarkovWeightedEFMs
using JLD2, PGFPlotsX, StatsBase, LinearAlgebra, BigCombinatorics
using CSV, DataFrames, Tables, XLSX, NaturalSort, LaTeXStrings
using Distances, Clustering, IterTools, Distributions, Random, GLM
using Catalyst, Dates

# For benchmarking, structural, state space
function load_aggregated_chmcs(path::String, datasets::Vector{String})
    C = Vector{Vector{Vector{CHMCAtomicSummary}}}(undef, length(datasets))
    jldopen(path, "r") do file
        for i in eachindex(datasets)
            @info datasets[i]
            C[i] = file[datasets[i]]
        end
    end
    return C
end

# For benchmarking
function load_atomic_chmc_times(dir1::String, datasets::Vector{String}, dir2::String, str::String)
    Cy = Vector{Vector{Vector{Float64}}}(undef, length(datasets))
    for i in 1:length(datasets)
        Cy[i] = jldopen(dir1 * datasets[i] * dir2)[str]
    end
    return Cy
end

# For benchmarking
function get_efm_total_vs_time_benchmark(#
    C::Vector{Vector{Vector{CHMCAtomicSummary}}},
    Cy::Vector{Vector{Vector{Float64}}}
)
    x_coords = Vector{Vector{Int64}}(undef, length(C))
    y_coords = Vector{Vector{Float64}}(undef, length(C))
    for i in 1:length(C)
        x_tmp = Vector{Int64}()
        y_tmp = Vector{Float64}()
        for j in 1:length(C[i])
            for k in 1:length(C[i][j])
                efm_total = length(C[i][j][k].e)
                push!(x_tmp, efm_total)
                push!(y_tmp, Cy[i][j][k])
            end
        end
        x_coords[i] = x_tmp
        y_coords[i] = y_tmp
    end
    return x_coords, y_coords
end

# For structural
function load_gem_metabolite_names(path::String, datasets::Vector{String})
    M = Vector{Vector{String}}(undef, length(datasets))
    for i in eachindex(datasets)
        M[i] = vec(#
            CSV.read(#
                path * datasets[i] * "/processed/metabolites-processed.csv",
                Tables.matrix,
                header = false,
            )
        )
    end
    return M
end

# For structural
function load_chmc_source_metabolite_indices(path::String, datasets::Vector{String})
    M = Vector{Vector{Int64}}(undef, length(datasets))
    for i in eachindex(datasets)
        S = Int16.(#
            CSV.read(#
                path * datasets[i] * "/processed/stoichiometry-matrix-processed.csv",
                Tables.matrix,
                header = false
            )
        )
        M[i] = get_source_metabolites(S)
    end
    return M
end

# For structural
function get_efm_length(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    g(x) = x[2] == 1 ? length(x) - 1 : length(x)
    G = Vector{Vector{Int64}}(undef, length(C))
    for i in eachindex(C)
        G[i] = reduce(vcat, [g.(first.(k.e)) for j in C[i] for k in j])
    end
    return G
end

# Used across all figures
function pgfplotsx_preamble()
    push!(#
          PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\usetikzlibrary{patterns}"""
         )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\usepackage{xcolor}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{cred}{HTML}{ED1C24}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{cgrey}{HTML}{7F7F7F}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{cblue}{HTML}{00A2E8}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{cgreen}{HTML}{22B14C}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{cyellow}{HTML}{FFF200}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{corange}{HTML}{EA7904}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{cpurple}{HTML}{9100FC}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{julia1}{HTML}{1F77B4}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{julia2}{HTML}{FF7F0E}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{julia3}{HTML}{2CA02C}"""
    )
    push!(#
        PGFPlotsX.CUSTOM_PREAMBLE,
        raw"""\definecolor{julia4}{HTML}{D62728}"""
    )
end

# For structural
function hepg2_names_to_bigg_dict()
    dict = Dict(#
        "valine[s]"                                      => "L-Valine_e",
        "urea[s]"                                        => "Urea_CH4N2O_e",
        "tyrosine[s]"                                    => "L-Tyrosine_e",
        "tryptophan[s]"                                  => "L-Tryptophan_e",
        "threonine[s]"                                   => "L-Threonine_e",
        "sulfate[s]"                                     => "Sulfate_e",
        "serine[s]"                                      => "L-Serine_e",
        "pyruvate[s]"                                    => "Pyruvate_e",
        "proline[s]"                                     => "L-Proline_e",
        "Pi[s]"                                          => "Phosphate_e",
        "phenylalanine[s]"                               => "L-Phenylalanine_e",
        "ornithine[s]"                                   => "Ornithine_e",
        "ornithine[s]"                                   => "Ornithine_e",
        "O2[s]"                                          => "O2_O2_e",
        "nicotinamide[s]"                                => "Nicotinamide_e",
        "NH3[s]"                                         => "Ammonium_e",
        "Na+[s]"                                         => "Sodium_e",
        "methionine[s]"                                  => "L-Methionine_e",
        "lysine[s]"                                      => "L-Lysine_e",
        "L-lactate[s]"                                   => "L-Lactate_e",
        "leucine[s]"                                     => "L-Leucine_e",
        "isoleucine[s]"                                  => "L-Isoleucine_e",
        "histidine[s]"                                   => "L-Histidine_e",
        "HCO3-[s]"                                       => "Bicarbonate_e",
        "H2O[s]"                                         => "H2O_H2O_e",
        "H+[s]"                                          => "H+_e",
        "glycine[s]"                                     => "Glycine_e",
        "glutamine[s]"                                   => "L-Glutamine_e",
        "glutamate[s]"                                   => "L-Glutamate_e",
        "glucose[s]"                                     => "D-Glucose_e",
        "formate[s]"                                     => "Formate_e",
        "CO2[s]"                                         => "CO2_CO2_e",
        "choline[s]"                                     => "Choline_C5H14NO_e",
        "chloride[s]"                                    => "Chloride_e",
        "aspartate[s]"                                   => "L-Aspartate_e",
        "arginine[s]"                                    => "L-Arginine_e",
        "alanine[s]"                                     => "L-Alanine_e",
        "CoA[c]"                                         => "Coenzyme_A_c",
        "glycogen[c]"                                    => "Glycogen_C6H10O5_c",
        "xanthosine-5-phosphate[c]"                      => "Xanthosine_5'-phosphate_c",
        "valine[c]"                                      => "L-Valine_c",
        "urocanate[c]"                                   => "Urocanate_C6H5N2O2_c",
        "urea[c]"                                        => "Urea_CH4N2O_c",
        "UMP[c]"                                         => "UMP_C9H11N2O9P_c",
        "UDP-N-acetylglucosamine[c]"                     => "UDP-N-acetyl-D-glucosamine_c",
        "UDP-glucuronate[c]"                             => "UDP-D-glucuronate_c",
        "UDP-glucose[c]"                                 => "UDPglucose_c",
        "UDP[c]"                                         => "UDP_C9H11N2O12P2_c",
        "tyrosine[c]"                                    => "L-Tyrosine_c",
        "tryptophan[c]"                                  => "L-Tryptophan_c",
        "threonine[c]"                                   => "L-Threonine_c",
        "sulfite[c]"                                     => "Sulfite_c",
        "sulfate[c]"                                     => "Sulfate_c",
        "squalene[c]"                                    => "Squalene_C30H50_c",
        "sn-glycerol-3-phosphate[c]"                     => "Glycerol_3-phosphate_c",
        "serine[c]"                                      => "L-Serine_c",
        "sedoheptulose-7-phosphate[c]"                   => "Sedoheptulose_7-phosphate_c",
        "SAM[c]"                                         => "S-Adenosyl-L-methionine_c",
        "SAICAR[c]"                                      => "(S)-2-[5-Amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxamido]succinate_c",
        "SAH[c]"                                         => "S-Adenosyl-L-homocysteine_c",
        "ribulose-5-phosphate[c]"                        => "D-Ribulose_5-phosphate_c",
        "ribose-5-phosphate[c]"                          => "Alpha-D-Ribose_5-phosphate_c",
        "ribose-1-phosphate[c]"                          => "Alpha-D-Ribose_1-phosphate_c",
        "pyruvate[c]"                                    => "Pyruvate_c",
        "proline[c]"                                     => "L-Proline_c",
        "PPi[c]"                                         => "Diphosphate_c",
        "Pi[c]"                                          => "Phosphate_c",
        "phenylpyruvate[c]"                              => "Phenylpyruvate_c",
        "phenylalanine[c]"                               => "L-Phenylalanine_c",
        "PEP[c]"                                         => "Phosphoenolpyruvate_c",
        "pantothenate[c]"                                => "(R)-Pantothenate_c",
        "palmitate[c]"                                   => "Hexadecanoate_(n-C16:0)_c",
        "orotidine-5-phosphate[c]"                       => "Orotidine_5'-phosphate_c",
        "orotate[c]"                                     => "Orotate_C5H3N2O4_c",
        "ornithine[c]"                                   => "Ornithine_c",
        "oleate[c]"                                      => "Octadecenoate_(n-C18:1)_c",
        "OAA[c]"                                         => "Oxaloacetate_c",
        "O2[c]"                                          => "O2_O2_c",
        "nicotinamide[c]"                                => "Nicotinamide_c",
        "NH3[c]"                                         => "Ammonium_c",
        "N-formyl-GAR[c]"                                => "N2-Formyl-N1-(5-phospho-D-ribosyl)glycinamide_c",
        "N-formimino-L-glutamate[c]"                     => "N-Formimidoyl-L-glutamate_c",
        "N-carbamoyl-L-aspartate[c]"                     => "N-Carbamoyl-L-aspartate_c",
        "NADPH[c]"                                       => "Nicotinamide_adenine_dinucleotide_phosphate_-_reduced_c",
        "NADP+[c]"                                       => "Nicotinamide_adenine_dinucleotide_phosphate_c",
        "NADH[c]"                                        => "Nicotinamide_adenine_dinucleotide_-_reduced_c",
        "NAD+[c]"                                        => "Nicotinamide_adenine_dinucleotide_c",
        "N-acetylglucosamine-6-phosphate[c]"             => "N-Acetyl-D-glucosamine_6-phosphate_c",
        "N-acetylglucosamine-1-phosphate[c]"             => "N-Acetyl-D-glucosamine_1-phosphate_c",
        "Na+[c]"                                         => "Sodium_c",
        "myristic acid[c]"                               => "Tetradecanoate_(n-C14:0)_c",
        "methionine[c]"                                  => "L-Methionine_c",
        "malonyl-CoA[c]"                                 => "Malonyl_CoA_C24H33N7O19P3S_c",
        "malate[c]"                                      => "L-Malate_c",
        "lysine[c]"                                      => "L-Lysine_c",
        "L-lactate[c]"                                   => "L-Lactate_c",
        "leucine[c]"                                     => "L-Leucine_c",
        "L-cystathionine[c]"                             => "L-Cystathionine_c",
        "K+[c]"                                          => "Potassium_c",
        "isopentenyl-pPP[c]"                             => "Isopentenyl_diphosphate_c",
        "isoleucine[c]"                                  => "L-Isoleucine_c",
        "isocitrate[c]"                                  => "Isocitrate_c",
        "inositol[c]"                                    => "Myo-Inositol_c",
        "homocysteine[c]"                                => "L-Homocysteine_c",
        "HMG-CoA[c]"                                     => "Hydroxymethylglutaryl_CoA_C27H39N7O20P3S_c",
        "histidine[c]"                                   => "L-Histidine_c",
        "HCO3-[c]"                                       => "Bicarbonate_c",
        "H2O[c]"                                         => "H2O_H2O_c",
        "H+[c]"                                          => "H+_c",
        "GTP[c]"                                         => "GTP_C10H12N5O14P3_c",
        "GSSG[c]"                                        => "Oxidized_glutathione_c",
        "GSH[c]"                                         => "Reduced_glutathione_c",
        "GMP[c]"                                         => "GMP_C10H12N5O8P_c",
        "glycolate[c]"                                   => "Glycolate_C2H3O3_c",
        "glycine[c]"                                     => "Glycine_c",
        "glutamine[c]"                                   => "L-Glutamine_c",
        "glutamate[c]"                                   => "D-Glutamate_c",
        "glucose-6-phosphate[c]"                         => "D-Glucose_6-phosphate_c",
        "glucose-1-phosphate[c]"                         => "D-Glucose_1-phosphate_c",
        "glucose[c]"                                     => "D-Glucose_c"                              ,
        "glucosamine-6-phosphate[c]"                     => "D-Glucosamine_6-phosphate_c",
        "glucono-1,5-lactone-6-phosphate[c]"             => "6-phospho-D-glucono-1,5-lactone_c",
        "geranyl-PP[c]"                                  => "Geranyl_diphosphate_c",
        "GDP[c]"                                         => "GDP_C10H12N5O11P2_c",
        "GAR[c]"                                         => "N1-(5-Phospho-D-ribosyl)glycinamide_c",
        "GAP[c]"                                         => "Glyceraldehyde_3-phosphate_c",
        "gamma-glutamyl-cysteine[c]"                     => "Gamma-L-Glutamyl-L-cysteine_c",
        "fumarate[c]"                                    => "Fumarate_c",
        "fructose-6-phosphate[c]"                        => "D-Fructose_6-phosphate_c",
        "fructose-1,6-bisphosphate[c]"                   => "D-Fructose_1,6-bisphosphate_c",
        "formate[c]"                                     => "Formate_c",
        "farnesyl-PP[c]"                                 => "Farnesyl_diphosphate_c",
        "FAICAR[c]"                                      => "5-Formamido-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide_c",
        "ethanolamine[c]"                                => "Ethanolamine_c",
        "erythrose-4-phosphate[c]"                       => "D-Erythrose_4-phosphate_c",
        "D-xylulose-5-phosphate[c]"                      => "D-Xylulose_5-phosphate_c",
        "D-xylulose[c]"                                  => "D-Xylulose_c",
        "dUDP[c]"                                        => "DUDP_C9H11N2O11P2_c",
        "dTMP[c]"                                        => "DTMP_C10H13N2O8P_c",
        "dimethylallyl-PP[c]"                            => "Dimethylallyl_diphosphate_c",
        "dihydrofolate[c]"                               => "7,8-Dihydrofolate_c",
        "DHAP[c]"                                        => "Dihydroxyacetone_c",
        "dGMP[c]"                                        => "DGMP_C10H12N5O7P_c",
        "dGDP[c]"                                        => "DGDP_C10H12N5O10P2_c",
        "dephospho-CoA[c]"                               => "Dephospho-CoA_c",
        "dCMP[c]"                                        => "DCMP_C9H12N3O7P_c",
        "dCDP[c]"                                        => "DCDP_C9H12N3O10P2_c",
        "dAMP[c]"                                        => "DAMP_C10H12N5O6P_c",
        "dADP[c]"                                        => "DADP_C10H12N5O9P2_c",
        "D-4-phosphopantothenate[c]"                     => "D-4'-Phosphopantothenate_c",
        "cysteine[c]"                                    => "L-Cysteine_c",
        "CTP[c]"                                         => "CTP_C9H12N3O14P3_c",
        "CMP[c]"                                         => "CMP_C9H12N3O8P_c",
        "citrate[c]"                                     => "Citrate_c",
        "choline[c]"                                     => "Choline_C5H14NO_c",
        "chloride[c]"                                    => "Chloride_c",
        "CDP[c]"                                         => "CDP_C9H12N3O11P2_c",
        "carbamoyl-phosphate[c]"                         => "Carbamoyl_phosphate_c",
        "ATP[c]"                                         => "ATP_C10H12N5O13P3_c",
        "aspartate[c]"                                   => "L-Aspartate_c",
        "asparagine[c]"                                  => "L-Asparagine_c",
        "arginine[c]"                                    => "L-Arginine_c",
        "AMP[c]"                                         => "AMP_C10H12N5O7P_c",
        "alanine[c]"                                     => "L-Alanine_c",
        "AKG[c]"                                         => "2-Oxoglutarate_c",
        "ADP[c]"                                         => "ADP_C10H12N5O10P2_c",
        "adenylyl sulfate[c]"                            => "Adenosine_5'-phosphosulfate_c",
        "adenosine[c]"                                   => "Adenosine_c",
        "acetyl-CoA[c]"                                  => "Acetyl-CoA_c",
        "acetoacetyl-CoA[c]"                             => "Acetoacetyl-CoA_c",
        "acetoacetate[c]"                                => "Acetoacetate_c",
        "acetate[c]"                                     => "Acetate_c",
        "6-phospho-D-gluconate[c]"                       => "6-Phospho-D-gluconate_c",
        "5-phosphoribosylamine[c]"                       => "5-Phospho-beta-D-ribosylamine_c",
        "5-phosphoribosyl-4-carboxy-5-aminoimidazole[c]" => "5-phosphoribosyl-5-carboxyaminoimidazole_c",
        "5-methyl-THF[c]"                                => "5-Methyltetrahydrofolate_c",
        "5,10-methylene-THF[c]"                          => "5,10-Methylenetetrahydrofolate_c",
        "5,10-methenyl-THF[c]"                           => "5,10-Methenyltetrahydrofolate_c",
        "4-methyl-2-oxopentanoate[c]"                    => "4-Methyl-2-oxopentanoate_c",
        "4-imidazolone-5-propanoate[c]"                  => "4-Imidazolone-5-propanoate_c",
        "4-hydroxyphenylpyruvate[c]"                     => "3-(4-Hydroxyphenyl)pyruvate_c",
        "3-phosphoserine[c]"                             => "O-Phospho-L-serine_c",
        "3-phospho-D-glycerate[c]"                       => "3-Phospho-D-glycerate_c",
        "3-hydroxyoctadecanoyl-CoA[c]"                   => "(R)-3-Hydroxyoctadecanoyl-[acyl-carrier_protein]_c",
        "2-phospho-D-glycerate[c]"                       => "D-Glycerate_2-phosphate_c",
        "2,3-bisphospho-D-glycerate[c]"                  => "2_3_Disphospho_D_glycerate_C3H3O10P2_c",
        "1-pyrroline-5-carboxylate[m]"                   => "1-Pyrroline-5-carboxylate_m",
        "10-formyl-THF[c]"                               => "10-Formyltetrahydrofolate_c",
        "(S)-dihydroorotate[c]"                          => "(S)-Dihydroorotate_c",
        "(R)-mevalonate[c]"                              => "_R__Mevalonate_C6H11O4_c",
        "(R)-5-phosphomevalonate[c]"                     => "_R__5_Phosphomevalonate_C6H10O7P_c",
        "(R)-4-phosphopantothenoyl-cysteine[c]"          => "N-((R)-4-Phosphopantothenoyl)-L-cysteine_c",
        "(2E)-octadecenoyl-CoA[c]"                       => "Octadecenoyl-CoA_(n-C18:1CoA)_c",
    )
    return dict
end

# For structural
function hepg2_names_to_bigg(mets::Vector{String})
    # hepg2 replacement names
    dict = hepg2_names_to_bigg_dict()

    for i in eachindex(mets)
        if haskey(dict, mets[i])
            mets[i] = dict[mets[i]]
        end
    end
    return mets
end

# For structural
function get_number_of_atomic_efms(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    # G[i] corresponds to GEM
    # G[i][j] corresponds to source metabolite
    # G[i][j][k] corresponds to atom position
    G = Vector{Vector{Vector{Int64}}}(undef, length(C))
    for i in 1:length(C)
        tmp2 = Vector{Vector{Int64}}(undef, length(C[i]))
        for j in 1:length(C[i])
            tmp1 = Vector{Int64}(undef, length(C[i][j]))
            for k in 1:length(C[i][j])
                tmp1[k] = length(C[i][j][k].e)
            end
            tmp2[j] = tmp1
        end
        G[i] = tmp2
    end
    return G
end

# For structural, state space
function get_atomic_chmc_state_space(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    # G[i] corresponds to GEM
    # G[i][j] corresponds to source metabolite
    # G[i][j][k] corresponds to atom position
    G = Vector{Vector{Vector{Int64}}}(undef, length(C))
    for i in 1:length(C)
        tmp2 = Vector{Vector{Int64}}(undef, length(C[i]))
        for j in 1:length(C[i])
            tmp1 = Vector{Int64}(undef, length(C[i][j]))
            for k in 1:length(C[i][j])
                tmp1[k] = length(C[i][j][k].dchmc)
            end
            tmp2[j] = tmp1
        end
        G[i] = tmp2
    end
    return G
end

# For state space
function get_atomic_mc_state_space(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    # G[i] corresponds to GEM
    # G[i][j] corresponds to source metabolite
    # G[i][j][k] corresponds to atom position
    G = Vector{Vector{Vector{Int64}}}(undef, length(C))
    for i in 1:length(C)
        tmp2 = Vector{Vector{Int64}}(undef, length(C[i]))
        for j in 1:length(C[i])
            tmp1 = Vector{Int64}(undef, length(C[i][j]))
            for k in 1:length(C[i][j])
                tmp1[k] = length(C[i][j][k].dmc)
            end
            tmp2[j] = tmp1
        end
        G[i] = tmp2
    end
    return G
end

# For structural
function get_number_of_unique_metabolite_atom_states(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    # G[i] corresponds to GEM
    # G[i][j] corresponds to source metabolite
    # G[i][j][k] corresponds to atom position
    G = Vector{Vector{Vector{Int64}}}(undef, length(C))
    for i in 1:length(C)
        tmp2 = Vector{Vector{Int64}}(undef, length(C[i]))
        for j in 1:length(C[i])
            tmp1 = Vector{Int64}(undef, length(C[i][j]))
            for k in 1:length(C[i][j])
                tmp1[k] = length(unique(Iterators.flatten(keys(C[i][j][k].dchmc))))
            end
            tmp2[j] = tmp1
        end
        G[i] = tmp2
    end
    return G
end

# For structural
function get_number_of_unique_metabolites(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    # G[i] corresponds to GEM
    # G[i][j] corresponds to source metabolite
    # G[i][j][k] corresponds to atom position
    G = Vector{Vector{Vector{Int64}}}(undef, length(C))
    for i in 1:length(C)
        tmp2 = Vector{Vector{Int64}}(undef, length(C[i]))
        for j in 1:length(C[i])
            tmp1 = Vector{Int64}(undef, length(C[i][j]))
            for k in 1:length(C[i][j])
                tmp1[k] = length(unique(first.(#
                    [#
                        C[i][j][k].dmc[l]
                        for l in unique(Iterators.flatten(keys(C[i][j][k].dchmc)))
                    ]
                )))
            end
            tmp2[j] = tmp1
        end
        G[i] = tmp2
    end
    return G
end

# For structural
function get_max_over_min_ratio_within_source(G::Vector{Vector{Vector{Int64}}})
    H = Vector{#
        Vector{@NamedTuple{min::Int64,max::Int64,ratio::Float64}}
    }(undef, length(G))
    for i in 1:length(G)
        tmp = Vector{#
            @NamedTuple{min::Int64, max::Int64, ratio::Float64}
        }(undef, length(G[i]))
        for j in 1:length(G[i])
            if !isempty(G[i][j]) && length(G[i][j]) > 1
                g_min = minimum(G[i][j])
                g_max = maximum(G[i][j])
                ratio = g_max / g_min
                tmp[j] = (min = g_min, max = g_max, ratio = ratio)
            else
                tmp[j] = (min = 0, max = 0, ratio = 0)
            end
        end
        H[i] = tmp
    end
    return H
end

# For pairwise
function get_all_max_over_min_ratio_within_source(G::Vector{Vector{Vector{Int64}}})
    H = Vector{#
        Vector{@NamedTuple{min::Int64,max::Int64,ratio::Float64}}
    }(undef, length(G))
    for i in 1:length(G)
        tmp = Vector{#
            @NamedTuple{min::Int64, max::Int64, ratio::Float64}
        }()
        for j in 1:length(G[i])
            if !isempty(G[i][j]) && length(G[i][j]) > 1
                g_min = minimum(G[i][j])
                ratios = G[i][j] ./ g_min
                for k in 1:length(ratios)
                    push!(tmp, (min = g_min, max = G[i][j][k], ratio = ratios[k])) 
                end
            else
                push!(tmp, (min = 0, max = 0, ratio = 0))
            end
        end
        H[i] = tmp
    end
    return H
end

# For cumulative mass flow, structural
function load_hepg2_chmcs(pth::String)
    C = Vector{Vector{CHMCAtomicSummary}}()
    dirs = readdir(pth)
    dirs = sort(dirs, lt = natural)
    for k in 1:length(dirs)
        tmp = Vector{CHMCAtomicSummary}()
        files = readdir(pth * "/" * dirs[k])
        files = sort(files, lt = natural)
        for l in 1:length(files)
            push!(tmp, jldopen(pth * "/" * dirs[k] * "/" * files[l])["res"])
        end
        push!(C, tmp)
    end
    return C
end

# For structural
function get_source_met_summary_dataframe(#
    D::Vector{Vector{@NamedTuple{min::Int64, max::Int64, ratio::Float64}}},
    source_names::Vector{Vector{String}},
    jit::Real
)
    df = DataFrame(x = Real[], y = Real[], name = String[], z = Real[])
    h(x) = x.ratio
    for i in eachindex(D)
        for j in 1:length(D[i])
            push!(df, [i, h(D[i][j]), source_names[i][j], jitter(i, jit)])
        end
    end
    return df
end

# For structural
function classify_atomic_efms(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    V = Vector{@NamedTuple{#
        n_pathed_no_rec::Int64, n_pathed_rec::Int64,
        n_looped_no_rec::Int64, n_looped_rec::Int64
    }}(undef, length(C))

    for i in 1:length(C)
        tmp = [0, 0, 0, 0]
        for j in 1:length(C[i])
            for k in 1:length(C[i][j])
                efms = efm_state_to_metabolite_atom_seq(C[i][j][k])
                c = 0
                for e in efms
                    c += 1
                    x = Int64.(first.(unique(e)))
                    has_recurrence = length(x) == length(unique(x)) ? false : true
                    if e[end] == (0, 0) && has_recurrence == false
                        tmp[1] += 1
                    elseif e[end] == (0, 0) && has_recurrence == true
                        tmp[2] += 1
                    elseif e[end] != (0, 0) && has_recurrence == false
                        tmp[3] += 1
                    elseif e[end] != (0, 0) && has_recurrence == true
                        push!(x, c)
                        tmp[4] += 1
                    else
                        @error("Something went wrong. Missing case.")
                    end
                end
            end
        end
        V[i] = (#
            n_pathed_no_rec = tmp[1],
            n_pathed_rec = tmp[2],
            n_looped_no_rec = tmp[3],
            n_looped_rec = tmp[4]
        )
    end
    return V
end

# For structural
function get_looped_revisitation_efm_indices(C::Vector{Vector{Vector{CHMCAtomicSummary}}})
    tmp2 = Any[]
    for i in 1:length(C)
        tmp = Any[]
        for j in 1:length(C[i])
            for k in 1:length(C[i][j])
                efms = efm_state_to_metabolite_atom_seq(C[i][j][k])
                c = 0
                for e in efms
                    c += 1
                    x = Int64.(first.(unique(e)))
                    has_recurrence = length(x) == length(unique(x)) ? false : true
                    if e[end] != (0, 0) && has_recurrence == true
                        push!(tmp, [i, j, k, c])
                    end
                end
            end
        end
        push!(tmp2, tmp)
    end
    return tmp2
end

# For structural
function group_looped_revisitations(res::Vector{Any}, met_names::Vector{Vector{String}}, A::Vector{Vector{Vector{CHMCAtomicSummary}}})
    unique_type = Vector{Vector{Vector{String}}}()
    unique_count = Vector{Vector{Int64}}()
    for resi in eachindex(res)
        unique_types = Vector{Vector{String}}()
        unique_counts = Vector{Int64}()
        for r in res[resi]
            x = efm_state_to_metabolite_atom_seq(A[r[1]][r[2]][r[3]])[r[4]]
            y = met_names[resi][first.(x)]

            if !isnothing(findfirst(==(y), unique_types))
                id = findfirst(==(y), unique_types)
                unique_counts[id] += 1
            else
                push!(unique_types, y)
                push!(unique_counts, 1)
            end
        end
        push!(unique_type, unique_types)
        push!(unique_count, unique_counts)
    end
    return unique_type, unique_count
end

# For structural
function sum_pathed_and_looped_efm_fluxes(D::Vector{Vector{CHMCAtomicSummary}})
    weight_pathed = 0.0
    weight_looped = 0.0
    for i in 1:length(D)
        for j in 1:length(D[i])
            efms = efm_state_to_metabolite_atom_seq(D[i][j])
            for k in eachindex(efms)
                x = Int64.(first.(unique(efms[k])))
                #has_recurrence = length(x) == length(unique(x)) ? false : true
                has_recurrence = x[end] == 0 ? false : true
                if has_recurrence == false
                    weight_pathed += D[i][j].w[k] * length(efms[k]) # number of reactions and nodes are equal for either pathed or looped EFMs
                else
                    weight_looped += D[i][j].w[k] * length(efms[k])
                end
            end
        end
    end
    return weight_pathed, weight_looped
end

# For cumulative mass flow
function output_percentage_explained_source_flux_by_efms(#
    C::Vector{CHMCAtomicSummary}
)
    res = percentage_explained_source_flux_by_efms(C)
    g(x) = x.efm_fluxes
    yc = reduce(vcat, cumsum.(g.(res)) ./ sum.(g.(res)))

    h(x) = x > 1 ? log10(0) : log10(1 - x)
    yclog = h.(yc)

    df = DataFrame(#
        x = reduce(vcat, [1:l for l in length.(g.(res))]),
        y = reduce(vcat, g.(res)),
        yc = yc,
        yclog = yclog,
        atom = reduce(vcat, [repeat([i], length(g(res[i]))) for i in 1:length(res)])
    )

    return df
end

# For HepG2
function efm_state_to_metabolite_atom_seq(res::CHMCAtomicSummary)
    e = first.(res.e)
    f = Vector{Vector{Tuple{Int16,Int16}}}(undef, length(e))
    for i in 1:length(e)
        if e[i][2] == 1
            k = 2
        else
            k = 1
        end
        f[i] = [res.dmc[e[i][l]] for l in k:length(e[i])]
    end
    return f
end

# For HepG2
function efm_state_to_metabolite_seq(res::CHMCAtomicSummary)
    e = first.(res.e)
    f = Vector{Vector{Int16}}(undef, length(e))
    for i in 1:length(e)
        if e[i][2] == 1
            k = 2
        else
            k = 1
        end
        f[i] = [first(res.dmc[e[i][l]]) for l in k:length(e[i])]
        filter!(!iszero, f[i])
    end
    return f
end

# For structural, state space
function jitter(x::Real, jit::Real=0.2)
    d = Normal(1, jit)
    return x - (1 .- rand(d, 1))[1]
end

# For cumulative mass flow
function percentage_explained_source_flux_by_efms(
    C::Vector{CHMCAtomicSummary}
)
    # Loop over atom type index
    res = Vector{@NamedTuple{sids::Vector{Int64}, efm_fluxes::Vector{Float64}}}(undef, length(C))
    for i in 1:length(C)

        # EFM indices ordered by greatest explained flux
        efm_lengths = length.(efm_state_to_metabolite_atom_seq(C[i]))
        efm_fluxes = C[i].w .* efm_lengths
        sorted_ids = sortperm(efm_fluxes, rev = true)

        res[i] = (sids = sorted_ids, efm_fluxes = efm_fluxes[sorted_ids])
    end

    return res
end

# For HepG2 linked lists
struct KEGGArc
    head::Int16
    tail::Int16
    weight::Float64
    subsystems::Vector{Int16}
end

# For HepG2 linked lists
function kegg_linked_list(#
    D::Vector{CHMCAtomicSummary},
    dc::Vector{Dict{Tuple{Int64, Int64}, Int64}},
    subsystems::Vector{String}
)
    # Root/source metabolite index
    root = first(D[1].dmc[1]) # MC state 1 is always the root

    # Collect atomic CHMC EFMs involving source metabolite as linked list
    arcs = Vector{KEGGArc}()

    for i in 1:length(D)
        # Identify MC state associated with the pseudo-metabolite state (0, 0)
        pm_idx = findfirst(==((0, 0)), collect(values(D[i].dmc)))
        pm_mc = collect(keys(D[i].dmc))[pm_idx]

        # Only generate linked lists for source-to-sink atomic EFMs
        for j in 1:length(D[i].e)
            if first(D[i].e[j])[end] == pm_mc
                # Pair up metabolite index sequence
                states = first.([D[i].dmc[e] for e in first(D[i].e[j])])
                cs_mets = [pair for pair in IterTools.partition(states, 2, 1)]

                # Identify reaction sequence
                cs_mc = [pair for pair in IterTools.partition(first(D[i].e[j]), 2, 1)]
                cs_rxns = [dc[i][c] for c in cs_mc]

                # Identify and remove sequences involving the same metabolites
                dups = nonunique(first.(cs_mets))
                while !isempty(dups)
                    ids = findall(==(dups[1]), first.(cs_mets))
                    deleteat!(cs_mets, UnitRange((ids[1]), ids[2]-1))
                    deleteat!(cs_rxns, UnitRange((ids[1]), ids[2]-1))
                    dups = nonunique(first.(cs_mets))
                end

                # Collect list with weights
                for c in 1:length(cs_mets)
                    push!(#
                        arcs,
                        KEGGArc(#
                            cs_mets[c][1][1],
                            cs_mets[c][2][1],
                            D[i].w[j],
                            [cs_rxns[c]]
                        )
                    )
                end
            end
        end
    end

    # Aggregate pairs involving transitions between the same metabolites
    hh(x) = [x.head; x.tail]
    jj(x) = x.weight
    ll(x) = x.subsystems
    ss(x) = sum(x[.!isnan.(x)])
    dups = nonunique(hh.(arcs))
    while !isempty(dups)
        ids = findall(==(popfirst!(dups)), hh.(arcs))
        topush = KEGGArc(#
            arcs[ids[1]].head,
            arcs[ids[1]].tail,
            ss(jj.(arcs[ids])),
            reduce(vcat, ll.(arcs[ids]))
        )
        deleteat!(arcs, ids)
        push!(arcs, topush)
    end

    # Remove source flux link
    kk(x) = x.head
    deleteat!(arcs, findall(==(0), kk.(arcs)))

    # Directed, weighted, but maybe cyclic graph
    return arcs, root
end

# For HepG2 linked lists
function kegg_linked_list(#
    D::Vector{CHMCAtomicSummary},
    dc::Vector{Dict{Tuple{Int64, Int64}, Int64}},
    subsystems::Vector{String},
    ur::UnitRange{Int64}
)
    # Root/source metabolite index
    root = first(D[1].dmc[1]) # MC state 1 is always the root

    # Collect atomic CHMC EFMs involving source metabolite as linked list
    arcs = Vector{KEGGArc}()

    for i in 1:length(D)
        # Identify MC state associated with the pseudo-metabolite state (0, 0)
        pm_idx = findfirst(==((0, 0)), collect(values(D[i].dmc)))
        pm_mc = collect(keys(D[i].dmc))[pm_idx]

        # Extract EFM sequences/weights within the given bin range
        sort_idx = sortperm(D[i].w, rev = true)
        range_idx = sort_idx[ur]
        ws = D[i].w[range_idx]
        es = D[i].e[range_idx]

        # Loop over range
        for j in 1:length(es)
            # Pair up metabolite index sequence
            states = first.([D[i].dmc[e] for e in first(es[j])])
            cs_mets = [pair for pair in IterTools.partition(states, 2, 1)]

            # Identify reaction sequence
            cs_mc = [pair for pair in IterTools.partition(first(es[j]), 2, 1)]
            cs_rxns = [dc[i][c] for c in cs_mc]

            # Identify and remove sequences involving the same metabolites
            dups = nonunique(first.(cs_mets))
            while !isempty(dups)
                ids = findall(==(dups[1]), first.(cs_mets))
                deleteat!(cs_mets, UnitRange((ids[1]), ids[2]-1))
                deleteat!(cs_rxns, UnitRange((ids[1]), ids[2]-1))
                dups = nonunique(first.(cs_mets))
            end

            # Collect list with weights
            for c in 1:length(cs_mets)
                push!(#
                    arcs,
                    KEGGArc(#
                        cs_mets[c][1][1],
                        cs_mets[c][2][1],
                        ws[j],
                        [cs_rxns[c]]
                    )
                )
            end
        end
    end

    # Aggregate pairs involving transitions between the same metabolites
    hh(x) = [x.head; x.tail]
    jj(x) = x.weight
    ll(x) = x.subsystems
    ss(x) = sum(x[.!isnan.(x)])
    dups = nonunique(hh.(arcs))
    while !isempty(dups)
        ids = findall(==(popfirst!(dups)), hh.(arcs))
        topush = KEGGArc(#
            arcs[ids[1]].head,
            arcs[ids[1]].tail,
            ss(jj.(arcs[ids])),
            reduce(vcat, ll.(arcs[ids]))
        )
        deleteat!(arcs, ids)
        push!(arcs, topush)
    end

    # Remove source flux link
    kk(x) = x.head
    deleteat!(arcs, findall(==(0), kk.(arcs)))

    # Directed, weighted, but maybe cyclic graph
    return arcs, root
end

# For HepG2 linked lists
function compress_matrix(flux_matrix::Matrix, kegg_matrix::Matrix; root::Int16=1, keep::Vector{Int64}=[1])
    # Identify states with a single incoming and outgoing transition
    n = size(flux_matrix,1)
    middle_states = Int64[]
    for state in 1:n
        if length(findall(x -> x > 0, flux_matrix[state,:])) == 1 &&
            length(findall(x -> x > 0, flux_matrix[:,state])) == 1 &&
            findfirst(x -> x > 0, flux_matrix[state,:]) != root
            push!(middle_states, state)
        end
    end

    # Store middle states between each (from, to) pair
    middle_states_removed = Dict{Tuple{Int16, Int16}, Vector{Int16}}()

    # Update weights
    for middle in middle_states
        from = findfirst(x -> x > 0, flux_matrix[:,middle])
        to = findfirst(x -> x > 0, flux_matrix[middle,:])
        if from != to
            flux_matrix[from, to] += flux_matrix[from, middle]
            flux_matrix[from, middle] = 0
            flux_matrix[middle, to] = 0

            kegg_matrix[from, to] = join([kegg_matrix[from, to]; kegg_matrix[from, middle]], " ")
            kegg_matrix[from, middle] = ""
            kegg_matrix[middle, to] = ""

            # Record middle states for each (from, to) pair
            if !haskey(middle_states_removed, (from, to))
                middle_states_removed[(from, middle)] = [to]
            else
                push!(middle_states_removed[(from, middle)], to)
            end
        end
    end

    # Assemble removed middle states into sequences
    function ismember(x::Pair{Tuple{Int16,Int16}, Vector{Int16}}, y::Vector{Int16})
        if x[1][1] == y[1] && x[1][2] == y[end]
            return true
        else
            return false
        end
    end
    assembly = Vector{Vector{Int16}}()
    while !isempty(middle_states_removed)
        x = pop!(middle_states_removed)
        i = findfirst(a -> ismember(x, a), assembly)
        if isnothing(i)
            push!(assembly, [x[1][1], x[1][2], x[2][1]])
        else
            push!(assembly[i], x[2][1])
        end
    end

    # Convert n+1 sink node indices back to 0
    for a in assembly
        a[findall(x -> x == n, a)] .= 0
    end

  return flux_matrix, kegg_matrix, assembly
end

# For HepG2 linked lists
function compress(arcs::Vector{KEGGArc}, root::Int16; keep::Vector{Int64}=[1])
    # Convert linked list to square atomic flux matrix
    src, dst, wts, keg = to_arrays(arcs)
    n = maximum([src; dst]) + 1
    dst[findall(==(0), dst)] .= n
    flux_matrix = zeros(n, n)
    kegg_matrix = fill("", n, n)
    for i in 1:length(src)
        flux_matrix[src[i], dst[i]] = wts[i]
        kegg_matrix[src[i], dst[i]] = join(keg[i], " ")
    end

    T, U, assembly = compress_matrix(flux_matrix, kegg_matrix; root = root, keep = keep)
    #T, U, assembly = compress_matrix(flux_matrix, kegg_matrix; root = Int16(1), keep = keep)

    # Convert flux matrix back to arcs
    cids = findall(!=(0), T)
    keg = [parse.(Int64, filter!(x -> .!isempty.(x), s)) for s in split.(strip.(U[cids]), " ")]
    arcs = KEGGArc.(first.(Tuple.(cids)), last.(Tuple.(cids)), T[cids], keg)
    return arcs, assembly
end

# For HepG2 linked lists
function sankey_atomic_chmc(#
    arcs::Vector{KEGGArc},
    root::Int16,
    mets::Vector{String},
    subsystems::Vector{String};
    cutoff::Float64   = 0.01,
    compression::Bool = false,
    font_size::Int64  = 20,
    atom_type::Symbol = :Carbon
)
    # Compress linear pathways
    if compression == true
        arcs, assembly = compress(arcs, root)
    end

    # Split arcs to source/destination/weights/KEGG index vectors
    src, dst, wts, keg = to_arrays(arcs)
    sink_idx = length(mets) + 1 # index for sink labels
    dst[dst .== 0] .= sink_idx # 0 corresponds to the sink pseudo-metabolite state

    # KEGG subsystem names from indices
    link_ids = unique.(keg)
    link_labels = [subsystems[i][1] for i in link_ids]
    for i in 1:length(link_labels)
        if isempty(link_labels[i])
            link_labels[i] = [""]
        end
    end

    # Remove arcs/weights/KEGG labels based on threshold
    idx = findall(x -> x >= cutoff, wts)
    src = src[idx]
    dst = dst[idx]
    wts = wts[idx]
    link_labels = link_labels[idx]

    # Must re-index arcs starting from 0 for PlotlyJS.sankey
    ks = unique([src; dst])
    d = Dict(zip(ks, 0:(length(ks)-1)))  # arc index to sankey index
    d′ = Dict(zip(0:(length(ks)-1), ks)) # sankey index to arc index

    # Set source/destination indices to sankey indices
    src = [d[i] for i in src]
    dst = [d[i] for i in dst]

    # Convert source/destination sankey indices to metabolite names
    names = Vector{String}(undef, length(unique([src; dst])))
    k = 0
    for i in sort(unique([src; dst]))
        k += 1
        if d′[i] == sink_idx
            names[k] = "Sink"
        else
            names[k] = mets[d′[i]]
        end
    end

    atom_type = string(atom_type)
    # Sankey plot
    p = PlotlyJS.plot(#
        PlotlyJS.sankey(#
            valueformat = ".3f",
            valuesuffix = "μmol gdw⁻¹h⁻¹",
            #arrangement = "perpendicular",
            node = attr(#
                pad = 15,
                thickness = 15,
                line = attr(color = "black", width = 0.5),
                label = names,
                color = repeat(#
                    [:blue, :green, :red, :purple, :yellow, :orange, :brown],
                    length(ks)
                )
            ),
            link = attr(
                source = src,
                target = dst,
                value  = abs.(wts),
                label  =  link_labels
            )
        ),
        Layout(#
            hovermode = "x",
            title = attr(#
                text = "$(atom_type) flux of source metabolite: $(mets[root])",
                x = 0.5,
                y = 1.0,
                xanchor = "center",
                yanchor = "top",
                font_size = font_size
            ),
            font = attr(size = 10, color = "black"),
            plot_bgcolor  = "white",
            paper_bgcolor = "white"
        )
    )
  return p
end

# For HepG2 linked lists
function to_arrays(X::Vector{KEGGArc})
    src = [x.head for x in X]
    dst = [x.tail for x in X]
    wgt = [x.weight for x in X]
    keg = [x.subsystems for x in X]
    return src, dst, wgt, keg
end

# For pairwise
function get_source_met_all_summary_dataframe(#
    D::Vector{Vector{@NamedTuple{min::Int64, max::Int64, ratio::Float64}}},
)
    h(x) = x.ratio
    df = DataFrame(#
        x = Real[], y = Real[]
    )
    for i in 1:length(D)
        for j in 1:length(D[i])
            push!(df, [i, h(D[i][j]),        ])
        end
    end
    return df
end

# For state space
function get_atomic_state_count_summary_dataframe(#
    F::Vector{Vector{Float64}},
    N::Vector{Vector{Int64}},
    jit::Real
)
    df = DataFrame(x = Real[], y = Real[], n = Real[], z = Real[])
    for i in eachindex(F)
        for j in 1:length(F[i])
            push!(#
                df,
                (i, F[i][j], N[i][j], jitter(i, jit))
            )
        end
    end
   
    return df
end

# For HepG2 linked lists
function export_top_X_efms_as_edge_list(#
    C::Vector{CHMCAtomicSummary},
    x::Int64,
    d::Vector{Dict{Tuple{Int64,Int64},Int64}},
    mets::Vector{String},
    rxns::Vector{String},
)

    df = DataFrame(#
        Source     = String[],  # source node
        Target     = String[],  # destination node
        Label      = String[],  # reaction name
        Weight     = Float64[], # atomic flux
        NumberEFMs = Int64[],   # number of EFMs involving that edge
    )
    
    # Loop over atom type index
    for i in 1:length(C)

        # EFM indices with greatest weights 
        top_ids = sortperm(C[i].w, rev = true)[1:x]

        for j in top_ids
            # Pair up metabolite index sequence
            states = first.([C[i].dmc[e] for e in first(C[i].e[j])])
            source_to_sink = false
            if states[1] == 0
                source_to_sink = true
                popfirst!(states)
                pop!(states)
            end
            cs_mets = [pair for pair in IterTools.partition(states, 2, 1)]

            # Identify reaction sequence
            cs_mc = [pair for pair in IterTools.partition(first(C[i].e[j]), 2, 1)]
            cs_rxns = [d[i][c] for c in cs_mc]
            if source_to_sink == true
                popfirst!(cs_rxns)
                pop!(cs_rxns)
            end

            if source_to_sink == true
                push!(#
                    df,
                    (#
                        Source     = "Source",
                        Target     = mets[cs_mets[1][1]],
                        Label      = "Source_reaction",
                        Weight     = C[i].w[j],
                        NumberEFMs = 1
                    )
                )
                push!(#
                    df,
                    (#
                        Source     = mets[cs_mets[end][2]],
                        Target     = "Sink",
                        Label      = "Sink_reaction",
                        Weight     = C[i].w[j],
                        NumberEFMs = 1
                    )
                )
            end

            for k in 1:length(cs_rxns)
                push!(#
                    df,
                    (#
                        Source     = mets[cs_mets[k][1]],
                        Target     = mets[cs_mets[k][2]],
                        Label      = rxns[cs_rxns[k]],
                        Weight     = C[i].w[j],
                        NumberEFMs = 1
                    )
                )
            end
        end
    end

    gdf = DataFrames.groupby(df, [:Source, :Target])
    gdfs = combine(gdf, [:Weight, :NumberEFMs] .=> sum; renamecols=false)

    return gdfs
end

# For HepG2 linked lists
function export_top_X_efms_as_edge_list(#
    C::Vector{CHMCAtomicSummary},
    x::UnitRange{Int64},
    d::Vector{Dict{Tuple{Int64,Int64},Int64}},
    mets::Vector{String},
    rxns::Vector{String},
)

    df = DataFrame(#
        Source     = String[],  # source node
        Target     = String[],  # destination node
        Label      = String[],  # reaction name
        Weight     = Float64[], # atomic flux
        NumberEFMs = Int64[],   # number of EFMs involving that edge
    )
    
    # Loop over atom type index
    for i in 1:length(C)

        # EFM indices with greatest weights 
        top_ids = sortperm(C[i].w, rev = true)[x]

        for j in top_ids
            # Pair up metabolite index sequence
            states = first.([C[i].dmc[e] for e in first(C[i].e[j])])
            source_to_sink = false
            if states[1] == 0
                source_to_sink = true
                popfirst!(states)
                pop!(states)
            end
            cs_mets = [pair for pair in IterTools.partition(states, 2, 1)]

            # Identify reaction sequence
            cs_mc = [pair for pair in IterTools.partition(first(C[i].e[j]), 2, 1)]
            cs_rxns = [d[i][c] for c in cs_mc]
            if source_to_sink == true
                popfirst!(cs_rxns)
                pop!(cs_rxns)
            end

            if source_to_sink == true
                push!(#
                    df,
                    (#
                        Source     = "Source",
                        Target     = mets[cs_mets[1][1]],
                        Label      = "Source_reaction",
                        Weight     = C[i].w[j],
                        NumberEFMs = 1
                    )
                )
                push!(#
                    df,
                    (#
                        Source     = mets[cs_mets[end][2]],
                        Target     = "Sink",
                        Label      = "Sink_reaction",
                        Weight     = C[i].w[j],
                        NumberEFMs = 1
                    )
                )
            end

            for k in 1:length(cs_rxns)
                push!(#
                    df,
                    (#
                        Source     = mets[cs_mets[k][1]],
                        Target     = mets[cs_mets[k][2]],
                        Label      = rxns[cs_rxns[k]],
                        Weight     = C[i].w[j],
                        NumberEFMs = 1
                    )
                )
            end
        end
    end

    gdf = DataFrames.groupby(df, [:Source, :Target])
    gdfs = combine(gdf, [:Weight, :NumberEFMs] .=> sum; renamecols=false)

    return gdfs
end

# For HepG2 linked lists
function kegg_subsystems_from_hmr(path::String, rxns::Vector{String})
    xf = XLSX.readdata(path, "RXNS", "A:K")
    ls = Vector{String}(undef, length(rxns))
    for i in 1:length(rxns)
        ii = findfirst(==(rxns[i]), xf[:,2])
        if isnothing(ii)
            if !isnothing(match(r"sink|Sink", rxns[i]))
                ls[i] = "Sink"
            elseif !isnothing(match(r"source|Source", rxns[i]))
                ls[i] = "Source"
            else
                ls[i] = "Sink" # assumed to be sink pseudoreactions
            end
        else
            ls[i] = xf[ii, 11]
        end
    end
    return ls
end

# For HepG2 linked lists
function nonunique(x::AbstractArray{T}) where T
    xs = sort(x)
    dups = T[]
    for i in 2:length(xs)
        if (isequal(xs[i], xs[i-1]) && (length(dups) == 0 || !isequal(dups[end], xs[i])))
            push!(dups, xs[i])
        end
    end
    return dups
end

# For benchmarking
function lower_upper_whisker(x::Vector{<:Number})
    lq = percentile(x, 25)
    uq = percentile(x, 75)
    x_iqr = iqr(x)
    lw = minimum(x[findall(x .> lq - 1.5 * x_iqr)])
    uw = maximum(x[findall(x .< uq + 1.5 * x_iqr)])
    return lw, uw
end

# For HepG2 linked lists
function precompute_linked_list_reactions(D::Vector{CHMCAtomicSummary})
    # Precompute dictionary of MC states i -> i+1 and reaction index
    d = Vector{Dict{Tuple{Int64, Int64}, Int64}}(undef, length(D))
    for i in 1:length(D)
        ks = last.(collect(keys(D[i].dchmc)))
        vs = first.(collect(values(D[i].dchmc)))
        tmp = Dict{Tuple{Int64, Int64}, Int64}()
        for k in D[i].R
            push!(tmp, (ks[findfirst(==(k.i), vs)], ks[findfirst(==(k.j), vs)]) => k.k)
        end
        d[i] = tmp
    end
    return d
end

# For HepG2 linked lists
function efms_to_r_node_dataframe(x::Vector{Vector{Int16}}, mets::Vector{String})
    i = unique(reduce(vcat, x))

    g(x) = x[1] != x[end] ? true : false

    if any(g.(x) .== true)
        sink_node = true
        sink_ids = last.(x[g.(x) .== true])
    else
        sink_node = false
        sink_ids = nothing
    end

    ID = ["ID" * string(ii) for ii in i]
    
    df = DataFrame(#
        ID = ID,
        label = mets[i],
        label_pos = repeat(["left"], length(i)),
        x = repeat(["0"], length(i)),
        y = repeat(["0"], length(i)),
        dir = repeat(["right"], length(i))
    )

    if sink_node == true
        push!(#
            df,
            (#
                "0",
                "Sink",
                "left",
                "0",
                "0",
                "right"
            )
        )
    end

    return df
end

# For HepG2 linked lists
function efms_to_r_flow_dataframe(x::Vector{Vector{Int16}}, w::Vector{Float64})
    # Add sink node to source-to-sink EFMs 
    x = deepcopy(x)
    g(x) = x[1] != x[end] ? true : false
    src_to_sink_ids = g.(x)
    if any(src_to_sink_ids)
        push!.(x[src_to_sink_ids], 0)
    end

    df = DataFrame(#
        from = String[],
        to = String[],
        substance = String[],
        quantity = Float64[]
    )

    for i in eachindex(x)
        for j in 1:(length(x[i])-1)
            push!(#
                df,
                (#
                    "ID" * string(x[i][j]),
                    "ID" * string(x[i][j+1]),
                    "EFM $i",
                    w[i],
                )
            )
        end
    end

    return df
end

