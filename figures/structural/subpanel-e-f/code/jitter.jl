## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the aggregated carbon and nitrogen CHMC output files.
# (2) For carbons within each source metabolite of each GEM, calculate the
#     max over min ratio of:
#     * Number of carbon EFMs.
#     * CHMC state space.
#     * CHMC states (metabolite/atom position)
#     * Metabolite (ignoring atom position)
# (3) Partition source metabolites into categories.
# (4) Export the statistics in three long long dataframe
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
datasets = ["e_coli_core", "iAB_RBC_283", "iIT341", "iSB619", "hepg2"]

# Import and export directories
load1 = "data/"
load2 = "/gems/"
save1 = "figures/structural/subpanel-e-f/"

# Import filenames
im_carbon    = "/gems-aggregated/gems-carbon.jld2"
im_nitrogen  = "/gems-aggregated/gems-nitrogen.jld2"

# Export filenames
ex_carbon_number_efms        = "figure-jitter-max-over-min-log10-number-efms-carbon.dat"
ex_carbon_state_space        = "figure-jitter-max-over-min-log10-state-space-carbon.dat"
ex_carbon_unique_states      = "figure-jitter-max-over-min-log10-unique-states-carbon.dat"
ex_carbon_unique_metabolites = "figure-jitter-max-over-min-log10-unique-metabolites-carbon.dat"
ex_nitrogen_number_efms        = "figure-jitter-max-over-min-log10-number-efms-nitrogen.dat"
ex_nitrogen_state_space        = "figure-jitter-max-over-min-log10-state-space-nitrogen.dat"
ex_nitrogen_unique_states      = "figure-jitter-max-over-min-log10-unique-states-nitrogen.dat"
ex_nitrogen_unique_metabolites = "figure-jitter-max-over-min-log10-unique-metabolites-nitrogen.dat"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))

# Setting seed for datapoint jittering!
Random.seed!(0)
# ------------------------------------------------------------------------------

## LOAD AGGREGATED DATASETS ----------------------------------------------------
@info("Loading aggregated datasets.")
C = load_aggregated_chmcs(load1 * im_carbon, datasets)
N = load_aggregated_chmcs(load1 * im_nitrogen, datasets)
# ------------------------------------------------------------------------------

## COMPUTE SUMMARY STATISTICS FOR EACH GEM -------------------------------------
# Carbon
C1 = get_number_of_atomic_efms(C)
C2 = get_atomic_chmc_state_space(C)
C3 = get_number_of_unique_metabolite_atom_states(C)
C4 = get_number_of_unique_metabolites(C)

D1 = get_max_over_min_ratio_within_source(C1)
D2 = get_max_over_min_ratio_within_source(C2)
D3 = get_max_over_min_ratio_within_source(C3)
D4 = get_max_over_min_ratio_within_source(C4)

# Nitrogen
N1 = get_number_of_atomic_efms(N)
N2 = get_atomic_chmc_state_space(N)
N3 = get_number_of_unique_metabolite_atom_states(N)
N4 = get_number_of_unique_metabolites(N)

O1 = get_max_over_min_ratio_within_source(N1)
O2 = get_max_over_min_ratio_within_source(N2)
O3 = get_max_over_min_ratio_within_source(N3)
O4 = get_max_over_min_ratio_within_source(N4)
# ------------------------------------------------------------------------------

## ASSIGN ROW INDICES CORRESPONDING TO UNIQUE METABOLITES ACROSS ALL GEMS ------
met_names = load_gem_metabolite_names(load1 * load2, datasets)

# List of source metabolites in each GEM
source_ids = load_chmc_source_metabolite_indices(load1 * load2, datasets)

# Map hepg2 names to bigg nomenclature in the other datasets
idx_hepg2 = findfirst(==("hepg2"), datasets)
met_names[idx_hepg2] = hepg2_names_to_bigg(met_names[idx_hepg2])

# Source metabolites
source_mets = [met_names[i][source_ids[i]] for i in eachindex(source_ids)]
met_remove = "3-beta-D-glucuronosyl-3-beta-D-galactosyl-4-beta-D-galactosyl-O-beta-D-xylosylprotein[g]"
idx = findfirst(==(met_remove), source_mets[idx_hepg2])
popat!(source_mets[idx_hepg2], idx)
popat!(source_ids[idx_hepg2], idx)
# ------------------------------------------------------------------------------

## RESHAPE DATA INTO DATAFRAMES ------------------------------------------------
# Set jitter (variance of a normal distribution)
jit = 0.1

# Assign labels based on source metabolite super class defined by HMDB
# x = dataset index, y = log10 statistic, z = jittered dataset index
dfc1 = get_source_met_summary_dataframe(D1, source_mets, jit)
dfc2 = get_source_met_summary_dataframe(D2, source_mets, jit)
dfc3 = get_source_met_summary_dataframe(D3, source_mets, jit)
dfc4 = get_source_met_summary_dataframe(D4, source_mets, jit)

dfn1 = get_source_met_summary_dataframe(O1, source_mets, jit)
dfn2 = get_source_met_summary_dataframe(O2, source_mets, jit)
dfn3 = get_source_met_summary_dataframe(O3, source_mets, jit)
dfn4 = get_source_met_summary_dataframe(O4, source_mets, jit)

# Log10 transform statistics (NaNs will be omitted in pgfplots)
dfc1.y = log10.(dfc1.y)
dfc2.y = log10.(dfc2.y)
dfc3.y = log10.(dfc3.y)
dfc4.y = log10.(dfc4.y)

dfn1.y = log10.(dfn1.y)
dfn2.y = log10.(dfn2.y)
dfn3.y = log10.(dfn3.y)
dfn4.y = log10.(dfn4.y)
# ------------------------------------------------------------------------------

## EXPORT SUPER CLASS MATRICES TO DATA FILES -----------------------------------
CSV.write(save1 * ex_carbon_number_efms,        dfc1, header = true, delim = " ")
CSV.write(save1 * ex_carbon_state_space,        dfc2, header = true, delim = " ")
CSV.write(save1 * ex_carbon_unique_states,      dfc3, header = true, delim = " ")
CSV.write(save1 * ex_carbon_unique_metabolites, dfc4, header = true, delim = " ")

CSV.write(save1 * ex_nitrogen_number_efms,        dfn1, header = true, delim = " ")
CSV.write(save1 * ex_nitrogen_state_space,        dfn2, header = true, delim = " ")
CSV.write(save1 * ex_nitrogen_unique_states,      dfn3, header = true, delim = " ")
CSV.write(save1 * ex_nitrogen_unique_metabolites, dfn4, header = true, delim = " ")
# ------------------------------------------------------------------------------

## MISCELLANEOUS STATISTICS ----------------------------------------------------
# Median ratio of number of carbon/nitrogen AEFMs without ratios of one
f(x) = x[findall(!isinf, x)]
g(x) = x[x .!= 0]

# E coli core
median(g(f( dfc1[(dfc1.x .== 1), :y] ))) # 0.58 ratio on a log10 scale

# iAB RBC 283
median(g(f( dfc1[(dfc1.x .== 2), :y] ))) # 0.31 ratio on a log10 scale

# iIT341
median(g(f( dfc1[(dfc1.x .== 3), :y] ))) # 0.59 ratio on a log10 scale

# iSB619
median(g(f( dfc1[(dfc1.x .== 4), :y] ))) # 0.66 ratio on a log10 scale

# HepG2
median(g(f( dfc1[(dfc1.x .== 5), :y] ))) # 1.62 ratio on a log10 scale

# E coli core
#median(g(f( dfn1[(dfn1.x .== 1), :y] ))) # N/A since all ratios are 1/0

# iAB RBC 283
median(g(f( dfn1[(dfn1.x .== 2), :y] ))) # 0.31 ratio on a log10 scale

# iIT341
median(g(f( dfn1[(dfn1.x .== 3), :y] ))) # 0.29 ratio on a log10 scale

# iSB619
median(g(f( dfn1[(dfn1.x .== 4), :y] ))) # 0.13 ratio on a log10 scale

# HepG2
median(g(f( dfn1[(dfn1.x .== 5), :y] ))) # 1.05 ratio on a log10 scale
# ------------------------------------------------------------------------------


## For manually assigning pins to jitter plot
regex = r"Phenylalanine_e" # change!
regex = r"Methionine_e" # change!
regex = r"Histidine_e" # change!
regex = r"Isoleucine_e" # change!
regex = r"Lysine_e" # change!
regex = r"Leucine_e" # change!
regex = r"Tryptophan_e" # change!
regex = r"cholesta" # change!
regex = r"Glycine_e" # change!
regex = r"2-Oxoglutarate_e" # change!
regex = r"Glycerophosphoserine_c" # change!
regex = r"N-Acetyl-D-muramoyl-L-alanine_c" # change!
regex = r"S-Adenosyl-L-homocysteine_c" # change!
regex = r"Orotidine_5'-phosphate_c" # change!
regex = r"2-Oxoglutarate_c" # change!

df = dfc1 # change! dfc1/dfn1 are the statistics for max/min AEFMs
df.name[findall(.!isnothing.(match.(regex, df.name)))]
df.y[findall(.!isnothing.(match.(regex, df.name)))]
df.z[findall(.!isnothing.(match.(regex, df.name)))]

## For manually assigning pins to jitter plot
sort(dfn1, [:y])
regex = r"Histidine_e" # change!
regex = r"Arginine_e" # change!
regex = r"N-Acetyl-D-muramoyl-L-alanine_c" # change!
regex = r"N-Succinyl-LL-2,6-diaminoheptanedioate_c" # change!
regex = r"S-Adenosyl-L-homocysteine_c" # change!

df = dfn1 # change! dfc1/dfn1 are the statistics for max/min AEFMs
df.name[findall(.!isnothing.(match.(regex, df.name)))]
df.y[findall(.!isnothing.(match.(regex, df.name)))]
df.z[findall(.!isnothing.(match.(regex, df.name)))]

## E coli
sum(.!isinf.(dfc1.y)) # total number of source metabolite/carbons 443
sum(.!isinf.(dfn1.y)) # total number of source metabolites/nitrogen 202
sum(dfc1.y .> 0) # total number of source metabolites with non-zero carbon statistic 209
sum(dfn1.y .> 0) # total number of source metabolites with non-zero nitrogen statistic 70

# All GEMs
idx_carbon = 0 # equals 16 source metabolites with only 1 carbon
for i in eachindex(C1)
    for j in 1:length(C1[i])
        if length(C1[i][j]) == 1
            idx_carbon += 1
        end
    end
end
idx_carbon

idx_nitrogen = 0 # equals 108 source metabolites with only 1 nitrogen
for i in eachindex(N1)
    for j in 1:length(N1[i])
        if length(N1[i][j]) == 1
            idx_nitrogen += 1
        end
    end
end
idx_nitrogen
