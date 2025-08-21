## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the carbon and nitrogen CHMC output files from the HepG2 dataset.
# (2) Computes cumulative, sorted, explained flux for individual atomic CHMCs.
# (3) Exports the data for plotting.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
dataset = "hepg2"

# Import and export directories
load1 = "data/gems/hepg2/"
save1 = "figures/hepg2-weights/subpanel-a/"

# Import filenames
im_carbon_hepg2   = "atomic-efm-weights-carbon"

# Export filenames (without -$i.csv)
ex_dataframes = save1 * "cumulative-source-met-"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD HEPG2 DATASETS ---------------------------------------------------------
Ch = load_hepg2_chmcs(load1 * im_carbon_hepg2)
# ------------------------------------------------------------------------------

# All amino acids plus glucose
i = 14
yc_threshold = 0.99
n_threshold = Vector{Vector{Int64}}(undef, 1)

output = Vector{DataFrame}(undef, length(Ch[i]))
df = output_percentage_explained_source_flux_by_efms(Ch[i])

# Cleaning up dataframe to remove "additional" points close to y=1
tmp = Int64[]
for j in 1:length(Ch[i])
    output[j] = df[(df.atom .== j) .& (df.yc .<= yc_threshold), :]
    push!(tmp, size(df[(df.atom .== j) .& (df.yc .> yc_threshold), :])[1])
end
output = reduce(vcat, output)

countmap(df.atom) # total number of EFMs
#   6,554
# 221,933
# 184,823
#  28,765
#   6,554

n_threshold[findfirst(==(i), [i])] = tmp # number of EFMs omitted from plots
#   6,316
# 221,420
# 184,318
#  28,005
#   4,552
maximum(output.x) # for specifying xmax on plots; value is 760

sum.(n_threshold) # 444,611 glutamine carbon AEFMs that explain less than 1%
sum.(n_threshold) / nrow(df) # 99.5% of pathways explain less than 1% of total source metabolite weights

CSV.write(ex_dataframes * "$i.dat", output, delim = ' ')
CSV.write(ex_dataframes * "$i-full.dat", df, delim = ' ')

# Explained carbon mass flow (accounting for AEFM length)
gln5_1 = output[(output.atom .== 1), :][5, :yc] # 54.05% explained
gln5_2 = output[(output.atom .== 2), :][5, :yc] # 45.46% explained
gln5_3 = output[(output.atom .== 3), :][5, :yc] # 45.64% explained
gln5_4 = output[(output.atom .== 4), :][5, :yc] # 35.17% explained
gln5_5 = output[(output.atom .== 5), :][5, :yc] # 39.15% explained

(gln5_1 + gln5_2 + gln5_3 + gln5_4 + gln5_5) / 5 # 43.89% explained

gln10_1 = output[(output.atom .== 1), :][10, :yc] # 64.89% explained
gln10_2 = output[(output.atom .== 2), :][10, :yc] # 58.86% explained
gln10_3 = output[(output.atom .== 3), :][10, :yc] # 59.10% explained
gln10_4 = output[(output.atom .== 4), :][10, :yc] # 47.00% explained
gln10_5 = output[(output.atom .== 5), :][10, :yc] # 51.95% explained

(gln10_1 + gln10_2 + gln10_3 + gln10_4 + gln10_5) / 5 # 56.36% explained


