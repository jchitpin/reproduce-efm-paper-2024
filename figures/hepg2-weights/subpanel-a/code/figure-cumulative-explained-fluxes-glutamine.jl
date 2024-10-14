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
ii = [14]
yc_threshold = 0.99
n_threshold = Vector{Vector{Int64}}(undef, length(ii))

for i in ii
    df = output_percentage_explained_source_flux_by_efms(Ch[i])

    # Cleaning up dataframe to remove "additional" points close to y=1
    output = Vector{DataFrame}(undef, length(Ch[i]))
    tmp = Int64[]
    for j in 1:length(Ch[i])
        output[j] = df[(df.atom .== j) .& (df.yc .<= yc_threshold), :]
        push!(tmp, size(df[(df.atom .== j) .& (df.yc .> yc_threshold), :])[1])
    end
    output = reduce(vcat, output)
    n_threshold[findfirst(==(i), ii)] = tmp # number of EFMs omitted from plots
    maximum(output.x) # for specifying xmax on plots

    CSV.write(ex_dataframes * "$i.dat", output, delim = ' ')
    CSV.write(ex_dataframes * "$i-full.dat", df, delim = ' ')
end

