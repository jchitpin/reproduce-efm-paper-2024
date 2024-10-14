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
save1 = "figures/cumulative-mass-flow/subpanels/"

# Import filenames
im_carbon_hepg2   = "atomic-efm-weights-carbon"
im_nitrogen_hepg2 = "atomic-efm-weights-nitrogen"

# Export filenames (without -$i.csv)
ex_dataframes = save1 * "cumulative-source-met-"
ex_remaining = save1 * "remaining-99-aefms-source-met.csv"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD HEPG2 DATASETS ---------------------------------------------------------
Ch = load_hepg2_chmcs(load1 * im_carbon_hepg2)
Nh = load_hepg2_chmcs(load1 * im_nitrogen_hepg2)
# ------------------------------------------------------------------------------

# All amino acids plus glucose
ii = [6, 10, 13, 14, 15, 18, 20, 21, 22, 23, 29, 33, 35, 40]
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

CSV.write(ex_remaining, Tables.table(n_threshold), header=false)
# ------------------------------------------------------------------------------

glutamine_mets = unique([#
    unique(first.(collect(values(Ch[14][1].dmc))));
    unique(first.(collect(values(Ch[14][2].dmc))));
    unique(first.(collect(values(Ch[14][3].dmc))));
    unique(first.(collect(values(Ch[14][4].dmc))));
    unique(first.(collect(values(Ch[14][5].dmc))));
])
filter!(!=(0), glutamine_mets) # glutamine can remodel into 203 distinct metabolites (including glutamine itself and accounting for compartment separation)


