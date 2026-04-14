## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Loads the aggregated carbon and nitrogen CHMC output files.
# (2) Computes summary statistics for each atomic EFM.
# (3) Exports these statistics into a LaTeX (.tex) table.
# ------------------------------------------------------------------------------

## DIRECTORIES AND FILE NAMES --------------------------------------------------
datasets = ["e_coli_core", "iAB_RBC_283", "iIT341", "iSB619", "hepg2"]
dataset_names = ["E. coli core", "iAB RBC 283", "iIT341", "iSB619", "HepG2"]

# Import and export directories
load1 = "data/"
load2 = "/gems/"

# Import filenames
im_carbon   = "/gems-aggregated/gems-carbon.jld2"
im_nitrogen = "/gems-aggregated/gems-nitrogen.jld2"
im_stoich   = "/processed/stoichiometry-matrix-processed.csv"
im_smiles   = "/processed/smiles-isomeric-processed.csv"
im_mets     = "/processed/metabolites-processed.csv"

# Export filenames
ex_df    = "tables/table-gem-chmc-summary/table-gem-chmc-summary.csv"
ex_table = "tables/table-gem-chmc-summary/table-gem-chmc-summary.tex"
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using MarkovWeightedEFMs, CSV, Tables, JLD2, DataFrames
dir_functions = pwd() * "/functions/"
include.(filter(contains(r".jl$"), readdir(dir_functions; join = true)))
# ------------------------------------------------------------------------------

## LOAD AGGREGATED DATASETS ----------------------------------------------------
@info("Loading aggregated datasets.")
C = load_aggregated_chmcs(load1 * im_carbon, datasets)
N = load_aggregated_chmcs(load1 * im_nitrogen, datasets)
# ------------------------------------------------------------------------------

## INITIALIZE AGGREGATED OUTPUT DATAFRAME --------------------------------------
df = DataFrame(#
    GEM                       = String[],
    Metabolites               = Int64[],
    Reactions                 = Int64[],
    Source_Metabolites        = Int64[],
    Source_Carbons            = Int64[],
    Source_Nitrogens          = Int64[],
    Sink_Metabolites          = Int64[],
    Sink_Carbons              = Int64[],
    Sink_Nitrogens            = Int64[],
    CHMC_State_Space_Carbon   = Int64[],
    CHMC_Transitions_Carbon   = Int64[],
    CHMC_State_Space_Nitrogen = Int64[],
    CHMC_Transitions_Nitrogen = Int64[],
    EFMs_Carbon               = Int64[],
    EFMs_Nitrogen             = Int64[],
)
# ------------------------------------------------------------------------------

## LOOP OVER DATASETS AND PUSH SUMMARY STATS TO DATAFRAME ----------------------
for i in eachindex(datasets)
    set = datasets[i]
    setname = dataset_names[i]
    @info "Completed: $set."

    # Load stoichiometry matrix
    S = Int16.(#
        CSV.read(#
            load1 * load2 * set * im_stoich,
            Tables.matrix,
            header = false
        )
    )

    # Load SMILES strings matching the stoichiometry rows
    smiles = vec(#
        CSV.read(#
            load1 * load2 * set * im_smiles,
            Tables.matrix,
            header = false
        )
    )

    # Load metabolites
    mets = vec(#
        CSV.read(#
            load1 * load2 * set * im_mets, Tables.matrix, header = false
        )
    )

    # Source/sink metabolite indices
    srcs = get_source_metabolites(Int16.(S))
    sncs = get_source_metabolites(Int16.(S * -1))

    # Remove 3-beta-D-glucuronsyl because it undergoes no internal reaction
    if set == "hepg2"
        met_remove = "3-beta-D-glucuronosyl-3-beta-D-galactosyl-4-beta-D-galactosyl-O-beta-D-xylosylprotein[g]"
        idx = findfirst(==(met_remove), mets[srcs])
        popat!(srcs, idx)
        idx = findfirst(==(met_remove), mets[sncs])
        popat!(sncs, idx)
    end

    # Total number of carbons/nitrogens in each metabolite
    amaxC = get_max_atoms(smiles, :C)
    amaxN = get_max_atoms(smiles, :N)

    sc = sum([length(k.dchmc) for j in C[i] for k in j])
    tc = sum([length(k.R)     for j in C[i] for k in j])
    sn = sum([length(k.dchmc) for j in N[i] for k in j])
    tn = sum([length(k.R)     for j in N[i] for k in j])
    ec = sum([length(k.e)     for j in C[i] for k in j])
    en = sum([length(k.e)     for j in N[i] for k in j])

    # Append info to dataframe
    push!(#
        df,
        (#
            setname,
            size(S, 1),
            size(S, 2),
            length(srcs),
            sum(amaxC[srcs]),
            sum(amaxN[srcs]),
            length(sncs),
            sum(amaxC[sncs]),
            sum(amaxN[sncs]),
            sc,
            tc,
            sn,
            tn,
            ec,
            en,
        )
    )
end
# ------------------------------------------------------------------------------

## PROGRAMMATICALLY CONSTRUCT TABLE --------------------------------------------
io = open(ex_table, "w");
write(io, "\\begin{tabular}{l rr rrrrr rrrrr}\n")
write(io, "\t\\toprule\n")
write(io, "\t& & & \\multicolumn{5}{c}{Carbon ACHMCs} & \\multicolumn{5}{c}{Nitrogen ACHMCs}\\\\\n")
write(io, "\t\\cmidrule(lr){4-8} \\cmidrule(lr){9-13}\n")
write(io, "\tGEM & Metabolites & Reactions & Sources & Sinks & States & Transitions & AEFMs & Sources & Sinks & States & Transitions & AEFMs\\\\\n")
write(io, "\t\\midrule\n")
for i in 1:nrow(df)
    nme  = replace(df[i, :GEM], "_" => "\\_")
    met  = string( df[i, :Metabolites])
    rxn  = string( df[i, :Reactions])
    srcC = string( df[i, :Source_Carbons])
    sncC = string( df[i, :Sink_Carbons])
    staC = string( df[i, :CHMC_State_Space_Carbon])
    traC = string( df[i, :CHMC_Transitions_Carbon])
    efmC = string( df[i, :EFMs_Carbon])
    srcN = string( df[i, :Source_Nitrogens])
    sncN = string( df[i, :Sink_Nitrogens])
    staN = string( df[i, :CHMC_State_Space_Nitrogen])
    traN = string( df[i, :CHMC_Transitions_Nitrogen])
    efmN = string( df[i, :EFMs_Nitrogen])
    line = join(#
        [nme, met, rxn, srcC, sncC, staC, traC, efmC, srcN, sncN, staN, traN, efmN],
        " & "
    )
    write(io, join(["\t", line, "\\\\\n"]))
end
write(io, "\t\\bottomrule\n")
write(io, "\\end{tabular}\n")
close(io)
# ------------------------------------------------------------------------------

## EXPORT SUMMARY DATAFRAME ----------------------------------------------------
CSV.write(ex_df, df)
# ------------------------------------------------------------------------------

