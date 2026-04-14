#!/usr/bin/env bash
THREAD=1

printf "This script enumerates the carbon and nitrogen EFMs for the following\n"
printf "genome-scale models (GEMs). For benchmarking purposes, Bash variable\n"
printf "THREAD is set to $THREAD. Set the thread count higher if you'd like\n"
printf "to reproduce the results faster.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)
ADD_PATH="/src/gems/"
CONSTRUCTION="01-dataset-construction"
ENUMERATION="02-efm-enumeration"
FNAME1="01-sbml-to-csv.jl"
FNAME2="02-construct-dataset.jl"

FNAME3="01-enumerate-atomic-efms.jl"
FNAME4="01-construct-dataset.jl"
# ------------------------------------------------------------------------------

## E_coli_core -----------------------------------------------------------------
SET="e_coli_core"
# Pre-process dataset
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME1
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME2
printf "Pre-processed model: e_coli_core.\n"

# Enumerate EFMs
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME3
printf "Enumerated EFMs for model: e_coli_core.\n\n"
# ------------------------------------------------------------------------------

## HepG2 -----------------------------------------------------------------------
SET="hepg2"
# Pre-process dataset
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME4
printf "Pre-processed model: hepg2.\n"

# Enumerate EFMs
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME3
printf "Enumerated EFMs for model: hepg2.\n\n"
# ------------------------------------------------------------------------------

## iAB_RBC_283 -----------------------------------------------------------------
SET="iAB_RBC_283"
# Pre-process dataset
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME1
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME2
printf "Pre-processed model: iAB_RBC_283.\n"

# Enumerate EFMs
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME3
printf "Enumerated EFMs for model: iAB_RBC_283.\n\n"
# ------------------------------------------------------------------------------

## iIT341 ----------------------------------------------------------------------
SET="iIT341"
# Pre-process dataset
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME1
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME2
printf "Pre-processed model: iIT341.\n"

# Enumerate EFMs
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME3
printf "Enumerated EFMs for model: iIT341.\n\n"
# ------------------------------------------------------------------------------

## iSB619 ----------------------------------------------------------------------
SET="iSB619"
# Pre-process dataset
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME1
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$CONSTRUCTION/$FNAME2
printf "Pre-processed model: iSB619.\n"

# Enumerate EFMs
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME3
printf "Enumerated EFMs for model: iSB619.\n\n"
# ------------------------------------------------------------------------------

printf "\n"
printf "ENUMERATION SCRIPT COMPLETED.\n\n"

