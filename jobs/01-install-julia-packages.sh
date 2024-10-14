#!/usr/bin/env bash

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

julia --project=$BASE_PATH src/install-julia-packages.jl

printf "\n"
printf "JULIA PACKAGE INSTALLATION SCRIPT COMPLETED.\n\n"


