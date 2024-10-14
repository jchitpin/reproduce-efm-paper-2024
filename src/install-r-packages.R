## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Installs all R packages necessary to run the scripts.
# ------------------------------------------------------------------------------

## INSTALLING PACKAGE DEPENDENCIES ---------------------------------------------
chooseCRANmirror(ind = 1)
list.of.packages <- c("grid", "PantaRhei")
new.packages <- list.of.packages[#
    !(list.of.packages %in% installed.packages()[,"Package"])
]
if(length(new.packages)) install.packages(new.packages)
# ------------------------------------------------------------------------------

