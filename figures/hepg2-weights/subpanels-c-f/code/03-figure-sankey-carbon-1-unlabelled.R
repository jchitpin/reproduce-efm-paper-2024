library("grid")
library("PantaRhei")

dir_inputs <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-1/01-for-R/"
dir_output <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-1/02-R-output-figures/sankey-unlabelled.pdf"

# Sankey diagram for glutamine bin 1 (5 carbon EFMs)
nodes1 <- read.csv(paste0(dir_inputs, "/nodes.csv"))
nodes1$x <- as.character(nodes1$x)
nodes1$y <- as.character(nodes1$y)
flows1 <- read.csv(paste0(dir_inputs, "/flows.csv"))

# Manually tweak nodes to prevent them from overlapping
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
nodes1$label <- capFirst(nodes1$label)

nodes1[3, ]$label_pos <- "above" # ID320. OAA[m]
nodes1[4, ]$x <- "ID320+2"       # ID137. Aspartate[m]
nodes1[4, ]$y <- "ID320"
nodes1[4, ]$label_pos <- "above"
nodes1[7, ]$x <- "ID320+2"       # ID154. Citrate[m]
nodes1[7, ]$y <- "ID320-1.0"
nodes1[7, ]$label_pos <- "above"
nodes1[5, ]$x <- "ID137+2"       # ID136. Aspartate[c]
nodes1[5, ]$y <- "ID137"
nodes1[5, ]$label_pos <- "above"
nodes1[8, ]$x <- "ID154+2"       # ID153. Citrate[c]
nodes1[8, ]$y <- "ID154"
nodes1[8, ]$label_pos <- "above"
nodes1[6, ]$x <- "ID320+6"      # ID319. OAA[c]
nodes1[6, ]$y <- "ID320+0.15"
nodes1[6, ]$label_pos <- "above"
nodes1[19,]$x <- "ID153+2"     # ID262. Isocitrate[c]
nodes1[19,]$y <- "ID153"
nodes1[19,]$label_pos <- "above"
nodes1[20,]$x <- "ID262+2"     # ID156. CO2[c]
nodes1[20,]$y <- "ID262"
nodes1[20,]$label_pos <- "above"
nodes1[21,]$x <- "ID156+2"     # ID158. CO2[s]
nodes1[21,]$y <- "ID156"
nodes1[21,]$label_pos <- "above"
nodes1[22,]$x <- "ID158+2"     # ID0. Sink
nodes1[22,]$y <- "ID158"
nodes1[22,]$label_pos <- "above"
nodes1[22,]$ID <- "ID0"
nodes1[13,]$x <- "ID320"         # ID219. Glutamine[s]
nodes1[13,]$y <- "ID320-3.0"
nodes1[13,]$label_pos <- "above"
nodes1[14,]$x <- "ID219+1.5"       # ID217. Glutamine[c]
nodes1[14,]$y <- "ID219"
nodes1[14,]$label_pos <- "above"
nodes1[15,]$x <- "ID217+1.5"       # ID218. Glutamine[m]
nodes1[15,]$y <- "ID217"
nodes1[15,]$label_pos <- "above"
nodes1[10,]$x <- "ID218+1.5"       # ID215. Glutamate[m]
nodes1[10,]$y <- "ID218-2.05"
nodes1[10,]$label_pos <- "above"
nodes1[11,]$x <- "ID215+2"       # ID125. AKG[m]
nodes1[11,]$y <- "ID215+0.23"
nodes1[11,]$label_pos <- "above"
nodes1[12,]$x <- "ID125+2"       # ID124. AKG[c]
nodes1[12,]$y <- "ID125"
nodes1[12,]$label_pos <- "above"
nodes1[9, ]$x <- "ID217-1"         # ID214. Glutamate[c]
nodes1[9, ]$y <- "ID217-2.16"
nodes1[9, ]$label_pos <- "above"
nodes1[2, ]$x <- "ID281-2.5"       # ID282. Malate[m]
nodes1[2, ]$y <- "ID281-5.50"
nodes1[2, ]$label_pos <- "above"
nodes1[16,]$x <- "ID125+1.5"     # ID390. Succinyl-CoA[m]
nodes1[16,]$y <- "ID125+0.5"
nodes1[16,]$label_pos <- "above"
nodes1[17,]$x <- "ID390+1.5"       # ID389. Succinate[m]
nodes1[17,]$y <- "ID390+0"
nodes1[17,]$label_pos <- "above"
nodes1[18,]$x <- "ID389+1.5"       # ID200. Fumarate[m]
nodes1[18,]$y <- "ID389+0"
nodes1[18,]$label_pos <- "above"

nodes1[1, ]$x <- "ID319+10"      # ID281. Malate[c]
nodes1[1, ]$y <- "ID319"
nodes1[1, ]$label_pos <- "below"

# Bottom malate-aspartate shuttle loop
nodes1[nrow(nodes1) + 1,] = c(".R1.AKG.GLU", "", "none", "ID124+1.0", "ID124-1.25", "down")
nodes1[nrow(nodes1) + 1,] = c(".R2.AKG.GLU", "", "none", "ID215", "ID215-2", "left")
nodes1[nrow(nodes1) + 1,] = c(".R3.AKG.GLU", "", "none", "ID217-2.0", "ID217-2.95", "up")

# Right malate-malate loop
nodes1[nrow(nodes1) + 1,] = c(".R1.MAL.MAL", "", "none", "ID281+1.75", "ID281-1.5", "down")
nodes1[nrow(nodes1) + 1,] = c(".R2.MAL.MAL", "", "none", "ID281", "ID281-10", "left")
nodes1[nrow(nodes1) + 1,] = c(".R3.MAL.MAL", "", "none", "ID282-2.5", "ID282-2.0", "up")

# Right malate-OAA loop
nodes1[nrow(nodes1) + 1,] = c(".R1.MAL.OAA", "", "none", "ID282+2", "ID282-1.5", "down")
nodes1[nrow(nodes1) + 1,] = c(".R2.MAL.OAA", "", "none", "ID282", "ID281-13", "left")
nodes1[nrow(nodes1) + 1,] = c(".R3.MAL.OAA", "", "none", "ID320-2", "ID320-2.0", "up")


flows1[nrow(flows1) + 1,] = list("ID124",      "R1.AKG.GLU", "EFM 3", 47.15961275393425)
flows1[nrow(flows1) + 1,] = list("R1.AKG.GLU", "R2.AKG.GLU", "EFM 3", 47.15961275393425)
flows1[nrow(flows1) + 1,] = list("R2.AKG.GLU", "R3.AKG.GLU", "EFM 3", 47.15961275393425)
flows1[nrow(flows1) + 1,] = list("R3.AKG.GLU", "ID214",      "EFM 3", 47.15961275393425)
flows1 <- flows1[-c(16, 1, 2, 7, 8, 17, 27), ] # 1/2 is EFM1; 7/8 is EFM2; 17/18 is EFM 4

flows1[nrow(flows1) + 1,] = list("ID281",      "R1.MAL.MAL", "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("R1.MAL.MAL", "R2.MAL.MAL", "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("R2.MAL.MAL", "R3.MAL.MAL", "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("R3.MAL.MAL", "ID282",      "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("ID281",      "R1.MAL.MAL", "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("R1.MAL.MAL", "R2.MAL.MAL", "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("R2.MAL.MAL", "R3.MAL.MAL", "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("R3.MAL.MAL", "ID282",      "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("ID281",      "R1.MAL.MAL", "EFM 4", 23.51639)
flows1[nrow(flows1) + 1,] = list("R1.MAL.MAL", "R2.MAL.MAL", "EFM 4", 23.51639)
flows1[nrow(flows1) + 1,] = list("R2.MAL.MAL", "R3.MAL.MAL", "EFM 4", 23.51639)
flows1[nrow(flows1) + 1,] = list("R3.MAL.MAL", "ID282",      "EFM 4", 23.51639)

flows1[nrow(flows1) + 1,] = list("ID282",      "R1.MAL.OAA", "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("R1.MAL.OAA", "R2.MAL.OAA", "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("R2.MAL.OAA", "R3.MAL.OAA", "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("R3.MAL.OAA", "ID320",      "EFM 1", 61.40487)
flows1[nrow(flows1) + 1,] = list("ID282",      "R1.MAL.OAA", "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("R1.MAL.OAA", "R2.MAL.OAA", "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("R2.MAL.OAA", "R3.MAL.OAA", "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("R3.MAL.OAA", "ID320",      "EFM 2", 53.03662)
flows1[nrow(flows1) + 1,] = list("ID282",      "R1.MAL.OAA", "EFM 5", 18.02056)
flows1[nrow(flows1) + 1,] = list("R1.MAL.OAA", "R2.MAL.OAA", "EFM 5", 18.02056)
flows1[nrow(flows1) + 1,] = list("R2.MAL.OAA", "R3.MAL.OAA", "EFM 5", 18.02056)
flows1[nrow(flows1) + 1,] = list("R3.MAL.OAA", "ID320",      "EFM 5", 18.02056)

flows1$substance <- paste0(flows1$substance, ": w = ", round(flows1$quantity, 1))

# PALETTE COLOURS
palette1 = data.frame(#
  substance = unique(flows1$substance),
  color = c("#017ec3", "#f31e26", "#019e5e", "#fba61d", "#916237")
)

# REMAINING COLOURS
text_title_col <- "#000000"
text_node_col  <- "#000000"
text_flow_col  <- "#000000"
fill_node_col  <- "#404040"

# FONT SIZE
font_master <- 12

# TITLE
title1 <- ""
attr(title1, "gp") <- grid::gpar(#
  fontsize = 16,
  fontface = "bold",
  col = text_title_col
)

# NODE STYLE
ns <- list(#
  type     = "arrow",
  gp       = gpar(fill = fill_node_col, col = "white", lwd = 2),
  length   = 0.7,
  label_gp = gpar(fontsize = font_master, col = text_node_col),
  mag_pos  = "label",
  mag_fmt  = "%.0f mu",
  mag_gp   = gpar(fontsize = 0, fontface = "bold", col = text_flow_col)
)
nodes1$label <- ""

pdf(dir_output, width = 7, height = 7) # Set up PDF device
sankey(#
  nodes1,
  flows1,
  palette1,
  max_width   = 0.1,
  rmin        = 0.5,
  node_style  = ns,
  page_margin = c(0, 0.1, 0, 0.1), # left bottom right top
  #legend      = gpar(col = text_node_col, ncols = 4, fontsize = font_master),
  title       = title1,
  grill       = FALSE, # TRUE for debugging
  copyright   = ""
)
dev.off()


