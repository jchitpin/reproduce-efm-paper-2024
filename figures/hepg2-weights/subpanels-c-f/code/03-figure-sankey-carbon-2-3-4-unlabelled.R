library("grid")
library("PantaRhei")

dir_inputs <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-2-3-4/01-for-R/"
dir_output <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-2-3-4/02-R-output-figures/sankey-unlabelled.pdf"

# Sankey diagram for glutamine bin 1 (5 carbon EFMs)
nodes2 <- read.csv(paste0(dir_inputs, "/nodes.csv"))
nodes2$x <- as.character(nodes2$x)
nodes2$y <- as.character(nodes2$y)
flows2 <- read.csv(paste0(dir_inputs, "/flows.csv"))

# Manually tweak nodes to prevent them from overlapping
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
nodes2$label <- capFirst(nodes2$label)

nodes2[3, ]$label_pos <- "above" # ID320. OAA[m]
nodes2[4, ]$x <- "ID320+4"       # ID137. Aspartate[m]
nodes2[4, ]$y <- "ID320"
nodes2[4, ]$label_pos <- "above"
nodes2[7, ]$x <- "ID320+4"       # ID154. Citrate[m]
nodes2[7, ]$y <- "ID320-2.0"
nodes2[7, ]$label_pos <- "above"
nodes2[5, ]$x <- "ID137+4"       # ID136. Aspartate[c]
nodes2[5, ]$y <- "ID137"
nodes2[5, ]$label_pos <- "above"
nodes2[8, ]$x <- "ID154+4"       # ID153. Citrate[c]
nodes2[8, ]$y <- "ID154"
nodes2[8, ]$label_pos <- "above"
nodes2[6, ]$x <- "ID320+12"       # ID319. OAA[c]
nodes2[6, ]$y <- "ID320+0.0"
nodes2[6, ]$label_pos <- "above"
nodes2[1, ]$x <- "ID319+4"       # ID281. Malate[c]
nodes2[1, ]$y <- "ID319"
nodes2[1, ]$label_pos <- "above"
nodes2[2, ]$x <- "ID281-3.0"    # ID282. Malate[m]
nodes2[2, ]$y <- "ID281-2.5"
nodes2[2, ]$label_pos <- "above"

nodes2[9, ]$x <- "ID320"         # ID124. AKG[c]
nodes2[9, ]$y <- "ID320-12.0"
nodes2[9, ]$label_pos <- "above"
nodes2[10,]$x <- "ID124+5.75"     # ID214. Glutamate[c]
nodes2[10,]$y <- "ID124"
nodes2[10,]$label_pos <- "above"
nodes2[11,]$x <- "ID214+5.75"     # ID215. Glutamate[m]
nodes2[11,]$y <- "ID214"
nodes2[11,]$label_pos <- "above"
nodes2[12,]$x <- "ID215+5.75"     # ID125. AKG[m]
nodes2[12,]$y <- "ID215"
nodes2[12,]$label_pos <- "above"

# Bottom malate-aspartate shuttle loop
nodes2[nrow(nodes2) + 1,] = c(".R1.AKG.AKG", "", "none", "ID125+1.5", "ID125-1.5", "down")
nodes2[nrow(nodes2) + 1,] = c(".R2.AKG.AKG", "", "none", "ID215", "ID215-2", "left")
nodes2[nrow(nodes2) + 1,] = c(".R3.AKG.AKG", "", "none", "ID124-1.5", "ID124-1.5", "up")

# Right malate-malate loop
nodes2[nrow(nodes2) + 1,] = c(".R1.MAL.MAL", "", "none", "ID281+2.5", "ID281-1.5", "down")
nodes2[nrow(nodes2) + 1,] = c(".R2.MAL.MAL", "", "none", "ID281", "ID281-7.5", "left")
nodes2[nrow(nodes2) + 1,] = c(".R3.MAL.MAL", "", "none", "ID282-2.5", "ID282-1.5", "up")

# Right malate-OAA loop
nodes2[nrow(nodes2) + 1,] = c(".R1.MAL.OAA", "", "none", "ID282+2.5", "ID282-1.5", "down")
nodes2[nrow(nodes2) + 1,] = c(".R2.MAL.OAA", "", "none", "ID282", "ID281-10", "left")
nodes2[nrow(nodes2) + 1,] = c(".R3.MAL.OAA", "", "none", "ID320-1.5", "ID320-1.75", "up")

flows2 <- flows2[-c(1, 2, 7, 8, 17, 18, 23, 16), ] # 1/2 is EFM1; 7/8 is EFM2; 17/18 is EFM4

flows2[nrow(flows2) + 1,] = list("ID125",      "R1.AKG.AKG", "EFM 3", 54.50813)
flows2[nrow(flows2) + 1,] = list("R1.AKG.AKG", "R2.AKG.AKG", "EFM 3", 54.50813)
flows2[nrow(flows2) + 1,] = list("R2.AKG.AKG", "R3.AKG.AKG", "EFM 3", 54.50813)
flows2[nrow(flows2) + 1,] = list("R3.AKG.AKG", "ID124",      "EFM 3", 54.50813)

flows2[nrow(flows2) + 1,] = list("ID281",      "R1.MAL.MAL", "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("R1.MAL.MAL", "R2.MAL.MAL", "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("R2.MAL.MAL", "R3.MAL.MAL", "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("R3.MAL.MAL", "ID282",      "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("ID281",      "R1.MAL.MAL", "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("R1.MAL.MAL", "R2.MAL.MAL", "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("R2.MAL.MAL", "R3.MAL.MAL", "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("R3.MAL.MAL", "ID282",      "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("ID281",      "R1.MAL.MAL", "EFM 4", 28.44492)
flows2[nrow(flows2) + 1,] = list("R1.MAL.MAL", "R2.MAL.MAL", "EFM 4", 28.44492)
flows2[nrow(flows2) + 1,] = list("R2.MAL.MAL", "R3.MAL.MAL", "EFM 4", 28.44492)
flows2[nrow(flows2) + 1,] = list("R3.MAL.MAL", "ID282",      "EFM 4", 28.44492)
flows2[nrow(flows2) + 1,] = list("ID281",      "R1.MAL.MAL", "EFM 5", 27.48700)
flows2[nrow(flows2) + 1,] = list("R1.MAL.MAL", "R2.MAL.MAL", "EFM 5", 27.48700)
flows2[nrow(flows2) + 1,] = list("R2.MAL.MAL", "R3.MAL.MAL", "EFM 5", 27.48700)
flows2[nrow(flows2) + 1,] = list("R3.MAL.MAL", "ID282",      "EFM 5", 27.48700)

flows2[nrow(flows2) + 1,] = list("ID282",      "R1.MAL.OAA", "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("R1.MAL.OAA", "R2.MAL.OAA", "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("R2.MAL.OAA", "R3.MAL.OAA", "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("R3.MAL.OAA", "ID320",      "EFM 1", 71.77272)
flows2[nrow(flows2) + 1,] = list("ID282",      "R1.MAL.OAA", "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("R1.MAL.OAA", "R2.MAL.OAA", "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("R2.MAL.OAA", "R3.MAL.OAA", "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("R3.MAL.OAA", "ID320",      "EFM 2", 61.99154)
flows2[nrow(flows2) + 1,] = list("ID282",      "R1.MAL.OAA", "EFM 4", 28.44492)
flows2[nrow(flows2) + 1,] = list("R1.MAL.OAA", "R2.MAL.OAA", "EFM 4", 28.44492)
flows2[nrow(flows2) + 1,] = list("R2.MAL.OAA", "R3.MAL.OAA", "EFM 4", 28.44492)
flows2[nrow(flows2) + 1,] = list("R3.MAL.OAA", "ID320",      "EFM 4", 28.44492)

flows2$substance <- paste0(flows2$substance, ": w = ", round(flows2$quantity, 1))

# PALETTE COLOURS
palette1 = data.frame(#
  substance = unique(flows2$substance),
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
nodes2$label <- ""

pdf(dir_output, width = 7, height = 7) # Set up PDF device
sankey(#
  #nodes2[c(3, 4, 5, 6, 7, 8, 1, 2, 13:18, 9, 10, 11, 12), ],
  nodes2,
  flows2,
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


