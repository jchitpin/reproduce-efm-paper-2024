library("grid")
library("PantaRhei")

dir_inputs <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-1-2/01-for-R/"
dir_output <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-1-2/02-R-output-figures/sankey.pdf"

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
nodes1[4, ]$x <- "ID320+4"       # ID137. Aspartate[m]
nodes1[4, ]$y <- "ID320+2.0"
nodes1[4, ]$label_pos <- "above"
nodes1[7, ]$x <- "ID320+4"       # ID154. Citrate[m]
nodes1[7, ]$y <- "ID320"
nodes1[7, ]$label_pos <- "above"
nodes1[5, ]$x <- "ID137+4"       # ID136. Aspartate[c]
nodes1[5, ]$y <- "ID137"
nodes1[5, ]$label_pos <- "above"
nodes1[8, ]$x <- "ID154+4"       # ID153. Citrate[c]
nodes1[8, ]$y <- "ID154"
nodes1[8, ]$label_pos <- "above"
nodes1[6, ]$x <- "ID320+12"      # ID319. OAA[c]
nodes1[6, ]$y <- "ID320+0.27"
nodes1[6, ]$label_pos <- "above"
nodes1[16, ]$x <- "ID153+2"      # ID262. Isocitrate[c]
nodes1[16, ]$y <- "ID153-4.0"
nodes1[16, ]$label_pos <- "below"
nodes1[12, ]$x <- "ID262+4"      # ID124. AKG[c]
nodes1[12, ]$y <- "ID262-2.0"
nodes1[12, ]$label_pos <- "above"
nodes1[9, ]$x <- "ID124+3"       # ID214. Glutamate[c]
nodes1[9, ]$y <- "ID124"
nodes1[9, ]$label_pos <- "below"
nodes1[10, ]$x <- "ID214+3"      # ID215. Glutamate[m]
nodes1[10, ]$y <- "ID214"
nodes1[10, ]$label_pos <- "above"
nodes1[nrow(nodes1) + 1,] = c(".R1.Glut.AKG", "", "none", "ID215+1.5", "ID215-1.5", "down")
nodes1[nrow(nodes1) + 1,] = c(".R2.Glut.AKG", "", "none", "ID262", "ID262-6", "left")
nodes1[nrow(nodes1) + 1,] = c(".R3.Glut.AKG", "", "none", "ID262-9.5", "ID262-3.5", "up")
nodes1[11, ]$x <- "ID124-12"     # ID125. AKG[m]
nodes1[11, ]$y <- "ID124"
nodes1[11, ]$label_pos <- "above"
nodes1[13, ]$x <- "ID153-3"      # ID390. Succinyl-CoA[m]
nodes1[13, ]$y <- "ID153-2.5"
nodes1[13, ]$label_pos <- "above"
nodes1[14, ]$x <- "ID390+2"      # ID389. Succinate[m]
nodes1[14, ]$y <- "ID390"
nodes1[14, ]$label_pos <- "below"
nodes1[1, ]$x <- "ID319+6"       # ID281. Malate[c]
nodes1[1, ]$y <- "ID319"
nodes1[1, ]$label_pos <- "below"
nodes1[15, ]$x <- "ID389+5"      # ID200. Fumarate[m]
nodes1[15, ]$y <- "ID389"
nodes1[15, ]$label_pos <- "below"
nodes1[2, ]$x <- "ID319+3.0"     # ID282. Malate[m]
nodes1[2, ]$y <- "ID319+2.39"
nodes1[2, ]$label_pos <- "above"
nodes1[nrow(nodes1) + 1,] = c(".L1.MalC.MalM", "", "none", "ID281+3.0", "ID281+2", "up")
nodes1[nrow(nodes1) + 1,] = c(".L2.MalC.MalM", "", "none", "ID281", "ID281+6.5", "left")
nodes1[nrow(nodes1) + 1,] = c(".L3.MalC.MalM", "", "none", "ID282-3", "ID282+2", "down")
nodes1[nrow(nodes1) + 1,] = c(".L1.MalM.OAA", "", "none", "ID281-0.24", "ID281+4", "up")
nodes1[nrow(nodes1) + 1,] = c(".L2.MalM.OAA", "", "none", "ID282", "ID282+7", "left")
nodes1[nrow(nodes1) + 1,] = c(".L3.MalM.OAA", "", "none", "ID320-3", "ID320+2", "down")

flows1[nrow(flows1) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 3", 73.07798877271968)
flows1[nrow(flows1) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 4", 48.56112282771928)
flows1[nrow(flows1) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 3", 73.07798877271968)
flows1[nrow(flows1) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 4", 48.56112282771928)
flows1[nrow(flows1) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 3", 73.07798877271968)
flows1[nrow(flows1) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 4", 48.56112282771928)
flows1[nrow(flows1) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 3", 73.07798877271968)
flows1[nrow(flows1) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 4", 48.56112282771928)
flows1 <- flows1[-c(14, 18), ]

flows1[nrow(flows1) + 1,] = list("ID281", "L1.MalC.MalM", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("ID281", "L1.MalC.MalM", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("ID281", "L1.MalC.MalM", "EFM 5", 36.396198953756176)
flows1[nrow(flows1) + 1,] = list("L1.MalC.MalM", "L2.MalC.MalM", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("L1.MalC.MalM", "L2.MalC.MalM", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("L1.MalC.MalM", "L2.MalC.MalM", "EFM 5", 36.396198953756176)
flows1[nrow(flows1) + 1,] = list("L2.MalC.MalM", "L3.MalC.MalM", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("L2.MalC.MalM", "L3.MalC.MalM", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("L2.MalC.MalM", "L3.MalC.MalM", "EFM 5", 36.396198953756176)
flows1[nrow(flows1) + 1,] = list("L3.MalC.MalM", "ID282", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("L3.MalC.MalM", "ID282", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("L3.MalC.MalM", "ID282", "EFM 5", 36.396198953756176)
flows1 <- flows1[-c(2, 8, 21), ]

flows1[nrow(flows1) + 1,] = list("ID282", "L1.MalM.OAA", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("ID282", "L1.MalM.OAA", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("ID282", "L1.MalM.OAA", "EFM 4", 48.56112282771928)
flows1[nrow(flows1) + 1,] = list("L1.MalM.OAA", "L2.MalM.OAA", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("L1.MalM.OAA", "L2.MalM.OAA", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("L1.MalM.OAA", "L2.MalM.OAA", "EFM 4", 48.56112282771928)
flows1[nrow(flows1) + 1,] = list("L2.MalM.OAA", "L3.MalM.OAA", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("L2.MalM.OAA", "L3.MalM.OAA", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("L2.MalM.OAA", "L3.MalM.OAA", "EFM 4", 48.56112282771928)
flows1[nrow(flows1) + 1,] = list("L3.MalM.OAA", "ID320", "EFM 1", 95.03598813297019)
flows1[nrow(flows1) + 1,] = list("L3.MalM.OAA", "ID320", "EFM 2", 82.08449514089918)
flows1[nrow(flows1) + 1,] = list("L3.MalM.OAA", "ID320", "EFM 4", 48.56112282771928)
flows1 <- flows1[-c(1, 6, 24), ]

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
title1 <- "Top 5 ranked AEFMs of glutamine carbons 1 or 2"
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

pdf(dir_output, width = 7, height = 7) # Set up PDF device
sankey(#
    nodes1,
    flows1,
    palette1,
    max_width   = 0.1,
    rmin        = 0.5,
    node_style  = ns,
    page_margin = c(0, 0.1, 0, 0.1), # left bottom right top
    legend      = gpar(col = text_node_col, ncols = 4, fontsize = font_master),
    title       = title1,
    grill       = FALSE, # TRUE for debugging
    copyright   = ""
)
dev.off()


