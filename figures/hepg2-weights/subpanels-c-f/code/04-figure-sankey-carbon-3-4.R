library("grid")
library("PantaRhei")

dir_inputs <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-3-4/01-for-R/"
dir_output <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-3-4/02-R-output-figures/sankey.pdf"

# Sankey diagram for glutamine bin 1 (5 carbon EFMs)
nodes3 <- read.csv(paste0(dir_inputs, "/nodes.csv"))
nodes3$x <- as.character(nodes3$x)
nodes3$y <- as.character(nodes3$y)
flows3 <- read.csv(paste0(dir_inputs, "/flows.csv"))

# Manually tweak nodes to prevent them from overlapping
capFirst <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
nodes3$label <- capFirst(nodes3$label)

nodes3[3, ]$label_pos <- "above" # ID320. OAA[m]
nodes3[3, ]$y <- "0.0"
nodes3[3, ]$label_pos <- "above"
nodes3[4, ]$x <- "ID320+4"       # ID137. Aspartate[m]
nodes3[4, ]$y <- "ID320+2.0"
nodes3[4, ]$label_pos <- "above"
nodes3[7, ]$x <- "ID320+4"       # ID154. Citrate[m]
nodes3[7, ]$y <- "ID320"
nodes3[7, ]$label_pos <- "above"
nodes3[5, ]$x <- "ID137+4"       # ID136. Aspartate[c]
nodes3[5, ]$y <- "ID137"
nodes3[5, ]$label_pos <- "above"
nodes3[8, ]$x <- "ID154+4"       # ID153. Citrate[c]
nodes3[8, ]$y <- "ID154"
nodes3[8, ]$label_pos <- "above"
nodes3[6, ]$x <- "ID320+12"      # ID319. OAA[c]
nodes3[6, ]$y <- "ID320+0.0"
nodes3[6, ]$label_pos <- "below"
nodes3[1, ]$x <- "ID319+6"       # ID281. Malate[c]
nodes3[1, ]$y <- "ID319"
nodes3[1, ]$label_pos <- "below"
nodes3[2, ]$x <- "ID319+3.0"     # ID282. Malate[m]
nodes3[2, ]$y <- "ID319+2.8"
nodes3[2, ]$label_pos <- "above"
nodes3[9, ]$x <- "ID320+0.0"     # ID214. Glutamate[c]
nodes3[9, ]$y <- "ID320-3.5"
nodes3[9, ]$label_pos <- "above"
nodes3[10, ]$x <- "ID281+0.0"     # ID215. Glutamate[m]
nodes3[10, ]$y <- "ID281-3.5-0.27"
nodes3[10, ]$label_pos <- "above"
nodes3[11, ]$x <- "ID154+0.0"     # ID125. AKG[m]
nodes3[11, ]$y <- "ID154-7.0"
nodes3[11, ]$label_pos <- "above"
nodes3[12, ]$x <- "ID319+0.0"     # ID124. AKG[c]
nodes3[12, ]$y <- "ID319-7.2"
nodes3[12, ]$label_pos <- "above"

nodes3[nrow(nodes3) + 1,] = c(".L1.MalC.MalM", "", "none", "ID281+3.0", "ID281+2", "up")
nodes3[nrow(nodes3) + 1,] = c(".L2.MalC.MalM", "", "none", "ID281", "ID281+7.0", "left")
nodes3[nrow(nodes3) + 1,] = c(".L3.MalC.MalM", "", "none", "ID282-3", "ID282+2", "down")

nodes3[nrow(nodes3) + 1,] = c(".L1.MalM.OAA", "", "none", "ID281-0.24", "ID281+5", "up")
nodes3[nrow(nodes3) + 1,] = c(".L2.MalM.OAA", "", "none", "ID282", "ID282+7", "left")
nodes3[nrow(nodes3) + 1,] = c(".L3.MalM.OAA", "", "none", "ID320-3", "ID320+2", "down")

nodes3[nrow(nodes3) + 1,] = c(".R1.Glut.AKG", "", "none", "ID215+2.5", "ID215-1.5", "down")
nodes3[nrow(nodes3) + 1,] = c(".R2.Glut.AKG", "", "none", "ID215", "ID215-6", "left")
nodes3[nrow(nodes3) + 1,] = c(".R3.Glut.AKG", "", "none", "ID125-2.5", "ID125-1.5", "up")

nodes3[nrow(nodes3) + 1,] = c(".R1.AKG.Glut", "", "none", "ID124+2.5", "ID124-1.5", "down")
nodes3[nrow(nodes3) + 1,] = c(".R2.AKG.Glut", "", "none", "ID124", "ID124-4", "left")
nodes3[nrow(nodes3) + 1,] = c(".R3.AKG.Glut", "", "none", "ID214-2.5", "ID214-1.5", "up")


flows3[nrow(flows3) + 1,] = list("ID281", "L1.MalC.MalM", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("ID281", "L1.MalC.MalM", "EFM 2", 52.372277776541964)
flows3[nrow(flows3) + 1,] = list("ID281", "L1.MalC.MalM", "EFM 4", 23.221825733891258)
flows3[nrow(flows3) + 1,] = list("L1.MalC.MalM", "L2.MalC.MalM", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("L1.MalC.MalM", "L2.MalC.MalM", "EFM 2", 52.372277776541964)
flows3[nrow(flows3) + 1,] = list("L1.MalC.MalM", "L2.MalC.MalM", "EFM 4", 23.221825733891258)
flows3[nrow(flows3) + 1,] = list("L2.MalC.MalM", "L3.MalC.MalM", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("L2.MalC.MalM", "L3.MalC.MalM", "EFM 2", 52.372277776541964)
flows3[nrow(flows3) + 1,] = list("L2.MalC.MalM", "L3.MalC.MalM", "EFM 4", 23.221825733891258)
flows3[nrow(flows3) + 1,] = list("L3.MalC.MalM", "ID282", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("L3.MalC.MalM", "ID282", "EFM 2", 52.372277776541964)
flows3[nrow(flows3) + 1,] = list("L3.MalC.MalM", "ID282", "EFM 4", 23.221825733891258)
flows3 <- flows3[-c(1, 7, 17), ]

flows3[nrow(flows3) + 1,] = list("ID282", "L1.MalM.OAA", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("ID282", "L1.MalM.OAA", "EFM 2", 52.372277776541964)
flows3[nrow(flows3) + 1,] = list("L1.MalM.OAA", "L2.MalM.OAA", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("L1.MalM.OAA", "L2.MalM.OAA", "EFM 2", 52.372277776541964)
flows3[nrow(flows3) + 1,] = list("L2.MalM.OAA", "L3.MalM.OAA", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("L2.MalM.OAA", "L3.MalM.OAA", "EFM 2", 52.372277776541964)
flows3[nrow(flows3) + 1,] = list("L3.MalM.OAA", "ID320", "EFM 1", 60.63570423043401)
flows3[nrow(flows3) + 1,] = list("L3.MalM.OAA", "ID320", "EFM 2", 52.372277776541964)
flows3 <- flows3[-c(1, 6), ]

flows3[nrow(flows3) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 3", 46.59569606300289)
flows3[nrow(flows3) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 5", 16.07600833312059)
flows3[nrow(flows3) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 3", 46.59569606300289)
flows3[nrow(flows3) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 5", 16.07600833312059)
flows3[nrow(flows3) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 3", 46.59569606300289)
flows3[nrow(flows3) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 5", 16.07600833312059)
flows3[nrow(flows3) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 3", 46.59569606300289)
flows3[nrow(flows3) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 5", 16.07600833312059)
flows3 <- flows3[-c(10, 15), ]

flows3[nrow(flows3) + 1,] = list("ID124", "R1.AKG.Glut", "EFM 3", 46.59569606300289)
flows3[nrow(flows3) + 1,] = list("R1.AKG.Glut", "R2.AKG.Glut", "EFM 3", 46.59569606300289)
flows3[nrow(flows3) + 1,] = list("R2.AKG.Glut", "R3.AKG.Glut", "EFM 3", 46.59569606300289)
flows3[nrow(flows3) + 1,] = list("R3.AKG.Glut", "ID214", "EFM 3", 46.59569606300289)
flows3 <- flows3[-c(11), ]

flows3$substance <- paste0(flows3$substance, ": w = ", round(flows3$quantity, 1))

# PALETTE COLOURS
palette3 = data.frame(#
    substance = unique(flows3$substance),
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
title3 <- "Top 5 ranked AEFMs of glutamine carbons 3 or 4"
attr(title3, "gp") <- grid::gpar(#
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
    nodes3,
    flows3,
    palette3,
    max_width   = 0.1,
    rmin        = 0.5,
    node_style  = ns,
    page_margin = c(0, 0.1, 0, 0.1), # left bottom right top
    legend      = gpar(col = text_node_col, ncols = 4, fontsize = font_master),
    title       = title3,
    grill       = FALSE,
    copyright   = ""
)
dev.off()


