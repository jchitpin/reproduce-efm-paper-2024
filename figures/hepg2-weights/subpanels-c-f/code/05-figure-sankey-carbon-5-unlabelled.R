library("grid")
library("PantaRhei")

dir_inputs <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-5/01-for-R/"
dir_output <- "figures/hepg2-weights/subpanels-c-f/glutamine/carbon-5/02-R-output-figures/sankey-unlabelled.pdf"

# Sankey diagram for glutamine bin 1 (5 carbon EFMs)
nodes5 <- read.csv(paste0(dir_inputs, "/nodes.csv"))
nodes5$x <- as.character(nodes5$x)
nodes5$y <- as.character(nodes5$y)
flows5 <- read.csv(paste0(dir_inputs, "/flows.csv"))
nodes5$ID[11] <- "ID0"

# Manually tweak nodes to prevent them from overlapping
capFirst <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
nodes5$label <- capFirst(nodes5$label)

nodes5[5, ]$label_pos <- "above"   # ID219. Glutamine[s]
nodes5[5, ]$y <- "0.0"
nodes5[5, ]$label_pos <- "above"
nodes5[6, ]$x <- "ID219+10"        # ID217. Glutamine[c]
nodes5[6, ]$y <- "ID219+0.0"
nodes5[6, ]$label_pos <- "above"
nodes5[11, ]$x <- "ID217+10"       # ID0. Sink
nodes5[11, ]$y <- "ID217+0.0"
nodes5[11, ]$label_pos <- "above"
nodes5[4, ]$x <- "ID219+3.0"       # ID125. AKG[m]
nodes5[4, ]$y <- "ID219-2.55"
nodes5[4, ]$label_pos <- "above"
nodes5[10, ]$x <- "ID0-3.0"        # ID158. CO2[s]
nodes5[10, ]$y <- "ID0-2.0"
nodes5[10, ]$label_pos <- "above"
nodes5[9, ]$x <- "ID158-3.0"       # ID156. CO2[c]
nodes5[9, ]$y <- "ID158-0.0"
nodes5[9, ]$label_pos <- "above"
nodes5[8, ]$x <- "ID156-6.0"       # ID157. CO2[m]
nodes5[8, ]$y <- "ID156-0.0"
nodes5[8, ]$label_pos <- "above"
nodes5[7, ]$x <- "ID157+4.0"       # ID218. Glutamine[m]
nodes5[7, ]$y <- "ID157-2.0"
nodes5[7, ]$label_pos <- "above"
nodes5[2, ]$x <- "ID157+3.29"       # ID214. Glutamate[c]
nodes5[2, ]$y <- "ID157-4.5"
nodes5[2, ]$label_pos <- "above"
nodes5[3, ]$x <- "ID214+4.0"       # ID215. Glutamate[m]
nodes5[3, ]$y <- "ID214+0.4"
nodes5[3, ]$label_pos <- "above"
nodes5[1, ]$x <- "ID125+5.0"       # ID124. AKG[c]
nodes5[1, ]$y <- "ID214-0.1"
nodes5[1, ]$label_pos <- "above"

nodes5[nrow(nodes5) + 1,] = c(".R1.Glut.AKG", "", "none", "ID215+3.0", "ID215-2", "down")
nodes5[nrow(nodes5) + 1,] = c(".R2.Glut.AKG", "", "none", "ID215", "ID215-6.0", "left")
nodes5[nrow(nodes5) + 1,] = c(".R3.Glut.AKG", "", "none", "ID125-3.0", "ID125-2", "up")

flows5[nrow(flows5) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 1", 30.322355112618716)
flows5[nrow(flows5) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 2", 18.603762568203496)
flows5[nrow(flows5) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 3", 10.463223200773367)
flows5[nrow(flows5) + 1,] = list("ID215", "R1.Glut.AKG", "EFM 5", 7.0041805433506665)
flows5[nrow(flows5) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 1", 30.322355112618716)
flows5[nrow(flows5) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 2", 18.603762568203496)
flows5[nrow(flows5) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 3", 10.463223200773367)
flows5[nrow(flows5) + 1,] = list("R1.Glut.AKG", "R2.Glut.AKG", "EFM 5", 7.0041805433506665)
flows5[nrow(flows5) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 1", 30.322355112618716)
flows5[nrow(flows5) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 2", 18.603762568203496)
flows5[nrow(flows5) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 3", 10.463223200773367)
flows5[nrow(flows5) + 1,] = list("R2.Glut.AKG", "R3.Glut.AKG", "EFM 5", 7.0041805433506665)
flows5[nrow(flows5) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 1", 30.322355112618716)
flows5[nrow(flows5) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 2", 18.603762568203496)
flows5[nrow(flows5) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 3", 10.463223200773367)
flows5[nrow(flows5) + 1,] = list("R3.Glut.AKG", "ID125", "EFM 5", 7.0041805433506665)


flows5 <- flows5[-c(3, 8, 14, 20), ]

flows5$substance <- paste0(flows5$substance, ": w = ", round(flows5$quantity, 1))

# PALETTE COLOURS
palette5 = data.frame(#
    substance = unique(flows5$substance),
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
title5 <- ""
attr(title5, "gp") <- grid::gpar(#
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
nodes5$label <- ""

pdf(dir_output, width = 10, height = 7) # Set up PDF device
sankey(#
    nodes5,
    flows5,
    palette5,
    max_width   = 0.1,
    rmin        = 0.5,
    node_style  = ns,
    page_margin = c(0, 0.1, 0, 0.1), # left bottom right top
    #legend      = gpar(col = text_node_col, ncols = 4, fontsize = font_master),
    title       = title5,
    grill       = FALSE,
    copyright   = ""
)
dev.off()


