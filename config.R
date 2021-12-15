options(bitmapType = "cairo")

# mae
config$id.prefix <- "Animal"

# kegg files
dir.kegg <- "00_dataraw/kegg/"
file.kegg.brite <- paste0(dir.kegg, "br08901.keg")
file.kegg.mmu.cpd <- paste0(dir.kegg, "pathway_cpd_mmu.txt")
file.kegg.mmu.names <- paste0(dir.kegg, "names_mmu.txt")
file.kegg.mmu.list <- paste0(dir.kegg, "pathway_list_mmu.txt")
file.kegg.mmu.gene <- paste0(dir.kegg, "pathway_gene_mmu.txt")
file.kegg.cpd.names <- paste0(dir.kegg, "names_cpd.txt")

# FELLA analysis
fella.db <- "00_dataraw/fella/2020-12_mmu-metabpathwaysonly"
# if (!dir.exists(fella.db)) {
#   g <- FELLA::buildGraphFromKEGGREST(organism = "mmu", filter.path = "01100")
#   FELLA::buildDataFromGraph(
#     g, databaseDir = fella.db, 
#     internalDir = FALSE, 
#     matrices = c("hypergeom", "diffusion"), 
#     normality = c("hypergeom", "diffusion"), 
#     niter = 1000)  
# }


# se metab
dir.metab <- "01_se_metab/"
file.metab.se <- paste0(dir.metab, "se_metab_lung.rds")

# se lipid
dir.lipid <- "02_se_lipid/"
file.lipid.se <- paste0(dir.lipid, "se_lipid_lung.rds")

# id.prefix <- "Mouse"
dir.mae <- "03_mae"
out.mae <- paste0(dir.mae, "/mae.rds")

# descriptive plots
dir.descriptive <- "10_descriptive"


# differential analysis
fc.ref <- 0
q.ref <- .05

dir.differential <- "20_differential_analysis"
df.fc.trt <- paste0(dir.differential, "/df_logfcs_pvalues_treatment.csv")
df.fc.age <- paste0(dir.differential, "/df_logfcs_pvalues_age.csv")


dir.pathways <- "21_pathway_analysis"
dir.fella <- "25_network_analysis"

# plots
dir.diffplots <- "30_differential_plots"
dir.diffplotsage <- "31_differential_plots_age"

se2name <- c(
  se.lip.lung = "Lipids lung", 
  se.met.lung = "Metabolites lung"
)

se2nameshort <- c(
  se.met.lung = "Metabolites", 
  se.lip.lung = "Lipids"
)


# colors and palettes
# use always names in alphabetical order for simplicity
# 
# to get hex codes
# v2hex <- function(v) dput(setNames(gplots::col2hex(v), names(v)))
# v2hex(col.briteL1)


# Order here is important! Order of plotting labels
col.age <- c(
  Young = "#E5B61D", 
  Old = "#BC7EC9"
)

# with alpha
col.agetrt <- c(
  "Young NaCl" = "#F3E4B2", 
  "Young Bleo" = "#E5B61D", 
  "Old NaCl" = "#E8C0E8", 
  "Old Bleo" = "#BC7EC9"
)

# without alpha
col.agetrt.noalpha <- c(
  "Young NaCl" = "#E5B61D", 
  "Young Bleo" = "#E5B61D", 
  "Old NaCl" = "#BC7EC9", 
  "Old Bleo" = "#BC7EC9"
)

# labels&colors for significance
v.oldyoung.labels <- c("00" = "n.s.", "01" = "Young", "10" = "Old", "11" = "Both")
v.oldyoung.colors <- setNames(
  c("gray85", config$col.age["Young"], config$col.age["Old"], "#C09393"), 
  v.oldyoung.labels
)




col.omic <- setNames(
  RColorBrewer::brewer.pal(4, "Set2")[3:4], 
  se2nameshort
)
col.omicLM <- col.omic

# order matters
col.trt <- c(
  NaCl = "gray70", 
  Bleo = "gray20"
)

# from set1
col.trt2 <- c(
  NaCl = "#377EB8", 
  Bleo = "#E41A1C"
)

v.bleonacl.labels <- c("00" = "n.s.", "01" = "NaCl", "10" = "Bleo", "11" = "Both")
v.bleonacl.colors <- setNames(
  c("gray85", config$col.trt2["Bleo"], config$col.trt2["NaCl"], "darkslategray"), 
  v.bleonacl.labels
)

col.briteL1 <- setNames(
  RColorBrewer::brewer.pal(7, "Set3"), 
  c("Metabolism", "Genetic Information Processing", 
    "Environmental Information Processing", "Cellular Processes", 
    "Organismal Systems", "Human Diseases", 
    "Drug Development")
)

pch.agetrt <- c("Young NaCl" = 1, 
                "Young Bleo" = 16, 
                "Old NaCl" = 2, 
                "Old Bleo" = 17)

pardefault <- par()

# theme for plots
theme_set(theme_bw())

ggtheme <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 7, family = "ArialMT"), 
  axis.title = ggplot2::element_text(size = 6, family = "ArialMT"),
  axis.text = ggplot2::element_text(size = 5, family = "ArialMT"),
  text = ggplot2::element_text(size = 5, family = "ArialMT"), 
  panel.grid.minor = ggplot2::element_line(size = .2), 
  panel.grid.major = ggplot2::element_line(size = .35))


# plot dimensions
pca.inches <- 5/2.54
