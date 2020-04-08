
# V2
# labsymbol <- "\u25CF" # circle
# V3
labsymbol <- "\u25A0" # squares

runFolder <-  "../v2_Yuanlong_Cancer_HiC_data_TAD_DA/"

pipFolder <- file.path(runFolder, "PIPELINE/OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

plotCex <- 1.4
axisCex <- 1.4
labCex <- 1.4
mainCex <- 1.4
subCex <- 1.2

myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

myHeightGG <- 7
myWidthGG <- 7


step8fcc_folder <- "8cOnlyFCC_runAllDown"

fontFamily <- "Hershey"

# ggsci package palettes:
#pal_npg()(100) # 10
#pal_aaas()(100) # 10
#pal_nejm()(100) # 8
#pal_lancet()(100) # 9
#pal_jama()(100) # 7
#pal_jco()(100) # 10
#pal_ucscgb()(100) # 26
#pal_d3()(100) # 10
#pal_locuszoom()(100) # 7
#pal_igv()(100)  # 51
#pal_cosmic()(100) # NF
#pal_uchicago()(100) # 9
#pal_startrek()(100) # 7
#pal_tron()(100) # 7
#pal_futurama()(100) # 12
#pal_rickandmorty()(100) # 12
#pal_simpsons()(100) # 16
#pal_gsea()(100) # 12
#pal_material()(100) # 10

tad_signif_col <- "dodgerblue3"
gene_signif_col <- "firebrick3"

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
# V2:
#all_cols[all_cols == "red"] <- "brown3"  # wt vs mut
#all_cols[all_cols == "blue"] <- "darkblue" # norm vs tum
#all_cols[all_cols == "green"] <- "forestgreen" # subtypes
# V3:
#all_cols[all_cols == "red"] <- "brown3"  # wt vs mut
#all_cols[all_cols == "blue"] <- "darkblue" # norm vs tum
#all_cols[all_cols == "green"] <- "forestgreen" # subtypes


all_cols[all_cols == "red"] <- "violetred" #"chocolate"  # wt vs mut
all_cols[all_cols == "blue"] <- "slateblue" # norm vs tum
all_cols[all_cols == "green"] <- "slategray" # "yellow3" # subtypes

options(save.defaults = list(version=2), scipen=100)

geneSignifThresh <- 0.01
tadSignifThresh <- 0.01



