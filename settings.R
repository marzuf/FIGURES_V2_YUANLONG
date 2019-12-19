runFolder <-  "../v2_Yuanlong_Cancer_HiC_data_TAD_DA/"

pipFolder <- file.path(runFolder, "PIPELINE/OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

plotCex <- 1.4

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
