# Rscript cmp_fcc_heatmapDiff.R

# each column => dataset; color-coded by density

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "cmp_fcc_heatmapDiff_permut.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
require(ggpubr)
require(ggplot2)
registerDoMC(40)

plotType <- "png"
myHeight <- 400
myWidth <- 400
myHeightGG <- 7
myWidthGG <- 14

plotCex <- 1.4

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

pipFolder<- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/")
stopifnot(dir.exists(pipFolder))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CMP_FCC_HEATMAPDIFF_PERMUT" 
dir.create(outFolder, recursive = TRUE)

auc_file <- file.path("..", "FIGURES_V3_YUANLONG", paste0("BARPLOT_FCC_AUC_RATIO"),  "all_dt.Rdata")
tmp <- get(load(auc_file))
tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
ds_levels <- file.path(tmp$hicds, tmp$exprds)
dsCols <- all_cols[all_cmps[basename(ds_levels)]]

obsDT <- get(load("FCC_DENSITY_HEATMAP_V2/density_plot_OBSERVED.Rdata"))

rd_type ="RANDOMMIDPOSSTRICT"
for(rd_type in c("RANDOMMIDPOSSTRICT", "PERMUTG2T")) {
  
  permutDT <- get(load(paste0("FCC_DENSITY_HEATMAP_V2/density_plot_", rd_type, ".Rdata")))
  permutDT$dataset <- gsub("_RANDOMMIDPOSSTRICT_40kb", "_40kb", permutDT$dataset)
  permutDT$dataset <- gsub("_PERMUTG2T_40kb", "_40kb", permutDT$dataset)
  
  merged_dt <- merge(obsDT, permutDT, by=c("dataset", "density_x"), suffixes = c("_obs", "_permut"))
  stopifnot(!is.na(merged_dt))
  
  merged_dt$density_y_ratioObsPermut <- log10(merged_dt$density_y_obs/merged_dt$density_y_permut)
  
  nDS <- length(unique(merged_dt$dataset))
  
  merged_dt$dataset <- factor(merged_dt$dataset, levels = ds_levels)
  stopifnot(!is.na(merged_dt$dataset))
  
  density_plot <- ggplot(merged_dt, aes(x = dataset, y = density_x, fill = density_y_ratioObsPermut))+ 
    geom_tile() +
    ggtitle(paste0("Ratio FCC score distribution - obs/permut"),   
            subtitle = paste0(rd_type, " - all datasets (n=", nDS, ")")) +
    scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
    scale_y_continuous(name="FCC score",
                       breaks = scales::pretty_breaks(n = 20),  expand = c(0, 0))+
    labs(fill = "Density\nobs/permut\n[log10]")+
    scale_fill_gradient( high="red", low="blue", na.value = "white" )  +
    theme(
      axis.text.x = element_text(colour = dsCols, size=12),
      axis.text.y= element_text(colour = "black", size=12),
      axis.title.x = element_text(colour = "black", size=14, face="bold"),
      axis.title.y = element_text(colour = "black", size=14, face="bold"),
      plot.title = element_text(hjust=0.5, size=16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
      panel.background = element_rect(fill = "transparent")
      # legend.background =  element_rect()
    )
  
  density_plot <- density_plot + geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3)
  
  outFile <- file.path(outFolder, paste0("ratio_FCC_score_dist_allDS_obs_permut_", rd_type, "_densityheatmap.", plotType))
  ggsave(density_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))



  
}

