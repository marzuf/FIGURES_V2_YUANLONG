# Rscript fcc_density_heatmap_v2.R

# each column => dataset; color-coded by density

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "fcc_density_heatmap.R"
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

outFolder <- "FCC_DENSITY_HEATMAP_V2" 
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files(pipOutFolder)

# all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
# all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds)) ]
# all_hicds <- all_hicds[grepl("ENCSR489OCU_NCI-H460_40kb", all_hicds)]
# all_hicds = "ENCSR489OCU_NCI-H460_40kb"
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMMIDPOSDISC" , "RANDOMMIDPOSSTRICT", "RANDOMSHIFT", "PERMUTG2T")


resolution <- 0.01

buildData <- TRUE

if(buildData){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, " \n")
      
      fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")
      if(!file.exists(fcc_file)) {
        fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCConlyObs_runAllDown", "all_obs_prodSignedRatio.Rdata")  
      }
      if(!file.exists(fcc_file)) return(NULL)
      stopifnot(file.exists(fcc_file))
      all_fcc <- as.numeric(get(load(fcc_file)))
      # if(!file.exists(fcc_file)) {
      #   fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCConlyObs_runAllDown", "all_obs_prodSignedRatio.Rdata")  
      # }
      
      data.frame(
        hicds = hicds,
        exprds=exprds,
        FCC = all_fcc,
        stringsAsFactors = FALSE
      )
    }
    hicds_dt
  }
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- get(load(outFile))
}  

all_result_dt$dataset <- file.path(all_result_dt$hicds, all_result_dt$exprds)


all_result_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                                gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                     gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
                                          gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                               gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                                    gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_result_dt$hicds))))))
all_result_dt$hicds_lab[! all_result_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"

all_types <- unique(all_result_dt$hicds_lab)
a_t = all_types[1]
for(a_t in all_types) {
  
  sub_dt <- all_result_dt[all_result_dt$hicds_lab == a_t,]

  # retrieve dataset order
  # if(a_t == "OBSERVED") {
  #   tmp <- get(load("../FIGURES_V3_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))
  #   tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
  #   ds_levels <- file.path(tmp$hicds, tmp$exprds)
  #   
  # } else {
  #   auc_file <- file.path("..", "FIGURES_V3_YUANLONG", paste0("BARPLOT_FCC_AUC_RATIO_", a_t),  "all_dt.Rdata")
  #   if(!file.exists(auc_file)) next
  #   tmp <- get(load(auc_file))
  #   tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
  #   ds_levels <- file.path(tmp$hicds, tmp$exprds)
  #   
  # }
  
  # to have something for g2t permut
  # retrieve dataset order
  
  tmp <- get(load("../FIGURES_V3_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))
  tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
  ds_levels <- file.path(tmp$hicds, tmp$exprds)
    
  if(a_t != "OBSERVED") {
    auc_file <- file.path("..", "FIGURES_V3_YUANLONG", paste0("BARPLOT_FCC_AUC_RATIO_", a_t),  "all_dt.Rdata")
    if(file.exists(auc_file)) {
    tmp <- get(load(auc_file))
    tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
    ds_levels <- file.path(tmp$hicds, tmp$exprds)
    } else {
      
      sub_dt$dataset <- gsub("_RANDOMSHIFT_40kb", "_40kb", 
                                      gsub("_RANDOMMIDPOSDISC_40kb", "_40kb",
                                           gsub("_RANDOMMIDPOSSTRICT_40kb", "_40kb",
                                                gsub("_RANDOMNBRGENES_40kb", "_40kb",
                                                     gsub("_PERMUTG2T_40kb", "_40kb", 
                                                          gsub("_RANDOMMIDPOS_40kb", "_40kb", sub_dt$dataset))))))
      
      
    }
  }
  
  
  

  sub_dt$dataset <- factor(sub_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(sub_dt$dataset))
  sub_dt$heatmap_x <- as.numeric(sub_dt$dataset) 
  dsCols <- all_cols[all_cmps[basename(ds_levels)]]
  

  sub_density_dt <- do.call(rbind, by(sub_dt, sub_dt$dataset, function(x) {
    hicds <- unique(x$hicds)
    stopifnot(length(hicds) == 1)
    exprds <- unique(x$exprds)
    stopifnot(length(exprds) == 1)
    ds_dt <- density(x$FCC)
    data.frame(
      hicds=hicds,
      exprds=exprds,
      density_x  = ds_dt$x,
      density_y  = ds_dt$y,
      stringsAsFactors = FALSE
    )
  }))
  sub_density_dt$dataset <- file.path(sub_density_dt$hicds, sub_density_dt$exprds)
  
  if(!sub_density_dt$dataset %in% ds_levels) {
    sub_density_dt$dataset <- gsub("_RANDOMSHIFT_40kb", "_40kb", 
                           gsub("_RANDOMMIDPOSDISC_40kb", "_40kb",
                                gsub("_RANDOMMIDPOSSTRICT_40kb", "_40kb",
                                     gsub("_RANDOMNBRGENES_40kb", "_40kb",
                                          gsub("_PERMUTG2T_40kb", "_40kb", 
                                               gsub("_RANDOMMIDPOS_40kb", "_40kb", sub_density_dt$dataset))))))
    
  }
  
  sub_density_dt$dataset <- factor(sub_density_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(sub_density_dt$dataset))
  sub_density_dt$heatmap_x <- as.numeric(sub_density_dt$dataset)

  nDS <- length(unique(as.character(sub_dt$dataset)))

  # new_density_x <- seq(min(sub_density_dt$density_x),max(sub_density_dt$density_x),by=resolution)
  new_density_x <- seq(-1,1,by=resolution)

  plot_dt <- foreach(i = 1:nDS, .combine='rbind') %dopar% {
    sub_dt <- sub_density_dt[sub_density_dt$heatmap_x == i,]
    new_density <- approx(x=sub_dt$density_x,  y=sub_dt$density_y, xout = new_density_x, rule=2)
    data.frame(
      dataset=ds_levels[i],
      density_x = new_density[["x"]],
      density_y = new_density[["y"]],
      stringsAsFactors = FALSE
    )
  }
  plot_dt$dataset <- factor(plot_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(plot_dt$dataset))

  density_plot <- ggplot(plot_dt, aes(x = dataset, y = density_x, fill = density_y))+
    geom_tile() +
  ggtitle(paste0("FCC score distribution - ", a_t),
          subtitle = paste0("all datasets (n=", nDS, ")")) +
    scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
    scale_y_continuous(name="FCC score",
                       breaks = scales::pretty_breaks(n = 20),  expand = c(0, 0))+
    labs(fill = "Density")+
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
  
  outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_densityheatmap.", plotType))
  ggsave(density_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0("density_plot_", a_t, ".Rdata"))
  save(plot_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  sub_dt <- all_result_dt[all_result_dt$hicds_lab == a_t,]
  sub_dt$FCC[sub_dt$hicds == dirname(ds_levels[1]) & sub_dt$exprds == basename(ds_levels[1])] <- sub_dt$FCC[sub_dt$hicds == dirname(ds_levels[1]) & sub_dt$exprds == basename(ds_levels[1])] - 0.5
    
  sub_density_dt <- do.call(rbind, by(sub_dt, sub_dt$dataset, function(x) {
    hicds <- unique(x$hicds)
    stopifnot(length(hicds) == 1)
    exprds <- unique(x$exprds)
    stopifnot(length(exprds) == 1)
    ds_dt <- density(x$FCC)
    data.frame(
      hicds=hicds,
      exprds=exprds,
      density_x  = ds_dt$x,
      density_y  = ds_dt$y,
      stringsAsFactors = FALSE
    )
  }))
  
  sub_density_dt$dataset <- file.path(sub_density_dt$hicds, sub_density_dt$exprds)
  
  if(!sub_density_dt$dataset %in% ds_levels) {
    sub_density_dt$dataset <- gsub("_RANDOMSHIFT_40kb", "_40kb", 
                                   gsub("_RANDOMMIDPOSDISC_40kb", "_40kb",
                                        gsub("_RANDOMMIDPOSSTRICT_40kb", "_40kb",
                                             gsub("_RANDOMNBRGENES_40kb", "_40kb",
                                                  gsub("_PERMUTG2T_40kb", "_40kb", 
                                                       gsub("_RANDOMMIDPOS_40kb", "_40kb", sub_density_dt$dataset))))))
    
  }
  
  
  
  sub_density_dt$dataset <- factor(sub_density_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(sub_density_dt$dataset))
  sub_density_dt$heatmap_x <- as.numeric(sub_density_dt$dataset) 
  
  
  nDS <- length(unique(as.character(sub_dt$dataset)))
  
  new_density_x <- seq(min(sub_density_dt$density_x),max(sub_density_dt$density_x),by=resolution)
  
  plot_dt <- foreach(i = 1:nDS, .combine='rbind') %dopar% {
    sub_dt <- sub_density_dt[sub_density_dt$heatmap_x == i,]
    new_density <- approx(x=sub_dt$density_x,  y=sub_dt$density_y, xout = new_density_x, rule=2)
    data.frame(
      dataset=ds_levels[i],
      density_x = new_density[["x"]],
      density_y = new_density[["y"]],
      stringsAsFactors = FALSE
    )
  }
  plot_dt$dataset <- factor(plot_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(plot_dt$dataset))
  
  density_plot <- ggplot(plot_dt, aes(x = dataset, y = density_x, fill = density_y))+ 
    geom_tile() +
    ggtitle(paste0("FCC score distribution - ", a_t),   
            subtitle = paste0("all datasets (n=", nDS, ")")) +
    scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
    scale_y_continuous(name="FCC score",
                       breaks = scales::pretty_breaks(n = 20),  expand = c(0, 0))+
    labs(fill = "Density")+
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
  
  outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_densityheatmap_checkplot.", plotType))
  ggsave(density_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  

  
  
  
}


