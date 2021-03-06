#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

# Rscript fcc_wave_plot.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc
# Rscript fcc_wave_plot.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad
# FCC AUC top1                           LG1_40kb      TCGAluad_mutKRAS_mutEGFR 1.531834
# FCC AUC top12             ENCSR079VIJ_G401_40kb            TCGAkich_norm_kich 1.425827
# FCC AUC top32              GSE105381_HepG2_40kb         TCGAlihc_wt_mutCTNNB1 1.280421
# FCC AUC top58         ENCSR312KHQ_SK-MEL-5_40kb       TCGAskcm_lowInf_highInf 1.160567


# Rscript fcc_wave_plot.R LG1_40kb TCGAluad_mutKRAS_mutEGFR
# Rscript fcc_wave_plot.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich
# Rscript fcc_wave_plot.R GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1
# Rscript fcc_wave_plot.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript fcc_wave_plot.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1

################################################################################

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)
plotCex <- 1.4

require(flux)

# set colors:
polygonPermutCol <- rgb(188/255,188/255,188/255, 0.3)
meanPermutCol <- rgb(135/255,135/255,135/255, 0.3)
# qt95PermutCol <- rgb(140/255,140/255,140/255, 0.3)
qt95PermutCol <- rgb(140/255,140/255,140/255)
lwdObs <- 1.2
colObs <- "darkred"
pointObsCol <- colObs

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
hicds <- args[1]
exprds <- args[2]
settingF <- file.path(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
cat(paste0(settingF, "\n"))
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

pipOutFold <- file.path(pipScriptDir, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds)

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8cOnlyFCC_runAllDown"
script_name <- "fcc_wave_plot.R"
#stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

#source("main_settings.R")
source(settingF)
#source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar


registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

outFolder <- file.path("FCC_WAVE_PLOT")
dir.create(outFolder, recursive=TRUE)

pipLogFile <- ""

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
cat(txt)

stopifnot(!any(duplicated(names(geneList))))

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
geneNbr <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))

### SET OUTPUT

# if plotSeparated == TRUE => 1 plot per ratio, otherwise all on the same figure (#x2 plots)
plotSeparated <- F

# "permThresh" the quantile of permutations to take is loaded from main_settings.R

###############
##### retrieve the direction of up/down
###############
# retrieve the direction of up/down
DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
exprDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
# samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
# samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
samp1 <- eval(parse(text=load(paste0(sample1_file))))
samp2 <- eval(parse(text=load(paste0(sample2_file))))
DE_topTable <- DE_topTable[DE_topTable$genes %in% names(DE_geneList),]
stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))
stopifnot(all(colnames(exprDT) %in% c(samp1, samp2)))
stopifnot(all(samp1 %in% colnames(exprDT)))
stopifnot(all(samp2 %in% colnames(exprDT)))
maxDownGene <- DE_topTable$genes[which.min(DE_topTable$logFC)]
stopifnot(maxDownGene %in% rownames(exprDT))
mean_expr1 <- mean(unlist(c(exprDT[maxDownGene, samp1])), na.rm=T)
mean_expr2 <- mean(unlist(c(exprDT[maxDownGene, samp2])), na.rm=T)

if(mean_expr1 > mean_expr2) {
  subtitDir <- paste0("down: ", toupper(cond1), " > ", toupper(cond2))
} else{
  subtitDir <- paste0("down: ", toupper(cond2), " > ", toupper(cond1))
}



################################****************************************************************************************
####################################################### ITERATE OVER RATIOS TO PLOT
################################****************************************************************************************

curr_ratio_type <- "prodSignedRatio"

#for(curr_ratio_type in allDown_limited) {
cat(paste0("*** START ", curr_ratio_type, "\n"))

obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
permut_currDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))

# ensure I used the same set of TADs for the permutation and for the calculation
# (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))

interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )

############################################################
############################################################ # filter the TADs and sort
############################################################
filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)

filter_permut_currDown_unsort <- permut_currDown[interReg,]
stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))

filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
rownames(filter_permut_currDown) <- NULL
stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))

# FOR ratioDown => plot ratioConcord, departure from 0.5

my_stat_curr_ratio <- "prodSignedRatio"
# => raw (departure 0)
# prodSignedRatio -> does not need to be transformed
filter_obs_curr_down_half <- filter_obs_curr_down
filter_permut_currDown_half <- filter_permut_currDown

# PLOT THE 1ST PLOT

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_genomeWide_FCC_cumsum_obs_permut.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

par(bty="l")

my_main <- paste0("Genome-wide intra-TAD fold-change concordance")
my_sub <- paste0(hicds, " - " , exprds)
my_xlab <- paste0("TADs ranked by FCC")
my_ylab <- paste0("FCC cumulative sum")

observ_vect= filter_obs_curr_down_half
permut_DT=filter_permut_currDown_half
my_stat = my_stat_curr_ratio


observ_vect <- sort(observ_vect, decreasing = T)
permut_DT <- apply(permut_DT, 2, sort, decreasing=T)

x_val <- c(1:length(observ_vect))
diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-0)))

meanPermut_cumsum <- apply(diff_05_permut, 1, mean)
qt95Permut_cumsum <- apply(diff_05_permut, 1, quantile, probs=0.95)

obs_cumsum <- cumsum(abs(observ_vect - 0))

plot(obs_cumsum ~ x_val,
   main= my_main,
   type = "l",
   pch = 16, cex = 0.7,
   xlab= my_xlab, 
   ylab= my_ylab,
   cex.main = plotCex,
   cex.lab = plotCex,
   cex.axis = plotCex,
   col = colObs,
   axes=FALSE,
   lwd = lwdObs,
   bty="l")
box()
axis(2, cex.axis=plotCex, cex.lab=plotCex)
axis(1, cex.axis=plotCex, cex.lab=plotCex, at = seq(from=0, to=2000, by=200))

mtext(my_sub, font=3)

polygon(x = c(x_val, rev(x_val)), 
      y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
      border=NA,
      col = polygonPermutCol)
# lines(
#   x = x_val,
#   y= meanPermut_cumsum,
#   col = meanPermutCol
# )

lines(
  x = x_val,
  y= qt95Permut_cumsum,
  col = qt95PermutCol
)

auc_obs <- auc(x = x_val, y = obs_cumsum)
auc_permut <- auc(x = x_val, y = meanPermut_cumsum)
auc_permutQt <- auc(x = x_val, y = qt95Permut_cumsum)


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "auc_obs.Rdata"))
save(auc_obs, file= outFile)
cat(paste0("... written: ", outFile, "\n"))
outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "auc_permutQt.Rdata"))
save(auc_permutQt, file= outFile)
cat(paste0("... written: ", outFile, "\n"))

# save(auc_obs, file= "auc_obs.Rdata")
# save(auc_permut, file="auc_permut.Rdata")
# save(x_val,file= "x_val.Rdata")
# save(obs_cumsum, file="obs_cumsum.Rdata")
# save(meanPermut_cumsum, file ="meanPermut_cumsum.Rdata")

pct_inc <- round(auc_obs/auc_permut,2)
pct_inc_qt <- round(auc_obs/auc_permutQt,2)

# 
# legend("topleft",
#      xjust=0.5, yjust=0,
#      pch = c(16, 15), 
#      legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
#      pt.cex = c(0.7, 2),
#      col = c(pointObsCol, polygonPermutCol),
#      bty="n")

# legend("topleft",
#        xjust=0.5, yjust=0,
#        pch = c(16, 15, -1, -1 ),
#        lty=c(-1,-1,1,1),
#        legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut.", "avg. permut.", "0.95-qt permut."), 
#        pt.cex = c(0.7, 2),
#        col = c(pointObsCol, polygonPermutCol, meanPermutCol, qt95PermutCol),
#        bty="n")
# legtxt1 <- as.expression(bquote(frac(AUC[obs.], AUC[avg.permut.]) == .(pct_inc)))
# legtxt2 <- as.expression(bquote(frac(AUC[obs.], AUC[qt.permut.]) == .(pct_inc_qt)))
# legend("bottomright", legend=c(legtxt1, legtxt2), bty="n")
legend("topleft",
       xjust=0.5, yjust=0,
       pch = c(16, 15, -1 ),
       lty=c(-1,-1,1),
       legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut.", "0.95-qt permut."), 
       pt.cex = c(0.7, 2),
       col = c(pointObsCol, polygonPermutCol, qt95PermutCol),
       bty="n")
legtxt <- as.expression(bquote(frac(AUC[obs.], AUC[permut.]) == .(pct_inc_qt)))
# legend("bottomright", legend=c(legtxt), bty="n")
legend("right", legend=c(legtxt), bty="n")

# text(x = length(observ_vect), y = 0+50,
#      pos=2,
#      bquote(frac(AUC[obs.], AUC[permut.]) == .(pct_inc)))


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
#printAndLog(txt, pipLogFile)
cat(txt)

cat(paste0("*** DONE: ", script_name, "\n"))


