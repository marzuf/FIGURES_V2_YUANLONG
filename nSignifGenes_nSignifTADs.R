
plotType <- "png"


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("settings.R")
require(ggsci)

# Rscript nSignifGenes_nSignifTADs.R

outFolder <- file.path("NSIGNIFGENES_NSIGNIFTADS")
dir.create(outFolder, recursive = TRUE)

geneSignifThresh <- 0.05
tadSignifThresh <- 0.01

gene_col <- pal_d3()(2)[1]
tad_col <- pal_d3()(2)[2]

dotSymb <- "\u25CF"


inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
geneDT <- inDT[,c("hicds", "exprds", "entrezID", "adj.P.Val")]
geneDT <- unique(geneDT)
nSignifGenes_dt <- aggregate(adj.P.Val~hicds + exprds, data = geneDT, FUN=function(x) sum(x<=geneSignifThresh))
colnames(nSignifGenes_dt)[colnames(nSignifGenes_dt) == "adj.P.Val"] <- "nSignifGenes"

tadDT <- inDT[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tadDT <- unique(tadDT)
nSignifTADs_dt <- aggregate(tad_adjCombPval~hicds + exprds, data = tadDT, FUN=function(x) sum(x<=tadSignifThresh))
colnames(nSignifTADs_dt)[colnames(nSignifTADs_dt) == "tad_adjCombPval"] <- "nSignifTADs"

nSignif_dt <- merge(nSignifGenes_dt, nSignifTADs_dt, by=c("hicds", "exprds"), all=TRUE)
stopifnot(!is.na(nSignif_dt))

nSignif_dt <- nSignif_dt[order(nSignif_dt$nSignifTADs, decreasing = TRUE),]
nSignif_dt$dataset <- paste0(nSignif_dt$hicds, "\n", nSignif_dt$exprds)

labcols <- all_cols[all_cmps[nSignif_dt$exprds]]

maxTADs <- max(ceiling(nSignif_dt$nSignifTADs/10)*10)
maxGenes <- max(ceiling(nSignif_dt$nSignifGenes/1000)*1000)

outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withSymb.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2))
par(bty="U", family=fontFamily)
par(mar = c(5,5,2,5))
plot(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifTADs, 
  col = tad_col,
  ylim=c(0, maxTADs),
  ylab = "# signif. TADs",
  xlab = "",
  pch=3,
  main = paste0("# signif. features"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  col.lab = tad_col,
  axes = FALSE
)


mtext(side=3, line=-1, text=paste0("all datasets - n= ", length(unique(nSignif_dt$dataset))))
# mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex =0.5)
axis(2, col = tad_col, col.ticks = tad_col, col.axis=tad_col, at=seq(from=0, to = maxTADs, by=10))
axis(1, labels=F, lwd.ticks = -1)
mtext(side=1, col = labcols, text = rep(dotSymb, nrow(nSignif_dt)), at= 1:nrow(nSignif_dt), las=2, cex = 1.2)
par(new = T, family=fontFamily)
plot(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifGenes,
  col = gene_col,
  ylab = NA,
  ylim=c(0, maxGenes),
  xlab = NA,
  pch=3,
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  axes = FALSE
)
axis(side=4, col = gene_col, col.ticks = gene_col, col.axis=gene_col, at = seq(from=0, to=maxGenes, by=1000))
mtext(side = 4, line = 3, '# signif. genes', col=gene_col,  cex=plotCex)

# greycol <-  "#BEBEBE19"
greycol <- "grey"
abline(v=1:nrow(nSignif_dt), lty=3, col = greycol)

# mtext(side = 1, line = -2, 
#       text = paste0(dotSymb, " ", names(all_cols)),
#       col=all_cols, cex=plotCex)

legend("bottom", 
       # legend = paste0(dotSymb, " ", names(all_cols)),
       legend = paste0(names(all_cols)),
       col=all_cols,
       pch=16,
       # lty=c(1,2),
       horiz=TRUE,
       inset=c(0,-0.2),
       cex = plotCex,
       xpd=TRUE, bty="n"
)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withLeg.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2.5))
par(bty="U")
par(mar = c(5+2,5,2,5))
plot(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifTADs, 
  col = tad_col,
  ylim=c(0, maxTADs),
  ylab = "# signif. TADs",
  xlab = "",
  pch=3,
  main = paste0("# signif. features"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  col.lab = tad_col,
  axes = FALSE
)
mtext(side=3, line=-1, text=paste0("all datasets - n= ", length(unique(nSignif_dt$dataset))))
# mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex =0.5)
axis(2, col = tad_col, col.ticks = tad_col, col.axis=tad_col, at=seq(from=0, to = maxTADs, by=10))
axis(1, labels=F, lwd.ticks = -1)
mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex = 0.6)
par(new = T)
plot(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifGenes,
  col = gene_col,
  ylab = NA,
  ylim=c(0, maxGenes),
  xlab = NA,
  pch=3,
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  axes = FALSE
)
axis(side=4, col = gene_col, col.ticks = gene_col, col.axis=gene_col, at = seq(from=0, to=maxGenes, by=1000))
mtext(side = 4, line = 3, '# signif. genes', col=gene_col, cex=plotCex)

# greycol <-  "#BEBEBE19"
greycol <- "grey"
abline(v=1:nrow(nSignif_dt), lty=3, col = greycol)


mtext(side = 1, line = -2, 
      text = paste0(dotSymb, " ", names(all_cols)),
      col=all_cols, cex=plotCex)


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))






