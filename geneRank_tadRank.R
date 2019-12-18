hicds <- "LG1_40kb"
exprds <- "TCGAluad_norm_luad"

plotType <- "png"

source("settings.R")
require(ggsci)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript geneRank_tadRank.R

outFolder <- file.path("GENERANK_TADRANK")
dir.create(outFolder, recursive = TRUE)

inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

ds_dt <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

tad_dt <- ds_dt[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tad_dt <- unique(tad_dt)
tad_dt$tad_rank <- rank(tad_dt$tad_adjCombPval, ties="min")
tad_dt$rev_tad_rank <- rank(-tad_dt$tad_adjCombPval, ties="min")

ds_dt <- merge(ds_dt, tad_dt[, c("hicds", "exprds", "region", "rev_tad_rank")], all=TRUE, by=c("hicds", "exprds", "region"))

nTADs <- length(unique(ds_dt$region))

head(ds_dt[order(ds_dt$tad_adjCombPval),])

top_col <- pal_d3()(2)[1]
last_col <- pal_d3()(2)[2]

nTop <- 10

top_reg <- unique(ds_dt$region[ds_dt$tad_rank <= nTop])
last_reg <- unique(ds_dt$region[ds_dt$rev_tad_rank <= nTop])


plot_dt <- ds_dt[ds_dt$region %in% top_reg | ds_dt$region %in% last_reg ,]

plot_dt$dotCol <- ifelse(plot_dt$region %in% top_reg, top_col,
                         ifelse(plot_dt$region %in% last_reg, last_col, NA))
stopifnot(!is.na(plot_dt$dotCol))

outFile <- file.path(outFolder, paste0("geneRank_tadRank_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l")
densplot(
  y = ds_dt$gene_rank,
  x = ds_dt$tad_rank,
  ylab = "Gene rank",
  xlab = "TAD rank",
  main =  paste0(hicds, " - ", exprds),
  cex = 0.7,
  cex.lab=plotCex,
  cex.axis=plotCex
)
mtext(side=3, text = paste0("# genes = ", nrow(ds_dt)))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# boxplot(gene_rank ~ tad_rank, data=ds_dt)
# ds_dt <- ds_dt[order(ds_dt$tad_rank, ds_dt$region, ds_dt$entrezID), ]
# barplot(ds_dt$gene_rank)#at= 1:length(ds_dt$tad_rank))
# plot(gene_rank ~ tad_rank, data=ds_dt)
# plot(
#   x = plot_dt$gene_rank,
#   y = plot_dt$tad_rank,
#   main=paste0(hicds, " - ", exprds),
#   cex.axis = plotCex,
#   cex.lab = plotCex,
#   cex = 0.7,
#   col = plot_dt$dotCol
# )

plot_dt <- plot_dt[order(plot_dt$tad_rank, plot_dt$region, plot_dt$entrezID), ]


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_geneRank_tadRank_tadOrder_top", nTop, "_barplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
par(bty="l")
barpos <- barplot(plot_dt$gene_rank, col = plot_dt$dotCol,
                  main = paste0(hicds, " - ", exprds),
                  cex.axis = plotCex,
                  cex.lab = plotCex,
        xlab = paste0("(", nrow(plot_dt), " genes ranked by TAD rank)"),
        ylab = "Gene rank")
topLegPos <- mean(barpos[which(plot_dt$dotCol == top_col)])
# mtext(side = 3, text = paste0("Genes from the top ", nTop, " TADs"), at = topLegPos, col = top_col)
mtext(side = 3, text = paste0("Genes from the top ", nTop, " TADs (", length(top_reg), ")"), at = topLegPos, col = top_col, cex=plotCex)
lastLegPos <- mean(barpos[which(plot_dt$dotCol == last_col)])
# mtext(side = 3, text = paste0("Genes from the last ", nTop, " TADs"), at = lastLegPos, col = last_col)
mtext(side = 3, text = paste0("Genes from the last ", nTop, " TADs (", length(last_reg), ")"), at = lastLegPos, col = last_col, cex=plotCex)
# mtext(side=3, text = paste0("# genes = ", nrow(plot_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


plot_dt <- plot_dt[order(plot_dt$dotCol, plot_dt$gene_rank), ]


outFile <- file.path(outFolder, paste0("geneRank_tadRank_geneOrder_top", nTop, "_barplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
par(bty="l")
barpos <- barplot(plot_dt$gene_rank, col = plot_dt$dotCol,
                  main = paste0(hicds, " - ", exprds),
                  cex.axis = plotCex,
                  cex.lab = plotCex,
               xlab="(Genes ranked by gene rank)",
               ylab = "Gene rank")
topLegPos <- mean(barpos[which(plot_dt$dotCol == top_col)])
# mtext(side = 3, text = paste0("Genes from the top ", nTop, " TADs"), at = topLegPos, col = top_col)
mtext(side = 3, text = paste0("Genes from the top ", nTop, " TADs (", length(top_reg), ")"), at = topLegPos, col = top_col, cex=plotCex)
lastLegPos <- mean(barpos[which(plot_dt$dotCol == last_col)])
# mtext(side = 3, text = paste0("Genes from the last ", nTop, " TADs"), at = lastLegPos, col = last_col)
mtext(side = 3, text = paste0("Genes from the last ", nTop, " TADs (", length(last_reg), ")"), at = lastLegPos, col = last_col, cex=plotCex)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



