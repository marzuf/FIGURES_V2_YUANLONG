hicds <- "LG1_40kb"
exprds <- "TCGAluad_norm_luad"

plotType <- "png"

source("settings.R")
require(ggsci)
require(VennDiagram)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript geneRank_tadRank.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc
# Rscript geneRank_tadRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad

outFolder <- file.path("GENERANK_TADRANK")
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
hicds <- args[1]
exprds <- args[2]

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
# yarrr::transparent("grey", trans.val = .9)
# [1]  "#BEBEBE19"
mid_col <-  "#BEBEBE19"

geneSignifPval <- 0.05
tadSignifPval <- 0.01


nTop <- 10

top_reg <- unique(ds_dt$region[ds_dt$tad_rank <= nTop])
last_reg <- unique(ds_dt$region[ds_dt$rev_tad_rank <= nTop])


ds_dt$dotCol <- ifelse(ds_dt$region %in% top_reg, top_col,
                         ifelse(ds_dt$region %in% last_reg, last_col, mid_col))
stopifnot(!is.na(ds_dt$dotCol))
plot_dt <- ds_dt[ds_dt$region %in% top_reg | ds_dt$region %in% last_reg ,]


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_geneRank_tadRank_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l", family=fontFamily)
densplot(
  y = ds_dt$gene_rank,
  x = ds_dt$tad_rank,
  ylab = "Gene rank",
  xlab = "TAD rank",
  main =  paste0(hicds, " - ", exprds),
  cex = 0.7,
  cex.lab=plotCex,
  cex.axis=plotCex,
  cex.main=1
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
par(bty="l", family=fontFamily)
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


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_geneRank_tadRank_geneOrder_top", nTop, "_barplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
par(bty="l", family=fontFamily)
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



#########################################################################################################################
#### BARPLOTS - all gene ranks
#########################################################################################################################
maxGeneRank <- max(ds_dt$gene_rank)
maxTADrank <- max(ds_dt$tad_rank)

ds_dt$gene_rank_rel <- ds_dt$gene_rank/maxGeneRank
ds_dt$tad_rank_rel <- ds_dt$tad_rank/maxTADrank


geneBar_pos <- 1
tadBar_pos <- 2
axisOffset <- 0.5

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_geneRank_tadRank_topBars_nTop", nTop, "_segments.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(bty="n", family=fontFamily)
plot(NULL,
     xlab="",
     ylab="",
     axes=F,
     main=paste0("Top and last ", nTop, " TADs"),
     # sub=paste0(),
     cex.axis=plotCex,
     cex.lab = plotCex,
     cex.main = plotCex,
  xlim = c(geneBar_pos-axisOffset, tadBar_pos+axisOffset),
  # ylim = c(-axisOffset, 1+axisOffset)
  ylim = c(1+axisOffset, -axisOffset)  # reverse to have the top ones at the top !
)
# mtext(side=3, paste0(hicds, " - ", exprds), cex=plotCex)
mtext(side=3, paste0(hicds, " - ", exprds), cex=plotCex-0.2)

segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=ds_dt$gene_rank_rel[ds_dt$dotCol %in% mid_col],
         y1=ds_dt$tad_rank_rel[ds_dt$dotCol %in% mid_col],
         lty=3,
         tcl=-.2,
         col = mid_col
)
segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=ds_dt$gene_rank_rel[ds_dt$dotCol %in% top_col],
         y1=ds_dt$tad_rank_rel[ds_dt$dotCol %in% top_col],
         col = top_col
         )
segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=ds_dt$gene_rank_rel[ds_dt$dotCol %in% last_col],
         y1=ds_dt$tad_rank_rel[ds_dt$dotCol %in% last_col],
         col = last_col
)


# add bar for the genes
segments(x0=geneBar_pos, y0=0, x1=geneBar_pos, y1=1, lwd=5)

# add bar for the TADs
segments(x0=tadBar_pos, y0=0, x1=tadBar_pos, y1=1, lwd=5)

text(
  x = c(geneBar_pos, tadBar_pos),
  # y = 1 + 0.2,
  y = 0 - 0.2, 
  labels=c("Gene ranks", "TAD ranks"),
  cex = plotCex
)
legend(
  "bottom",
  # cex=0.6,
  # horiz=T,
  lty=c(1, -1, 1, -1),
  col = c(top_col,top_col, last_col, last_col),
  legend=c(paste0("# top TADs=", length(top_reg)),
           paste0("# genes Top TADs=", sum(ds_dt$dotCol %in% top_col)),
           paste0("# last TADs=", length(last_reg)),
           paste0("# genes Last TADs=", sum(ds_dt$dotCol %in% last_col))),
  bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#########################################################################################################################
#### BARPLOTS - all gene med ranks
#########################################################################################################################

med_ds_dt <- aggregate(gene_rank ~ hicds + exprds+ region, data=ds_dt, FUN=median)
colnames(med_ds_dt)[colnames(med_ds_dt) == "gene_rank"] <- "med_gene_rank"

med_plot_dt <- merge(med_ds_dt, tad_dt, all=TRUE, by=c("hicds", "exprds","region"))
stopifnot(!duplicated(med_plot_dt$region))


maxGeneRank <- max(med_plot_dt$med_gene_rank)
maxTADrank <- max(med_plot_dt$tad_rank)

med_plot_dt$gene_rank_rel <- med_plot_dt$med_gene_rank/maxGeneRank
med_plot_dt$tad_rank_rel <- med_plot_dt$tad_rank/maxTADrank

med_plot_dt$dotCol <- ifelse(med_plot_dt$region %in% top_reg, top_col,
                       ifelse(med_plot_dt$region %in% last_reg, last_col, mid_col))
stopifnot(!is.na(med_plot_dt$dotCol))

geneBar_pos <- 1
tadBar_pos <- 2
axisOffset <- 0.5

outFile <- file.path(outFolder, paste0(hicds,"_", exprds, "_medGeneRank_tadRank_topBars_nTop", nTop, "_segments.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(bty="n", family=fontFamily)
plot(NULL,
     xlab="",
     ylab="",
     axes=F,
     main=paste0("Top and last ", nTop, " TADs"),
     # sub=paste0(),
     cex.axis=plotCex,
     cex.lab = plotCex,
     cex.main = plotCex,
     xlim = c(geneBar_pos-axisOffset, tadBar_pos+axisOffset),
     # ylim = c(-axisOffset, 1+axisOffset)
     ylim = c(1+axisOffset, -axisOffset)
)
# mtext(side=3, paste0(hicds, " - ", exprds), cex=plotCex)
mtext(side=3, paste0(hicds, " - ", exprds), cex=plotCex-0.2)

segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=med_plot_dt$gene_rank_rel[med_plot_dt$dotCol %in% mid_col],
         y1=med_plot_dt$tad_rank_rel[med_plot_dt$dotCol %in% mid_col],
         lty=3,
         tcl=-.2,
         col = mid_col
)
segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=med_plot_dt$gene_rank_rel[med_plot_dt$dotCol %in% top_col],
         y1=med_plot_dt$tad_rank_rel[med_plot_dt$dotCol %in% top_col],
         col = top_col
)
segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=med_plot_dt$gene_rank_rel[med_plot_dt$dotCol %in% last_col],
         y1=med_plot_dt$tad_rank_rel[med_plot_dt$dotCol %in% last_col],
         col = last_col
)


# add bar for the genes
segments(x0=geneBar_pos, y0=0, x1=geneBar_pos, y1=1, lwd=5)

# add bar for the TADs
segments(x0=tadBar_pos, y0=0, x1=tadBar_pos, y1=1, lwd=5)

text(
  x = c(geneBar_pos, tadBar_pos),
  # y = 1 + 0.2,
  y = 0 - 0.2,
  labels=c("Median gene rank", "TAD ranks"),
  cex = plotCex
)
legend(
  "bottom",
  # cex=0.6,
  # horiz=T,
  lty=c(1,1),
  col = c(top_col, last_col),
  legend=c(paste0("# top TADs=", length(top_reg)),
           # paste0("# genes Top TADs=", sum(med_plot_dt$dotCol %in% top_col)),
           # paste0("# genes Last TADs=", sum(med_plot_dt$dotCol %in% last_col)),
          paste0("# last TADs=", length(last_reg))),
  bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




#########################################################################################################################
#### VENN DIAGRAM
#########################################################################################################################
stopifnot(!duplicated(ds_dt$entrezID))
stopifnot(!duplicated(tad_dt$region))

signif_genes <-  ds_dt$entrezID[ds_dt$adj.P.Val <= geneSignifPval ]
signif_tads <-  ds_dt$entrezID[ds_dt$tad_adjCombPval <= tadSignifPval ]

signif_genesOnly <- ds_dt$entrezID[ds_dt$adj.P.Val <= geneSignifPval &
                                     ds_dt$tad_adjCombPval > tadSignifPval ]
nSignif_genesOnly <- length(signif_genesOnly)

signif_genesAndTADs <- ds_dt$entrezID[ds_dt$adj.P.Val <= geneSignifPval &
                                     ds_dt$tad_adjCombPval <= tadSignifPval ]
nSignif_genesAndTADs <- length(signif_genesAndTADs)

signif_tadsOnly <- ds_dt$entrezID[ds_dt$adj.P.Val > geneSignifPval &
                                        ds_dt$tad_adjCombPval <= tadSignifPval ]
nSignif_tadsOnly <- length(signif_tadsOnly)

stopifnot(nSignif_genesOnly+nSignif_tadsOnly+nSignif_genesAndTADs == sum(ds_dt$tad_adjCombPval <= tadSignifPval | ds_dt$adj.P.Val <= geneSignifPval))

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_nbrSignifGenes_vennDiagram.", plotType))

require(ggplot2)
require(ggsci)

vd <- venn.diagram(
  x = list(signif_genes, signif_tads),
  main = paste0("# signif. genes"),
  sub = paste0(hicds, "\n", exprds),
  category.names = c(paste0("signif. gene-level\n(", nSignif_genesOnly+nSignif_genesAndTADs, ")") , paste0("signif. TAD level\n(", nSignif_tadsOnly+nSignif_genesAndTADs, ")")),
  
  fontfamily="Hershey",
  cat.fontfamily="Hershey",
  main.fontfamily="Hershey",
  sub.fontfamily="Hershey",
  
  sub.fontface="plain",
  main.fontface="bold",
  
  margin=c(0,0,0,0),
  
  main.cex = 2,
  sub.cex=1.6,
  cat.cex=1.4,
  cex = 2,
  
  cat.default.pos="outer",
  cat.pos=c(-25, 25),
  
  cat.col = pal_lancet()(2),
  fill = pal_lancet()(2),
  col = pal_lancet()(2),
  alpha=0.2,
  scaled=FALSE,
  filename = NULL
)
ggsave(vd, file=outFile,height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

system(paste0("rm -f VennDiagram2019*.log"))

