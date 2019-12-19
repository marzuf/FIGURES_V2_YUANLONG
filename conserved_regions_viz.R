# IGV style


# Rscript conserved_regions_viz.R
# Rscript conserved_regions_viz.R norm_vs_tumor
# Rscript conserved_regions_viz.R subtypes
# Rscript conserved_regions_viz.R wt_vs_mut

# Rscript conserved_regions_viz.R <cmpType>

cat("> START ", "conserved_regions_viz.R", "\n")

startTime <- Sys.time()

plotType <- "svg"

source("settings.R")

require(ggsci)
tad_col <- pal_d3()(2)[1]
gene_col <- pal_d3()(2)[2]

outFolder <- file.path("CONSERVED_REGIONS_VIZ")
dir.create(outFolder, recursive = TRUE)


args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  cmpType <- ""
  filePrefix <- ""
  cmpTit <- paste0("all")
} else if(length(args) == 1) {
  cmpType <- args[1]  
  filePrefix <- paste0(cmpType, "_")
  cmpTit <- cmpType
}else {
  stop("---error\n")
}
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3
gene_matching_fuse_threshold <- 0.8


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)


inFolder <- file.path(runFolder, "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2")
outFile <- file.path(inFolder, paste0(filePrefix, "conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_dt <- get(load(outFile))
conserved_dt$conserved_region <- as.character(conserved_dt$conserved_region)

outFile <- file.path(inFolder, paste0(filePrefix, "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_list <- get(load(outFile))
nConserved <- lengths(conserved_list)

stopifnot(length(nConserved) == nrow(conserved_dt))

maxConserved <- names(which.max(nConserved))


max_dt <- conserved_dt[conserved_dt$conserved_region == maxConserved,,drop=F]
stopifnot(nrow(max_dt) == 1)

all_max_regions <- unlist(strsplit(max_dt$corresp_tads, split=","))
stopifnot(length(all_max_regions) == nConserved[maxConserved])

all_max_entrez <- unlist(strsplit(max_dt$intersect_genes_entrez, split=","))
stopifnot(all_max_entrez %in% gff_dt$entrezID)
all_genes_starts_ends <- sapply(all_max_entrez, function(x) {
  c(start = gff_dt$start[gff_dt$entrezID == x],
    end = gff_dt$end[gff_dt$entrezID == x],
    symbol = gff_dt$symbol[gff_dt$entrezID == x]
    )
})


all_regions_starts_ends <- sapply(all_max_regions, function(x) {
    g2t_dt <- read.delim(file.path(runFolder, dirname(dirname(x)), "genes2tad", "all_assigned_regions.txt"), header=FALSE, stringsAsFactors=FALSE, col.names=c("chromo", "region", "start", "end"))
    stopifnot(sum(g2t_dt$region == basename(x)) == 1)
    g2t_dt$start[g2t_dt$region == basename(x)]
    c(start=g2t_dt$start[g2t_dt$region == basename(x)], end=g2t_dt$end[g2t_dt$region == basename(x)], chromo = g2t_dt$chromo[g2t_dt$region == basename(x)])
  })

all_exprds <- basename(dirname(colnames(all_regions_starts_ends)))
all_regions_starts_ends <- all_regions_starts_ends[,rev(colnames(all_regions_starts_ends)[order(all_exprds)])]

chromo <- unique(as.character(all_regions_starts_ends["chromo",]))
stopifnot(length(chromo) == 1)


dsSpace <- 0.5
geneSpace <- 0.2

tadOffset <- 50000
yOffset <- 0.3

textOffset <- 0.05

xStart <- min(as.numeric(as.character(all_regions_starts_ends["start",]))) - tadOffset
xEnd <- max(as.numeric(as.character(all_regions_starts_ends["end",]))) + tadOffset

dsPos <- seq(from=0, by=dsSpace, length.out=ncol(all_regions_starts_ends))
genePos <- seq(from = max(dsPos) + dsSpace, by=geneSpace, length.out=ncol(all_genes_starts_ends))
yStart <- min(c(dsPos, genePos)) - yOffset
yEnd <- max(c(dsPos, genePos)) + yOffset

outFile <- file.path(outFolder, paste0(maxConserved, "_viz.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
initMar <- par()$mar
par(mar=initMar+c(0,10,0,0))
par(family="Heshley")
par(xpd=TRUE)
plot(NULL,
     main = paste0(maxConserved),
     xlim = c(xStart, xEnd),
     ylim = c(yStart, yEnd),
     # xlab = "",
     xlab = paste0(gsub("chr", "chromosome ", chromo)),
     ylab = "",
     axes = FALSE,
     cex.main = plotCex
     )
axis(1,
     at = unique(sort(c(as.numeric(as.character(all_regions_starts_ends["start",])), as.numeric(as.character(all_regions_starts_ends["end",]))))), 
     cex = 0.6)
# draw the tads
segments(
  x0 = as.numeric(all_regions_starts_ends["start",]),
  y0 = dsPos,
  x1 = as.numeric(all_regions_starts_ends["end",]),
  y1 = dsPos,
  col=tad_col
)
text(
  x = as.numeric(all_regions_starts_ends["start",]),
  y = dsPos + textOffset,
  # labels = colnames(all_regions_starts_ends),
  labels = dirname(colnames(all_regions_starts_ends)),
  # cex = 0.5,
  cex = 0.7,
  pos=2,
  col = tad_col
)
segments(
  x0 = as.numeric(all_genes_starts_ends["start",]),
  y0 = genePos,
  x1 = as.numeric(all_genes_starts_ends["end",]),
  y1 = genePos,
  col=gene_col
)
text(
  x = 0.5*(as.numeric(all_genes_starts_ends["start",]) + as.numeric(all_genes_starts_ends["end",])),
  y = genePos - 5*textOffset,
  # y = genePos + 0,
  labels = as.character(all_genes_starts_ends["symbol",]),
  cex = 1,
  pos=3,
  col=gene_col
)

segments(x0=c(as.numeric(all_genes_starts_ends["start",]), as.numeric(all_genes_starts_ends["end",])),
         y0=0-yOffset,
         x1=c(as.numeric(all_genes_starts_ends["start",]), as.numeric(all_genes_starts_ends["end",])),
         y1 = rep(genePos,each=2),
         lty=2, col = gene_col)

cat(paste0("... written: ", outFile, "\n"))



























