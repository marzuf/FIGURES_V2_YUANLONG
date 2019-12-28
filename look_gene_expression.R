########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript look_TAD_expression.R\n"))

script_name <- "look_TAD_expression.R"


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

require(reshape2)
log10_offset <- 0.01


# Rscript look_TAD_expression.R <hicds> <exprds> <gene_symbol>
# Rscript look_TAD_expression.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390


hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAlusc_norm_lusc"
tad_to_plot="chr11_TAD390"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
tad_to_plot <- args[3]


# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 7
myWidthGG <- 7

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

mainFolder <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(mainFolder, "PIPELINE", "INPUT_FILES")

outFolder <- file.path("LOOK_TAD_EXPRESSION")
dir.create(outFolder, recursive = TRUE)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end, gff_dt$entrezID),]

g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)

stopifnot(tad_to_plot %in% g2t_dt$region)

settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

samp1 <- get(load(file.path(setDir, sample1_file)))
samp2 <- get(load(file.path(setDir, sample2_file)))


geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
stopifnot(geneList %in% g2t_dt$entrezID)

toplot_dt <- g2t_dt[g2t_dt$region == tad_to_plot & g2t_dt$entrezID %in% geneList,]
toplot_dt <- toplot_dt[order(toplot_dt$start, toplot_dt$end),]

stopifnot(toplot_dt$entrezID %in% geneList)
plot_geneList <- geneList[geneList %in% toplot_dt$entrezID]


count_file <- file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_rnaseqDT.Rdata")
stopifnot(file.exists(count_file))
fpkm_dt <- get(load(count_file))
stopifnot(names(geneList) %in% rownames(fpkm_dt))

fpkm_plot_dt <- fpkm_dt[rownames(fpkm_dt) %in% names(plot_geneList),]
fpkm_plot_dt$entrezID <- plot_geneList[paste0(rownames(fpkm_plot_dt))]
stopifnot(!is.na(fpkm_plot_dt$entrez))

stopifnot(samp1 %in% colnames(fpkm_plot_dt))
stopifnot(samp2 %in% colnames(fpkm_plot_dt))

fpkm_plot_dt <- fpkm_plot_dt[,c("entrezID", samp1, samp2)]
m_fpkm_dt <- melt(fpkm_plot_dt, id="entrezID")

stopifnot(m_fpkm_dt$entrezID %in% toplot_dt$entrezID)

toplot_dt <- merge(m_fpkm_dt, toplot_dt, merge="entrezID", all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(toplot_dt))
toplot_dt$cond <- ifelse(toplot_dt$variable %in% samp1, cond1, 
                         ifelse(toplot_dt$variable %in% samp2, cond2, NA ))
stopifnot(!is.na(toplot_dt$cond))
toplot_dt$symbol <- entrez2symb[paste0(toplot_dt$entrezID)]
stopifnot(!is.na(toplot_dt$symbol))



toplot_dt <- toplot_dt[order(toplot_dt$chromo, toplot_dt$start, toplot_dt$end, toplot_dt$value),]
toplot_dt$value_log10 <- log10(toplot_dt$value + log10_offset)

withRank_toplot_dt <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol, toplot_dt$cond), function(x) {
  dt <- x[order(x$value, decreasing = TRUE),]
  dt$samp_rank <- 1:nrow(dt)
  dt
}))

require(ggsci)
ggsci_pal <- "d3"
ggsci_subpal <- ""
# p_var <-  ggplot(withRank_toplot_dt, aes(x = samp_rank, y = value_log10, fill = cond)) + 
#   geom_point()+
#   # geom_boxplot()+
#   facet_grid(~symbol+cond, switch="x") + 
#   coord_cartesian(expand = FALSE) +
#   ggtitle(paste0(""), subtitle = paste0(""))+
#   scale_x_discrete(name="")+
#   scale_y_continuous(name=paste0(""),
#                      breaks = scales::pretty_breaks(n = 20))+
#   
#   eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
#   eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
# 
#   labs(fill  = paste0("")) +
#   theme( 
#     strip.text = element_text(size = 12),
#     plot.title = element_text(hjust = 0.5, face = "bold", size=16),
#     plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
#     panel.grid = element_blank(),
#     panel.grid.major.y = element_line(colour = "grey"),
#     panel.grid.minor.y = element_line(colour = "grey"),
#     strip.text.x = element_text(size = 10),
#     axis.line.x = element_line(size = .2, color = "black"),
#     axis.line.y = element_line(size = .3, color = "black"),
#     axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.y = element_text(color="black", size=12),
#     axis.title.x = element_text(color="black", size=12),
#     panel.border = element_blank(),
#     panel.background = element_rect(fill = "transparent"),
#     legend.background =  element_rect(),
#     legend.key = element_blank(),
#     legend.title = element_text(face="bold")
#   )
# 
# outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_dotplot.", plotType))
# ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
# cat(paste0("... written: ", outFile, "\n"))


withRank_toplot_dt2 <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol), function(x) {
  x$cond <- factor(x$cond, levels=c(cond1,cond2))
  dt <- x[order(-as.numeric(x$cond), x$value, decreasing = TRUE),]
  dt$samp_rank <- 1:nrow(dt)
  dt
}))


tmp <- withRank_toplot_dt2[,c("symbol", "start", "end")]
tmp <- unique(tmp)
tmp <- tmp[order(tmp$start, tmp$end),]

withRank_toplot_dt2$symbol <- factor(withRank_toplot_dt2$symbol, levels=tmp$symbol)
withRank_toplot_dt2$cond <- factor(withRank_toplot_dt2$cond, levels = c(cond1,cond2))

p_var <-  ggplot(withRank_toplot_dt2, aes(x = samp_rank, y = value_log10, fill = cond, color=cond)) + 
  geom_point()+
  # geom_boxplot()+
  facet_grid(~symbol, switch="x") + 
  coord_cartesian(expand = FALSE) +
  ggtitle(paste0(""), subtitle = paste0(""))+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(""),
                     breaks = scales::pretty_breaks(n = 20))+
  
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  
  labs(fill  = paste0(""), color=paste0("")) +
  theme( 
    strip.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_dotplot.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))








# 
# 
# 
# 
# 
# 
# 
# 
# 
# all_gene_entrez_0 <- all_gene_entrez
# 
# all_gene_entrez <- all_gene_entrez[all_gene_entrez %in% rownames(fpkm_dt)]
# if(length(all_gene_entrez) < length(all_gene_entrez_0)) {
#   missing_entrez <- all_gene_entrez_0[!all_gene_entrez_0 %in% all_gene_entrez]
#   cat(paste0("!!! WARNING: following genes discarded (missing von fpkm_dt):\t", paste0(entrez2symb[paste0(missing_entrez)], collapse=","), "\n"))
#   all_gene_symbols <- entrez2symb[paste0(all_gene_entrez)]
# }
# 
# stopifnot(all_gene_entrez %in% rownames(fpkm_dt))
# stopifnot(samp1 %in% colnames(fpkm_dt))
# stopifnot(samp2 %in% colnames(fpkm_dt))
# 
# plotTit <- paste0(hicds, " - ", exprds)
# plotSubTit <- paste0("~ ", paste0(all_gene_symbols, collapse=","), " ~")
# myylab <- "RNA-seq scaled estimates"
# 
# if(length(all_gene_symbols) == 0) {
#   stop(cat(paste0("... no data found for this gene\n")))
# } else if(length(all_gene_symbols) == 1) {
#   
#   
#   cond1_fpkm <- unlist(fpkm_dt[paste0(all_gene_entrez), samp1])
#   cond2_fpkm <- unlist(fpkm_dt[paste0(all_gene_entrez), samp2])
#   
#   boxplot_dt <- data.frame(
#     FPKM = c(cond1_fpkm, cond2_fpkm),
#     cond = c(rep(cond1, length(samp1)), rep(cond2, length(samp2))),
#     stringsAsFactors = FALSE
#   )
#   
#   outFile <- file.path(outFolder, paste0(all_gene_symbols, "_", hicds, "_", exprds, ".", plotType))
#   do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
#   par(bty="l")
#   boxplot(FPKM~cond, data=boxplot_dt,
#           main = plotTit,
#           ylab=myylab,
#           xlab="",
#           cex.axis = plotCex,
#           cex.lab=plotCex)
#   mtext(side=3, text = plotSubTit, cex=1.4)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))    
#   
# } else {
# 
#   count_dt <- fpkm_dt[paste0(all_gene_entrez), c(samp1, samp2)]
#   
#   count_dt_t <- as.data.frame(t(count_dt))
#   stopifnot(colnames(count_dt_t) %in% names( entrez2symb))
#   colnames(count_dt_t) <- entrez2symb[colnames(count_dt_t)]
#   
#   stopifnot(all_gene_symbols %in% colnames(count_dt_t))
#   
#   gene_order <- init_order[init_order %in% all_gene_symbols]
#   stopifnot(length(gene_order) == length(all_gene_symbols))
#   
#   count_dt_t$variable <- rownames(count_dt_t)
#   count_dt_t$cond <- ifelse(count_dt_t$variable %in% samp1, cond1, 
#          ifelse(count_dt_t$variable %in% samp2, cond2, NA))
#   
#   
#   expr_p <- ggboxplot(count_dt_t, x = "cond",
#             xlab = "",
#             scales='free_y',
#             y = gene_order,
#             combine = TRUE,
#             ylab = myylab,
#             color = "cond", palette = "jco") + 
#     guides(color=FALSE)+
#     # ggtitle(plotTit, sub=plotSubTit)+
#     ggtitle(plotTit)+
#     theme(
#       plot.title = element_text(size=16, hjust=0.5, face = "bold"),
#       plot.subtitle = element_text(size=14, hjust=0.5),
#       strip.text.x =  element_text(size=12),
#       axis.text.x = element_text(size=12),
#       axis.title.y = element_text(size=14)
#     )
#   
#   
#   outFile <- file.path(outFolder, paste0(paste0(all_gene_symbols, collapse="_"), "_", hicds, "_", exprds, ".", plotType))
#   ggsave(filename = outFile, plot = expr_p, height=myHeightGG, width=myWidthGG*length(all_gene_symbols)*0.8)
#   cat(paste0("... written: ", outFile, "\n"))    
#   
#   
# }
# 
# 
# 
# 

##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

