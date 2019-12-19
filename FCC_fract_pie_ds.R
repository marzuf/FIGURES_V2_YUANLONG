hicds <- "LG1_40kb"
exprds <- "TCGAluad_norm_luad"

plotType <- "png"

source("settings.R")
require(ggplot2)
require(ggsci)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript FCC_fract_pie_ds.R

outFolder <- file.path("FCC_FRACT_PIE_DS")
dir.create(outFolder, recursive = TRUE)

inDT <- get(load("BARPLOT_WITH_FCC_FRACT/all_dt.Rdata"))
ds_dt <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]
stopifnot(sum(ds_dt$countFCC) == 1)
ds_dt$pct <- ds_dt$countFCC * 100
ds_dt$pct_lab <- paste0(round(ds_dt$pct, 2), " %")
ds_dt$pct_lab_gg <- paste0(round(ds_dt$pct, 2), " %")
ds_dt$pct_lab_gg[ds_dt$pct < 1] <- ""
ds_dt$labs1 <- paste0(ds_dt$intervalFCC, "\n",ds_dt$pct_lab)
ds_dt$labs <- paste0(ds_dt$intervalFCC, " (",ds_dt$pct_lab, ")")
ds_dt$labs <- gsub("FCC", "", ds_dt$labs)
ds_dt$labs <- factor(ds_dt$labs, levels = rev(as.character(ds_dt$labs)))

ggsci_pal <- "lancet"
ggsci_subpal <- ""
mycols <- eval(parse(text = paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(ds_dt$labs))


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_FCCfract_pie_withLab.", plotType))
do.call(plotType, list(outFile, height=myHeight*2, width=myWidth*2.5))
par(family=fontFamily)
pie(x=ds_dt$pct,labels = ds_dt$labs1, col=mycols,
    cex=plotCex,
    cex.main=2,
    main=paste0("TAD FCC fract."))
mtext(side=3, text = paste0(hicds, " - ", exprds), cex=plotCex)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

ds_dt$labs_2 <- as.character(ds_dt$pct_lab)
ds_dt$labs_2[ds_dt$pct < 1] <- ""

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_FCCfract_pie_withLeg.", plotType))
do.call(plotType, list(outFile, height=myHeight*2, width=myWidth*2))
par(family=fontFamily)
pie(x=ds_dt$pct,labels = ds_dt$labs_2, col=mycols,
    cex.main=2,
    cex  =plotCex,
    main=paste0("TAD FCC fract."))
legend(
  "bottomright",
  legend = rev(ds_dt$labs),
  pch=15,
  col = mycols,
  # horiz=TRUE
  bty="n"
)
mtext(side=3, text = paste0(hicds, " - ", exprds), cex=plotCex)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



pie_fract <- ggplot(ds_dt, aes(x="", y=pct, fill=labs))+ 
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0) +
  ggtitle("TAD FCC fract.", subtitle =  paste0(hicds, " - ", exprds))+
  labs(fill="", color="")+
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size=14, face="bold", family=fontFamily),
        plot.subtitle = element_text(hjust = 0.5, size=13, family=fontFamily),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
        ) +
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) +
geom_text(aes(y = pct/3 + c(0, cumsum(pct)[-length(pct)]), 
              label = pct_lab_gg), size=5, color="black")

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_FCCfract_pie_ggplot.", plotType))
ggsave(plot = pie_fract, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))




