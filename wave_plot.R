
plot_cumsumDiff05 <- function(observ_vect, permut_DT, pointObsCol = "black", 
           my_stat = "ratioDown",
           my_main = NULL,
           my_ylab = NULL,
           my_xlab = NULL,
           cexMain = 1,
           polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),
           departureValue = 0.5, drawline=FALSE) {
  
  plotTrueType <- ifelse(drawline, "l", "p")

  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from ", departureValue)
  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)

  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  
  x_val <- c(1:length(observ_vect))
  diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  
  plot(cumsum(abs(observ_vect - departureValue)) ~ x_val,
       main= my_main,
       cex.main = cexMain,
       type = plotTrueType,
       pch = 16, cex = 0.7,
       xlab= my_xlab, 
       ylab= my_ylab,
       bty="l")
  polygon(x = c(x_val, rev(x_val)), 
          y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
          border=NA,
          col = polygonPermutCol)
  legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
         pt.cex = c(0.7, 2),
         col = c(pointObsCol, polygonPermutCol),
         bty="n")
}
