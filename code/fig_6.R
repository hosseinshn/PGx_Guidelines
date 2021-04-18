library(ggplot2); library(ggsci); library(ggpubr)
library(gtable); library(grid); #library(ggplotify)
library(wesanderson)
library(ComplexHeatmap)
library(readxl)
library(viridis)
library(circlize)
library(reshape2)
#library(fmsb)
source("R/radarChart.R")

coolors_co_pal <- function(url) {
  ## get color from https://coolors.co url
  #url %>%
  s = stringr::str_extract_all(url, "([a-f0-9]){6}", simplify = TRUE) #%>%
  stringr::str_pad(s, 7, "left", "#")
}

getFigData <- function()
{
  #dfall <- xlsx::read.xlsx("data/PGx_GuideLine_Figures_data.xlsx", sheetIndex = 1)
  dfall <- data.frame(read_excel("data/PGx_GuideLine_Figures_data.xlsx", sheet="Figure-4-AAC"))
  return(dfall)
  if(1==2){
  dfall$GDSCv2[dfall$GDSCv2=="nan"] <- NA
  dfall$gCSI[dfall$gCSI=="nan"] <- NA
  dfall$GDSCv2 <- as.numeric(dfall$GDSCv2)
  dfall$gCSI <- as.numeric(dfall$gCSI)

  dfall$method <- paste0(dfall$Method, "#", dfall$Metric)

  gdsc <- dfall[, c("Drug", "method", "GDSCv2")]
  gdsc <- acast(gdsc, Drug ~ method)
  colnames(gdsc) <- paste0("GDSC", "#", colnames(gdsc))

  gcsi <- dfall[, c("Drug", "method", "gCSI")]
  gcsi <- acast(gcsi, Drug ~ method)
  colnames(gcsi) <- paste0("gCSI", "#", colnames(gcsi))

  mat <- cbind(gdsc, gcsi[rownames(gdsc),])
  df <- data.frame(t(sapply(colnames(mat), function(i)strsplit(i, "#")[[1]])))
  colnames(df) <- c("Dataset", "Method", "Metric")

  return(list(mat=mat, df=df))
  }
}

###------------------------------------------
#dall = getFigData()
#mat=dall$mat; df=dall$df

df = getFigData(); df$Metric <- NULL
dt <- acast(df, Dataset~Drug, value.var="gCSI")
dt <- dt[, sort(colnames(dt))]
#dt <- dcast(df, Dataset~Drug, value.var="gCSI")
#colnames(dt) <- gsub("Dataset", "group", colnames(dt))

dt <- dt + 0.1 #abs(min(dt))
dk <- rbind.data.frame(rep(1,ncol(dt)) , rep(0.0,ncol(dt)) , dt)

colors_border <- pal_jama("default")(4)
colors_border <- pal_nejm("default")(4)
colors_border[3] <- "#fec44f"

colors_border <- c("#BC3C29FF", "#0072B5FF", "#fec44f", "#20854EFF")[c(2,1,3,4)]

drord <- c("Paclitaxel", "Gemcitabine", "Pictilisib", "Sirolimus",
           "Vorinostat", "Bortezomib", "Crizotinib", "Entinostat", "Docetaxel",
           "Erlotinib", "Lapatinib")

png("figures/fig-4.png", width = 9, height = 6.0, res=600, units = "in")
radarchart(dk[, drord], axistype=1 , pcol=colors_border ,
           #pfcol=c(colors_border[1], NA,NA,NA),# colors_in ,
           plwd=1 , plty=1, pty=16,
           cglcol="#d9d9d9", cglty=1, cglwd=0.5,
           axislabcol="black",
           seg=10, caxislabels=paste0(seq(-0.1,0.9,by = 0.1)," "),
           calcex=0.5,
           vlcex=0.8 #, vlabels="",
           ) #centerzero = T,
#legend(x=1, y= -0.5, legend = rownames(dk[-c(1,2),]), bty = "n", pch=20 ,
#       col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
legend(x=1.2, y= -0.5, #x = "bottom",
       legend = rownames(dk[-c(1,2),]), #horiz = TRUE,
       bty = "n", pch = 20 , col = colors_border, text.col = "black",
       cex = 0.85, pt.cex = 1.20)

dev.off()
