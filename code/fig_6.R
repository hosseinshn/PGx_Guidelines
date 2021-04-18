library(ggplot2); library(ggsci); library(ggpubr)
library(gtable); library(grid); #library(ggplotify)
library(wesanderson)
library(ComplexHeatmap)
library(readxl)
library(viridis)
library(circlize)
library(reshape2)
#library(fmsb)
library(GGally)

library(extrafont)
#font_import() # import all your fonts
#fonts() #get a list of fonts

getFigData <- function()
{
  #dfall <- xlsx::read.xlsx("data/PGx_GuideLine_Figures_data.xlsx", sheetIndex = 1)
  dfall <- data.frame(read_excel("data/PGx_GuideLine_Figures_data.xlsx", sheet="Figure-6"))
  dfall
}

df <- getFigData()


plotRadar <- function(dt, title)
{
dk <- rbind.data.frame(rep(1,ncol(dt)) , rep(0.0,ncol(dt)) , dt)

drord <- c("Paclitaxel", "Vorinostat","Sirolimus", "Pictilisib", 
           "Erlotinib", "Lapatinib", "Bortezomib", "Crizotinib", "Docetaxel", "Entinostat", 
           "Gemcitabine")
#colors_border <- pal_nejm("default")(4)
colors_border <- c('#636363', #"#BC3C29FF", "#0072B5FF", "#20854EFF","#fec44f")#[c(2:5)]
                   '#d95f02','#1b9e77','#7570b3')


source("R/radarChart.R")

radarchart(dk[, drord], axistype=1 , pcol=colors_border , 
           pfcol=c("#d9d9d980", NA,NA,NA),# colors_in , 
           plwd=1 , plty=1, pty=16, 
           cglcol="#d9d9d9", cglty=1, cglwd=0.5,
           axislabcol="black", 
           seg=10, caxislabels=paste0(seq(0,1,by = 0.1)," "), 
           calcex=0.5, 
           vlcex=0.8 #, vlabels="", 
           ,title =title) #centerzero = T, 
#legend(x=1, y= -0.5, legend = rownames(dk[-c(1,2),]), bty = "n", pch=20 , 
#       col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
legend(x=1.2, y= -0.5, #x = "bottom", 
       legend = rownames(dk[-c(1,2),]), #horiz = TRUE,
       bty = "n", pch = 20 , col = colors_border, text.col = "black", 
       cex = 1.2, pt.cex = 1.9)

}

df$Type <- factor(df$Type, levels = c("Baseline AUPR", "Solid samples",
                                      "Random samples", "All samples"))

png("figures/fig-6B.png", width = 9, height = 6.0, res=600, units="in")
dt <- acast(df, Type~Drug, value.var="gCSI")
plotRadar(dt, title ="gCSI")

dev.off()
png("figures/fig-6A.png", width = 9, height = 6.0, res=600, units="in")
dt <- acast(df, Type~Drug, value.var="GDSCv2")
plotRadar(dt, title ="GDSCv2")

dev.off()
