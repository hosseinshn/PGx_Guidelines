library(ggplot2); library(ggsci); library(ggpubr)
library(gtable); library(grid); #library(ggplotify)
library(wesanderson)
library(ComplexHeatmap)
library(readxl)
library(viridis)
library(circlize)
library(reshape2)
library(fmsb)
library(GGally)
library(extrafont)
#font_import() # import all your fonts
#fonts() #get a list of fonts

### CROSS-DOMAIN
getFigData <- function()
{
  #dfall <- xlsx::read.xlsx("data/PGx_GuideLine_Figures_data.xlsx", sheetIndex = 1)
  dfall <- data.frame(read_excel("data/PGx_GuideLine_Figures_data.xlsx", sheet="Figure-5-AAC"))
  dfall #melt(dfall, id=3:4)
}

dfall <- getFigData()
rng <- range(c(dfall$GDSCv2, dfall$gCSI))
rng <- c(-0.1, 0.8)

dfall <- BBmisc::sortByCol(dfall, col=c("GDSCv2", "gCSI"))
drgOrd <- unique(dfall$Drug)
dfall_gdsc <- dfall[, c("Method", "Drug", "GDSCv2")]
dfall_gcsi <- dfall[, c("Method", "Drug", "gCSI")]

dfall_gdsc2 <- dcast(dfall_gdsc, Method~Drug)
dfall_gdsc2$Method <- factor(dfall_gdsc2$Method, levels = c("Baseline", "Ridge Regression",
                                          "Elastic Net", "Random Forest", "Deep Neural Network"))
dfall_gdsc2$dataset <- "GDSCv2"
dfall_gdsc2 <- dfall_gdsc2[, c("Method", drgOrd, 'dataset')]

dfall_gcsi2 <- dcast(dfall_gcsi, Method~Drug)
dfall_gcsi2$Method <- factor(dfall_gcsi2$Method, levels = c("Baseline", "Ridge Regression",
                                          "Elastic Net", "Random Forest", "Deep Neural Network"))
dfall_gcsi2$dataset <- "gCSI"
dfall_gcsi2 <- dfall_gcsi2[, c("Method", drgOrd, 'dataset')]

dfall_final <- rbind(dfall_gcsi2, dfall_gdsc2)

col_pal <- c('#636363','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')


dfall_final_m <- melt(dfall_final)
#Bar Plot #1 (cross domain)

png("figures/fig-5A.png", width = 8, height = 4.50, units="in",res=600)
ggplot(dfall_final_m, aes(fill=Method, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity",colour='black') + facet_wrap(~dfall_final_m$dataset)  +
  scale_colour_manual(values=c('#636363', "#4DBBD5FF","#00A087FF",'#e41a1c','#377eb8')) +
  scale_fill_manual(values=c('#636363', "#4DBBD5FF","#00A087FF",'#e41a1c','#377eb8')) + xlab("Drug") +
  ylab("Pearson") + coord_flip() + labs(x = "", fill = "Methods (Cross-Domain)")

dev.off()


#WITHIN-DOMAIN

source("data/Within-domain results and plots--Protein coding genes/within_domain_plots.R")
ps <- get_Fig_Data()
colpal=c("#8491B4FF","#4DBBD5FF", "#00A087FF", "#3C5488FF",  "#E64B35FF")[c(2,5)]

df <- ps$pearson
df <- BBmisc::sortByCol(df, "value")

drugOrd <- unique(df$drug)
methodOrd <- c("Elastic Net","Random Forest","Ridge Regression")
matOrd <- c("IC50", "AAC")

df$drug <- factor(df$drug, levels = drugOrd)
df$method <- factor(df$method, levels = methodOrd)
df$Metric <- factor(df$metric, levels = matOrd)
#pearsonPlot <- get_barError_plot(df, colpal, xt="Pearson correlation")

df <- df[which(df$metric == "AAC"),]
df <- df[,c("method","drug","value")]
colnames(df) <- c("Method","Drug","Value")

head(dfall_final)

df2 <- dcast(df, Method~Drug)
df2$Method <- factor(df2$Method, levels = c("Ridge Regression","Elastic Net", "Random Forest"))
df2_copy <- df2
df2$dataset <- "gCSI"
df2_copy$dataset <- "GDSCv2"
df2 <- df2[, c("Method", drgOrd, 'dataset')]
df2_copy <- df2_copy[, c("Method", drgOrd, 'dataset')]

df2_final <- rbind(df2, df2_copy)


dfall_final$Method <- paste0(dfall_final$Method," (Cross-Domain)")
df2_final$Method <- paste0(df2_final$Method," (Within-Domain)")


final_df <- rbind(dfall_final, df2_final)

col_pal2 <- c('#636363','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00',"#ffa500","#c76b99")

final_df$Method <- factor(final_df$Method, levels = c("Baseline (Cross-Domain)","Deep Neural Network (Cross-Domain)", 
                                                      "Elastic Net (Cross-Domain)","Random Forest (Cross-Domain)",
                                                      "Ridge Regression (Cross-Domain)",
                                                      "Elastic Net (Within-Domain)","Random Forest (Within-Domain)",
                                                      "Ridge Regression (Within-Domain)"))




###Fig 5-B,C,D (Average per method for each drug)

df3 <- rbind(df2,df2_copy)
df3$section <- "Within-Domain"
df4 <- rbind(dfall_gcsi2, dfall_gdsc2)

df4$section <- "Cross-Domain"

df3 <- rbind(df3, df4)
df3 <- df3[which(!df3$Method == "Baseline"),]

df3_m <- melt(df3)

for (i in unique(df3_m$Method)){
  
  df3_m_x <- df3_m[which(df3_m$Method == i),]
  
  png(paste0("figures/fig5-part2-",i,"_pearson.png"), width = 8, height = 4.50, units="in", res=600)
  pp <- ggplot(df3_m_x, aes(fill=section, y=value, x=variable)) + 
    geom_bar(position="stack", stat="identity",colour='black') + facet_wrap(~df3_m_x$dataset) +
    scale_colour_manual(values=c("#E64B35FF", "#4DBBD5FF","#00A087FF")) +
    scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF","#00A087FF")) + xlab("Drug") +
    ylab("Pearson") + coord_flip() + ggtitle(i) + labs(x = "", fill = "Domain")
  
  print(pp)
  dev.off()
  
}
