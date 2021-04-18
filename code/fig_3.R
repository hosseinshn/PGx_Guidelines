library(ggplot2); library(ggsci); library(ggpubr)
library(gtable); library(grid); #library(ggplotify)
library(wesanderson)

library(readxl)
coolors_co_pal <- function(url) {
  ## get color from https://coolors.co url
  #url %>%
  s = stringr::str_extract_all(url, "([a-f0-9]){6}", simplify = TRUE) #%>%
  stringr::str_pad(s, 7, "left", "#")
}


getFigData <- function()
{
#dfall <- xlsx::read.xlsx("PGx_GuideLine_Figures_data.xlsx", sheetIndex = 1)
dfall <- data.frame(read_excel("data/PGx_GuideLine_Figures_data.xlsx", sheet="Figure-3"))
dfall <- dfall[!apply(dfall[,1:4], 1, function(i)all(is.na(i))), ]
updateCol <- function(v) { as.numeric(gsub("nan", "", v)) }
dfall$CTRPv2<- updateCol(dfall$CTRPv2)
dfall$GDSCv2<- updateCol(dfall$GDSCv2)
dfall$gCSI  <- updateCol(dfall$gCSI)

return(dfall)
}
####----------------------------------

getPlot <- function(df, valCol, pal, ylm=c(-1,1), title="")
{
  df$val <- df[, valCol]

  df$Gene <- droplevels(df$Gene)
  df$xpos <- as.numeric(df$Gene)
  #df$Type <- factor(df$Type)
  typeCol <- setNames(pal[1:5], levels(df$Type))
  geneCol <- setNames(rep_len(c("#d9d9d9", "#f7f7f7"), length.out = length(levels(df$Gene))),
                      levels(df$Gene))
  cpal <- c(typeCol, geneCol)

  plt <- ggplot(data = df, aes(x = xpos, y = val, fill = Type))
  plt <- plt +  geom_rect(aes(xmin = xpos - .5, xmax = xpos + .5,
                  ymin = -Inf, ymax = Inf, fill = Gene),alpha = 0.95)
  #plt + geom_col(position = position_dodge(.5), width = .5)
  plt <- plt + geom_bar(stat="identity", position=position_dodge(0.85),
                        width = 0.75#, colour = "#f7f7f7"
                        )
  plt <- plt+scale_fill_manual(values = cpal, breaks = rev(names(typeCol)))

  plt <- plt + theme_classic()
  plt <- plt + scale_x_continuous(labels = levels(df$Gene),
                           breaks = seq(1, length(unique(as.character(df$Gene)))),
                           expand = c(0, 0))

  #plt <- plt + ylim(ylm) #+ geom_hline(yintercept=0, color = "#737373", size=0.5)

  prettyZero <- function(l){
    max.decimals = max(nchar(stringr::str_extract(l, "\\.[0-9]+")), na.rm = T)-1
    lnew = formatC(l, replace.zero = T, zero.print = "0",
                   digits = max.decimals, format = "f", preserve.width=T)
    return(lnew)
  }

  plt <- plt + scale_y_continuous(limits =ylm, expand = c(0, 0),
                                  labels = prettyZero)

  plt <- plt + facet_grid(. ~ title)+ theme_bw() +
         theme(strip.background =element_rect(fill="#EDDEB8"),
               strip.text = element_text(colour = 'black', face = "bold",size=10))
  plt <- plt + ylab("Correlation")+ xlab("")
  plt + coord_flip()
}

##-----------------------------------------
dfall <- getFigData()
dfall <- do.call(rbind.data.frame,
               lapply(c("CTRPv2", "GDSCv2", "gCSI"), function(i){
                 dft <- dfall[, c("Gene", "Type", "Drug", i)]; dft$title<- i;
                 colnames(dft)[4]<-"value";dft}))
rownames(dfall)<- NULL

dfall$value[is.na(dfall$value)] <- 0
dfall$value <- abs(dfall$value)

lev <- c( "IC50", "log(IC50)", "trc. IC50", "log(trc. IC50)", "AAC")
dfall$Type <- factor(dfall$Type, levels = lev)

dfall <- BBmisc::sortByCol(dfall, c("Drug", "value"))

dfall$gd <- paste0(dfall$Gene, "\n(in ", dfall$Drug, ")")
dfall$gd <- factor(dfall$gd)#, levels = c("EGFR-Erlotinib","ERBB2-Lapatinib"))

dfall$Gene <- dfall$gd

pal=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
      "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF")[c(1,3,4,9,10)]

pal = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')
pal = coolors_co_pal("https://coolors.co/335c67-fff3b0-e09f3e-9e2a2b-540b0e")

pal = pal_npg("nrc")(5)
pal=rev(c(#"#DC0000FF",
          "#E64B35FF",
          "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#8491B4FF"))# "#91D1C2FF")

pal=c("#8491B4FF","#4DBBD5FF", "#00A087FF", "#3C5488FF",  "#E64B35FF")
ylm = c(0, 0.5) #max(dfall$value))
plt <- getPlot(dfall, valCol="value", pal, ylm)
plt <- plt+ theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
           axis.ticks.y = element_blank(), panel.spacing.x = unit(1, "lines"))

plt <- plt+ theme(axis.text = element_text(face = "bold", color = "black"))
#plt
##8.27 × 11.69 inches
png("figures/fig-3A.png", width = 8, height = 3.0, units = "in" ,res = 600)
print(plt)
dev.off()
##-----------------------------------------
##----for AAC and IC50 comparision fig. ---

get_barError_plot <- function(df, colpal, xt="Correlation")
{
df$xpos <- as.numeric(df$drug)
#df$Type <- factor(df$Type)
typeCol <- setNames(colpal[1:2], levels(df$Metric))
drgCol <- setNames(rep_len(c("#d9d9d9", "#f7f7f7"), length.out = length(levels(df$drug))),
                    levels(df$drug))
cpal <- c(typeCol, drgCol)

plt <- ggplot(data = df, aes(x = xpos, y = value, fill = Metric))
plt <- plt +  geom_rect(aes(xmin = xpos - .5, xmax = xpos + .5,
                            ymin = -Inf, ymax = Inf, fill = drug),alpha = 0.95)
#plt + geom_col(position = position_dodge(.5), width = .5)
plt <- plt + geom_bar(stat="identity", position=position_dodge(0.85),
                      width = 0.75)
plt <- plt+geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(.9), color="#525252")
plt <- plt+scale_fill_manual(values = cpal, breaks = rev(names(typeCol)))
plt <- plt + theme_classic()
plt <- plt + facet_grid(. ~ method, scales = "free_x")+ theme_bw() +
  theme(strip.background =element_rect(fill="#EDDEB8"),
        strip.text = element_text(colour = 'black', face = "bold",size=9))
plt <- plt + ylab(xt)+ xlab("")
plt <- plt + scale_x_continuous(labels = levels(df$drug),
                                breaks = seq(1, length(unique(as.character(df$drug)))),
                                expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

plt <- plt+ theme(#axis.text.y = element_text(angle = 90, hjust = 0.5),
                  axis.ticks.y = element_blank(), panel.spacing.x = unit(0.75, "lines"))
plt <- plt+ theme(axis.text = element_text(face = "bold", color = "black"))
plt + coord_flip()
}

source("R/within_domain_plots.R")
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
pearsonPlot <- get_barError_plot(df, colpal, xt="Pearson correlation")

##
ds <- ps$spearman
ds$drug <- factor(ds$drug, levels = drugOrd)
ds$method <- factor(ds$method, levels = methodOrd)
ds$Metric <- factor(ds$metric, levels = matOrd)
spearmanPlot <- get_barError_plot(ds, colpal, xt="Spearman correlation")

#####

##8.27 × 11.69 inches
png("figures/fig-3B_AAC-IC50_Pearson.png", width = 8, height = 4.50, units = "in", res = 600)
print(pearsonPlot)
dev.off()

####-------------------
####-------------------
#df <- dfall[dfall$Drug=="Lapatinib",]
#df <- BBmisc::sortByCol(df, c("value"))
#df$Gene <- factor(df$Gene, levels = c("ERBB2", "ERBB3", "EGFR", "IGF1R", "PIK3CA", "MET"))
