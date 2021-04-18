
get_Fig_Data <- function(){
  library(dplyr)
  library(ggplot2)
  require(data.table)
  library(readxl)
  library(RColorBrewer)
  library(dplyr)
  
  read_excel_pearson <- function(f, col_names = FALSE)
  {
    dir <- 'data/Within-domain results and plots--Protein coding genes/Pearson'
    read_excel(sprintf("%s/%s", dir, f), col_names = col_names)
  }
  #Pearson
  #setwd("data/Within-domain results and plots--Protein coding genes/Pearson")
  drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")
  
  #AAC
  
  EN <- read_excel_pearson("EN AAC.xlsx", col_names = FALSE)
  EN[,1] <- drugs
  RG <- read_excel_pearson("Ridge AAC.xlsx", col_names = FALSE)
  RG[,1] <- drugs
  RF <- read_excel_pearson("RF AAC.xlsx", col_names = FALSE)
  RF[,1] <- drugs
  
  Ridge <- replicate(11, "Ridge Regression")
  Elastic_Net <- replicate(11, "Elastic Net")
  Random_Forest <- replicate(11, "Random Forest")
  
  df1 <- cbind(EN, Elastic_Net)
  colnames(df1)[12]<-"method"
  colnames(df1)[1]<-"drug"
  
  df2 <- cbind(RG, Ridge)
  colnames(df2)[12]<-"method"
  colnames(df2)[1]<-"drug"
  
  df3 <- cbind(RF, Random_Forest)
  colnames(df3)[12]<-"method"
  colnames(df3)[1]<-"drug"
  
  
  df<- rbind(df1,df2,df3)
  df_new <- melt(df, id=c("drug", "method"))
  
  df_new2 <- df_new[,c("drug","method","value")]
  
  mean_drug <- as.data.frame(df_new2 %>%
                               group_by(drug, method) %>% 
                               summarise_at(vars("value"), mean))
  
  std_drug <- as.data.frame(df_new2 %>%
                              group_by(drug, method) %>% 
                              summarise_at(vars("value"), sd))
  
  mean_drug$sd <- std_drug$value
  
  
  #IC50
  ENI <- read_excel_pearson("EN IC50.xlsx", col_names = FALSE)
  ENI[,1] <- drugs
  RGI <- read_excel_pearson("Ridge IC50.xlsx", col_names = FALSE)
  RGI[,1] <- drugs
  RFI <- read_excel_pearson("RF IC50.xlsx", col_names = FALSE)
  RFI[,1] <- drugs
  
  df1I <- cbind(ENI, Elastic_Net)
  colnames(df1I)[12]<-"method"
  colnames(df1I)[1]<-"drug"
  
  df2I <- cbind(RGI, Ridge)
  colnames(df2I)[12]<-"method"
  colnames(df2I)[1]<-"drug"
  
  df3I <- cbind(RFI, Random_Forest)
  colnames(df3I)[12]<-"method"
  colnames(df3I)[1]<-"drug"
  
  
  dfI<- rbind(df1I,df2I,df3I)
  dfI_new <- melt(dfI, id=c("drug", "method"))
  
  
  dfI_new2 <- dfI_new[,c("drug","method","value")]
  
  mean_drugI <- as.data.frame(dfI_new2 %>%
                                group_by(drug, method) %>% 
                                summarise_at(vars("value"), mean))
  
  std_drug <- as.data.frame(dfI_new2 %>%
                              group_by(drug, method) %>% 
                              summarise_at(vars("value"), sd))
  
  mean_drugI$sd <- std_drug$value
  
  mean_drug$metric <- "AAC"
  mean_drugI$metric <- "IC50"
  
  mean_combined <- rbind(mean_drug, mean_drugI)
  mean_combined_Pearson <- mean_combined
  
  if(1==2){
    plt <- ggplot(mean_combined, aes(x=drug, y=value, fill=as.factor(metric))) +
      geom_bar(position=position_dodge(), stat="identity", colour='white') +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9), color="#525252") +
      labs(x = "", fill = "metric") + xlab("Drug") +  ylab("Mean Pearson")
    
    pal=c("#8491B4FF","#4DBBD5FF", "#00A087FF", "#3C5488FF",  "#E64B35FF")[c(5,2)]
    plt <- plt + scale_colour_manual(values=pal) + scale_fill_manual(values=pal)
    plt <- plt + coord_flip() 
    plt <- plt + facet_wrap(~ method, scales = "free_x")+
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         strip.background =element_rect(fill="#EDDEB8"),
                         strip.text = element_text(colour = 'black', face = "bold",size=10))
  }
  
  
  #####----------------------------------
  #Spearman
  
  #AAC
  
  read_excel_spearman <- function(f, col_names = FALSE)
  {
    dir <- 'data/Within-domain results and plots--Protein coding genes/Spearman'
    read_excel(sprintf("%s/%s", dir, f), col_names = col_names)
  }
  
  #setwd("~/Desktop/Within-domain results and plots--Protein coding genes/Spearman")
  EN <- read_excel_spearman("EN AACS.xlsx", col_names = FALSE)
  EN[,1] <- drugs
  RG <- read_excel_spearman("Ridge AACS.xlsx", col_names = FALSE)
  RG[,1] <- drugs
  RF <- read_excel_spearman("RF AACS.xlsx", col_names = FALSE)
  RF[,1] <- drugs
  
  Ridge <- replicate(11, "Ridge Regression")
  Elastic_Net <- replicate(11, "Elastic Net")
  Random_Forest <- replicate(11, "Random Forest")
  
  df1 <- cbind(EN, Elastic_Net)
  colnames(df1)[12]<-"method"
  colnames(df1)[1]<-"drug"
  
  df2 <- cbind(RG, Ridge)
  colnames(df2)[12]<-"method"
  colnames(df2)[1]<-"drug"
  
  df3 <- cbind(RF, Random_Forest)
  colnames(df3)[12]<-"method"
  colnames(df3)[1]<-"drug"
  
  
  df<- rbind(df1,df2,df3)
  df_new <- melt(df, id=c("drug", "method"))
  df_new2 <- df_new[,c("drug","method","value")]
  
  mean_drug <- as.data.frame(df_new2 %>%
                               group_by(drug, method) %>% 
                               summarise_at(vars("value"), mean))
  
  std_drug <- as.data.frame(df_new2 %>%
                              group_by(drug, method) %>% 
                              summarise_at(vars("value"), sd))
  
  mean_drug$sd <- std_drug$value
  
  
  #IC50
  
  ENI <- read_excel_spearman("EN IC50S.xlsx", col_names = FALSE)
  ENI[,1] <- drugs
  RGI <- read_excel_spearman("Ridge IC50S.xlsx", col_names = FALSE)
  RGI[,1] <- drugs
  RFI <- read_excel_spearman("RF IC50S.xlsx", col_names = FALSE)
  RFI[,1] <- drugs
  
  df1I <- cbind(ENI, Elastic_Net)
  colnames(df1I)[12]<-"method"
  colnames(df1I)[1]<-"drug"
  
  df2I <- cbind(RGI, Ridge)
  colnames(df2I)[12]<-"method"
  colnames(df2I)[1]<-"drug"
  
  df3I <- cbind(RFI, Random_Forest)
  colnames(df3I)[12]<-"method"
  colnames(df3I)[1]<-"drug"
  
  
  dfI<- rbind(df1I,df2I,df3I)
  dfI_new <- melt(dfI, id=c("drug", "method"))
  
  dfI_new2 <- dfI_new[,c("drug","method","value")]
  
  mean_drugI <- as.data.frame(dfI_new2 %>%
                                group_by(drug, method) %>% 
                                summarise_at(vars("value"), mean))
  
  std_drug <- as.data.frame(dfI_new2 %>%
                              group_by(drug, method) %>% 
                              summarise_at(vars("value"), sd))
  
  mean_drugI$sd <- std_drug$value
  
  mean_drug$metric <- "AAC"
  mean_drugI$metric <- "IC50"
  
  mean_combined <- rbind(mean_drug, mean_drugI)
  mean_combined_Spearman <- mean_combined
  if(1==2){
    ggplot(mean_combined, aes(x=drug, y=value, fill=as.factor(metric))) +
      geom_bar(position=position_dodge(), stat="identity", colour='black') +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) + coord_flip() +
      labs(x = "", fill = "metric") +
      scale_colour_manual(values=c("#A6CEE3", "#B2DF8A")) +
      scale_fill_manual(values=c("#A6CEE3", "#B2DF8A")) + xlab("Drug") +
      ylab("Mean Spearman") + facet_wrap(~ method, scales = "free_x")
  }
  
  return(list(pearson=mean_combined_Pearson, spearman=mean_combined_Spearman))
}