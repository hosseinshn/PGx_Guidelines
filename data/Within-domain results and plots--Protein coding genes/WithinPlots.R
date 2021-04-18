library(dplyr)
library(ggplot2)
require(data.table)
library(readxl)

drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")

EN <- read_excel("EN AAC.xlsx", col_names = FALSE)
EN[,1] <- drugs
RG <- read_excel("Ridge AAC.xlsx", col_names = FALSE)
RG[,1] <- drugs
RF <- read_excel("RF AAC.xlsx", col_names = FALSE)
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

tiff('AAC Within plot.png', units="in", width=10, height=7, res=600, compression = 'lzw')
ggplot(df_new, aes(x=drug, y=value, fill=method)) +
  geom_boxplot() + ggtitle("Within-domain AAC performance") +
  geom_point(position=position_jitterdodge(),alpha=0.2) +
  theme_bw(base_size = 16) + ylab("Pearson Correlation") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
  #stat_compare_means(aes(group = value), label = "p.signif")
dev.off()

ENI <- read_excel("EN IC50.xlsx", col_names = FALSE)
ENI[,1] <- drugs
RGI <- read_excel("Ridge IC50.xlsx", col_names = FALSE)
RGI[,1] <- drugs
RFI <- read_excel("RF IC50.xlsx", col_names = FALSE)
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

tiff('IC50 Within plot.png', units="in", width=10, height=7, res=600, compression = 'lzw')
ggplot(dfI_new, aes(x=drug, y=value, fill=method)) +
  geom_boxplot() + ggtitle("Within-domain IC50 performance") +
  geom_point(position=position_jitterdodge(),alpha=0.2) +
  theme_bw(base_size = 16) + ylab("Pearson Correlation") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
  #stat_compare_means(aes(group = value), label = "p.signif")
dev.off()



EN <- read_excel("EN AACS.xlsx", col_names = FALSE)
EN[,1] <- drugs
RG <- read_excel("Ridge AACS.xlsx", col_names = FALSE)
RG[,1] <- drugs
RF <- read_excel("RF AACS.xlsx", col_names = FALSE)
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

tiff('AACS Within plot.png', units="in", width=10, height=7, res=600, compression = 'lzw')
ggplot(df_new, aes(x=drug, y=value, fill=method)) +
  geom_boxplot() + ggtitle("Within-domain AAC performance") +
  geom_point(position=position_jitterdodge(),alpha=0.2) +
  theme_bw(base_size = 16) + ylab("Spearman Correlation") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
#stat_compare_means(aes(group = value), label = "p.signif")
dev.off()

ENI <- read_excel("EN IC50S.xlsx", col_names = FALSE)
ENI[,1] <- drugs
RGI <- read_excel("Ridge IC50S.xlsx", col_names = FALSE)
RGI[,1] <- drugs
RFI <- read_excel("RF IC50S.xlsx", col_names = FALSE)
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

tiff('IC50S Within plot.png', units="in", width=10, height=7, res=600, compression = 'lzw')
ggplot(dfI_new, aes(x=drug, y=value, fill=method)) +
  geom_boxplot() + ggtitle("Within-domain IC50 performance") +
  geom_point(position=position_jitterdodge(),alpha=0.2) +
  theme_bw(base_size = 16) + ylab("Spearman Correlation") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
#stat_compare_means(aes(group = value), label = "p.signif")
dev.off()

