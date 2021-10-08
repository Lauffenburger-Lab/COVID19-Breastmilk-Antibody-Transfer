## Figure 4C ##

library(ggplot2)
library(reshape2)

# load the data - you can request our data by emailing the corresponding authors

full <- read.csv(file ='/Users/kpullen/breastmilk2_no_asymp_Sx.csv',header= TRUE)
serum<-full[c(15:28),]
days<-serum$D.from.Sx
y<-days

# data preprocessing (QC control & log transformation. No normalization necessary for these 1-Feature Linear Regressions)

luminex_data <- full[,c(7:110)] #CHANGE COLUMN NUMBER
luminex_data[] <- lapply(luminex_data, function(x) as.numeric(as.character(x)))
luminex_data<-luminex_data[colMedians(data.matrix(luminex_data)) > 10000]
X<-data.frame(luminex_data)
X <- select(X, -contains("EBOV"))
serum <- X[c(15:28),]
serum[] <- lapply(serum, function(x) as.numeric(as.character(x)))
serum <- data.frame(serum,full[c(15:28),117])
s_logged<-log10(serum*20)
s <- data.frame(s_logged,full[c(15:28),c(111:116,3)])
s<-data.frame(s)
colnames(s)[71]<-'ADCD.Spike'
s<-select(s,c(contains("SARS"),contains("SASR"),contains("Spike"),contains("spike"),contains("Nucleocapsid"),contains("RBD"),contains("NT50")))
s<-select(s,c(contains("IgA1"),contains("IgG"),contains("IgM")))

bm <- X[c(1:14),]
bm[] <- lapply(bm, function(x) as.numeric(as.character(x)))
bm <- data.frame(bm,full[c(1:14),117])
bm_logged<-log10(bm)
bm <- data.frame(bm_logged,full[c(1:14),c(111:116,3)])
bm<-data.frame(bm)
colnames(bm)[71]<-'ADCD.Spike'
bm<-select(bm,c(contains("SARS"),contains("SASR"),contains("Spike"),contains("spike"),contains("Nucleocapsid"),contains("RBD"),contains("NT50")))
bm<-select(bm,c(contains("IgA1"),contains("IgG"),contains("IgM")))

# data formatting to be inputted into ggplot()

IgA<-data.frame(s$IgA1.SARS.2.Spike,bm$IgA1.SARS.2.Spike)
IgG<-data.frame(s$IgG1.SARS.2.Spike,bm$IgG1.SARS.2.Spike)
IgM<-data.frame(s$IgM.SARS.2.Spike,bm$IgM.SARS.2.Spike)
IgA['y']<-y
IgG['y']<-y
IgM['y']<-y
colnames(IgA)<-c('serum','bm','days')
colnames(IgG)<-c('serum','bm','days')
colnames(IgM)<-c('serum','bm','days')

#plot IgA Regressions
IgA_melt<-melt(IgA,id='days')
ggplot(IgA_melt, aes(x=days, y=value, color=variable)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ylim(3.75, 7.25)

#plot IgG Regressions
IgG_melt<-melt(IgG,id='days')
ggplot(IgG_melt, aes(x=days, y=value, color=variable)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ylim(3.75, 7.25)

#plot IgM Regressions
IgM_melt<-melt(IgM,id='days')
ggplot(IgM_melt, aes(x=days, y=value, color=variable)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ylim(3.75, 7.25)
