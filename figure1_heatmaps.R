library(gplots)
library(robustbase)
library(dplyr)
library(impute)
library(viridis)

# load the data - you can request our data by emailing the corresponding authors
full <- read.csv(file ='/Users/kpullen/breastmilk_covid_nt50.csv',header= TRUE)
luminex_data <- full[,c(5:108)] #CHANGE COLUMN NUMBER
luminex_data[] <- lapply(luminex_data, function(x) as.numeric(as.character(x)))

# account for 20x greater dilution of serum than breastmilk (assuming linearity)
luminex_data[full$sample.type=="Serum",]<-luminex_data[full$sample.type=="Serum",]*20

# quality control: remove luminex features with median below 10,000 MFI
luminex_data<-luminex_data[colMedians(data.matrix(luminex_data)) > 10000]
X<-data.frame(luminex_data)
X <- select(X, -contains("EBOV"),-contains('FcRn'))
X<-select(X,c(contains("SARS"),contains("SASR"),contains("Spike"),contains("spike"),contains("Nucleocapsid")))

#log transform data to make individual feature distributions closer to a normal distribution
X<-log10(X)

#separate serum and breastmilk (BM) samples
SERUM<-X[full$sample.type=="Serum",]
BM<-X[full$sample.type=="BM",]

#impute missing values using k-nearest neighbors
SERUM<-impute.knn(as.matrix(SERUM))
SERUM<-SERUM[["data"]]
BM<-impute.knn(as.matrix(BM))
BM<-BM[["data"]]

# normalize (z-score) the breastmilk and serum data separately 
SERUM<-scale(SERUM, center = TRUE, scale = TRUE)
BM<-scale(BM, center = TRUE, scale = TRUE)
BM$FcaR.SARS.2.RBD<-NULL

# the "y" variable refers to SARS-CoV-2 infection status of the nursing individual
y<-full$status[1:44]
y<-as.factor(y)
BM_pos<-BM[y=='pos',]
BM_neg<-BM[y=='neg',]
X<-data.frame(SERUM)
SERUM$FcaR.SARS.2.RBD<-NULL
SERUM_pos<-SERUM[y=='pos',]
SERUM_neg<-SERUM[y=='neg',]

#Create SARS-CoV-2 Negative Serum Heatmap
dev.off()
heatmap.2(data.matrix(SERUM_neg),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col = viridis(100,option = "magma"),margins = c(5,28),breaks = seq(-4, 6, length.out = 101))

#Create SARS-CoV-2 Positive Serum Heatmap
dev.off()
heatmap.2(data.matrix(SERUM_pos),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col = viridis(100,option = "magma"),margins = c(5,28),breaks = seq(-4, 6, length.out = 101))

#Create SARS-CoV-2 Negative Breastmilk Heatmap
dev.off()
heatmap.2(data.matrix(BM_neg),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col = viridis(100,option = "magma"),margins = c(5,28),breaks = seq(-4, 6, length.out = 101))

#Create SARS-CoV-2 Positive Breastmilk Heatmap
dev.off()
heatmap.2(data.matrix(BM_oos),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col = viridis(100,option = "magma"),margins = c(5,28),breaks = seq(-4, 6, length.out = 101))


