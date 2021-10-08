library(corrplot)
library(stats)
library(reshape2)

# load the data - you can request our data by emailing the corresponding authors

full <- read.csv(file ='/Users/kpullen/breastmilk_covid_nt50.csv',header= TRUE)
luminex_data <- full[,c(4:108)]
luminex_data[] <- lapply(luminex_data, function(x) as.numeric(as.character(x)))
X<-data.frame(luminex_data)
X <- select(X, -contains("EBOV"))

# correct for unequal dilutions for fc effector function assays (1:10 for all except for ADCD which is 1:5) ##

functions<-full[,c(109:115)]
functions[full$sample.type=="Serum",c(1:6)]<-functions[full$sample.type=="Serum",c(1:6)]*10
functions[full$sample.type=="Serum",7]<-functions[full$sample.type=="Serum",7]*5
X<-data.frame(X,functions)
X<-impute.knn(as.matrix(X))
X<-X[["data"]]
X<-data.frame(X)
X$FcaR.SARS.2.RBD<-NULL
breastmilk <- X[c(26:44),]
serum <- X[c(70:88),]
breastmilk<-data.frame(breastmilk)
serum<-data.frame(serum)

# log transform and z-score the data to normalize ##

breastmilk[] <- lapply(breastmilk, function(x) as.numeric(as.character(x)))
serum[] <- lapply(serum, function(x) as.numeric(as.character(x)))
BM_logged<-log10(breastmilk)
s_logged<-log10(serum)
BM<-scale(BM_logged)
BM<-data.frame(BM)
s<-scale(s_logged)
s<-data.frame(s)
BM<-select(BM,c(contains("SARS"),contains("SASR"),contains("Spike"),contains("spike"),contains("Nucleocapsid"),contains("nt50")))
s<-select(s,c(contains("SARS"),contains("SASR"),contains("Spike"),contains("spike"),contains("Nucleocapsid"),contains("nt50")))
BM<-select(BM,-contains("FcRn"))
s<-select(s,-contains("FcRn"))

# the following function will calculate the p-value of the spearman correlation between two features (derived from Stackoverflow)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method ="spearman", ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

## Figure 3A ##

serum<-cor(s, method ="spearman")
p.mat_serum <- cor.mtest(s)
p.mat_serum<-melt(p.mat_serum)

# perform multiple hypothesis correction
p.mat_serum$value<-p.adjust(p.mat_serum$value,method="bonferroni", n = 1764)
p.mat_serum<-dcast(p.mat_serum, Var1 ~ Var2)
rownames(p.mat_serum)<-p.mat_serum$Var1
p.mat_serum$Var1<-NULL
p.mat_serum<-p.mat_serum[row.names(serum),colnames(serum)]

#plot correlation matrix
corrplot(serum, method = 'color',col= viridis::magma(100, alpha = 1, begin = 0, end = 1, direction = 1),p.mat =data.matrix(p.mat_serum),
         sig.level = 0.05,
         insig = "label_sig")

## Figure 3B ##

BM_cor <-cor(BM, method ="spearman")
p.mat_BM <- cor.mtest(BM)
p.mat_BM<-melt(p.mat_BM)

# perform multiple hypothesis correction
p.mat_BM$value<-p.adjust(p.mat_BM$value,method="bonferroni", n = 1764)
p.mat_BM<-dcast(p.mat_BM, Var1 ~ Var2)
rownames(p.mat_BM)<-p.mat_BM$Var1
p.mat_BM$Var1<-NULL
p.mat_BM<-p.mat_BM[row.names(BM_cor),colnames(BM_cor)]

#plot correlation matrix
corrplot(BM_cor, method = 'color',col= viridis::magma(100, alpha = 1, begin = 0, end = 1, direction = 1),p.mat =data.matrix(p.mat_BM),
         sig.level = 0.05,
         insig = "label_sig")

## Figure 3C ##

colnames(BM) <- paste("BM", colnames(BM), sep = "_")
big<-data.frame(BM,s)

M<-cor(big, method ="spearman")
M<-M[c(1:42),c(43:84)]
p.mat <- cor.mtest(big)
p.mat<-p.mat[c(1:42),c(43:84)]
adjusted_p <- matrix(0, 42, 42)

# perform multiple hypothesis correction
for (i in 1:42){
  adjusted_p[,i]<-p.adjust(p.mat[,i], method ="bonferroni", n =1764) 
}
adjusted_p<-data.frame(adjusted_p)
colnames(adjusted_p)<-colnames(M)
row.names(adjusted_p)<-row.names(M)

#plot correlation matrix
corrplot(M, method = 'color',col= viridis::magma(100, alpha = 1, begin = 0, end = 1, direction = 1),p.mat =data.matrix(adjusted_p),
         sig.level = 0.05,
         insig = "label_sig")
