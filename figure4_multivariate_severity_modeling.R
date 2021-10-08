## Figure 4A ##

library(caret)

# load the data - you can request our data by emailing the corresponding authors
full <- read.csv(file ='/Users/kpullen/breastmilk_full_Sx.csv',header= TRUE)
serum<-full[full$sample.type=="Serum",]
y<-serum$severity
y_name<-"maternal disease severity"
luminex_data <- serum[,c(7:110,117)]
luminex_data[] <- lapply(luminex_data, function(x) as.numeric(as.character(x)))

# quality control: remove luminex features with median below 10,000 MFI
luminex_data<-luminex_data[colMedians(data.matrix(luminex_data)) > 10000]
X <- data.frame(luminex_data)

#log transform data to make individual feature distributions closer to a normal distribution
X<-log10(X)
X <- select(X, -contains("EBOV"))
X <- data.frame(X,serum[,c(111:116,3)])

#impute missing values using k-nearest neighbors
X<-impute.knn(as.matrix(X))
X<-X[["data"]]

# normalize (z-score) 
X<-data.frame(X)
X<-scale(X)
X<-data.frame(X)
X<-select(X,c(contains("SARS"),contains("SASR"),contains("Spike"),contains("spike"),contains("Nucleocapsid"),contains("RBD"),contains("NT50")))

# ElasticNet Feature Selection

model <-train(x = data.matrix(X), y = y,method="glmnet",trControl = trainControl("LOOCV"),tuneLength = 10)

# the select_repeat and select_lasso functions below are customized versions of functions from the systemsseRology package (https://github.com/LoosC/systemsseRology)
select_repeat <- function(X, y, selector, options = list()) {
  
  # ----------------- BEGIN OPTIONS ----------------- #
  # How often it should be repeated
  if (!("n_trials" %in% names(options))) {
    options$n_trials <- 100
  }
  if (!("threshold" %in% names(options))) {
    options$threshold <- 0.8
  }
  # returns the whole data frame of how often a feature was selected
  if (!("return_count" %in% names(options))) {
    options$return_count <- FALSE
  }
  if (!("force_select" %in% names(options))) {
    options$force_select <- TRUE
  }
  # ----------------- END OPTIONS ----------------- #
  
  # vector counting how often each feature is selected
  feature_count <- rep(0, ncol(X))
  names(feature_count) <- colnames(X)
  
  # run the feature selector trials times and increment the counters
  # for the features that are selected
  for (trial in 1:options$n_trials) {
    features <- selector(X, y)
    feature_count[features] <- 1 + feature_count[features]
  }
  
  # keep those features that were selected in more than threshold
  # percent of the trials
  selected <- feature_count[-which(feature_count <= options$threshold * options$n_trials)]
  
  # if a selection is forced, return the features which are selected the most
  if (length(selected) == 0 & options$force_select) {
    selected <- feature_count[which(feature_count == max(feature_count))]
  }
  
  if (options$return_count) {
    return(list(feature_count = feature_count, sel_features = names(selected)))
  } else {
    return(names(selected))
  }
}
select_lasso <- function(X, y, options = list()) {
  # check if X is already z-scored and warn if this is not the case
  # nz_col_means <- abs(colMeans(X)) > 1e-10
  # nu_col_vars <- abs(matrixStats::colVars(X) - 1) > 1e-10
  # if (sum(nz_col_means) + sum(nu_col_vars) > 0) {
  #   message("Warning in select_lasso():")
  #   message("    X is not z-scored")
  # }
  
  # decide on type of GLM depending on type of y
  # 2-level factor -> binomial, n-level factor -> multinomial
  # vector -> gaussian, anything else gives an error
  if (is.factor(y)) {
    if (nlevels(y) == 1) {
      stop("y is a factor with only one level")
    } else if (nlevels(y) == 2) {
      fam <- "binomial"
    } else {
      fam <- "multinomial"
    }
  } else if (is.numeric(y)) {
    fam <- "gaussian"
  } else {
    stop("y must be of type factor or numeric vector")
  }
  
  # default alpha to 1 if it is not set
  if (!("alpha" %in% names(options))) {
    options$alpha <- model$bestTune$alpha
  }
  
  # cv.glmnet needs at least 3 folds, so we need at least three features
  n_samples <- ncol(X)
  if (n_samples < 3) {
    stop("select_lasso() requires more than three samples for internal cross-validation")
  }
  
  # default subfolds if its is not set
  # if there are 5 or more samples, default to 5
  # otherwise default to 3
  # while were at it, ensure that subfolds is in [3, n_samples]
  # and print a warning if it wasn't
  if (!("subfolds" %in% names(options))) {
    if (n_samples >= 5) {
      options$subfolds <- 5
    } else {
      options$subfolds <- 3
    }
  } else {
    if (options$subfolds > n_samples) {
      message("Warning in select_lasso():")
      message("    options$subfolds greater than number of samples")
      message("    setting options$subfolds = number of samples")
      options$subfolds <- n_samples
    }
    if (options$subfolds < 3) {
      message("Warning in select_lasso():")
      message("    options$subfolds was less than 3")
      message("    setting options$subfolds = 3")
      options$subfolds <- 3
    }
  }
  
  # fit an appropriate lasso model with a number of trials corresponding to
  # different values of lambda
  lasso <- glmnet::cv.glmnet(X, y, type.measure = "mse", alpha = options$alpha,
                             family = fam, type.multinomial = "grouped",
                             nfolds = options$subfolds)
  
  # lasso$lambda[k] is the value of lambda in the k-th trial
  # lasso$nzero[k] is the number of non-zero coefficients in the fitted model
  # (= number of features not including the intercept) in the k-th trial
  # things are arranged such that the first entry of nzero is 0
  # and the final entry of nzero equals the number of total features
  # lasso$cvm is the cross-validated score of the k-th trial
  
  # find the model with the smallest error that has at least one non-zero
  # coefficient other than the intercept
  indices <- which(lasso$nzero > 0)
  lambdas <- lasso$lambda[indices]
  scores <- lasso$cvm[indices]
  # if there is more than one index attaining the minimum which.min picks
  # the smallest one - this corresponds to choosing best score with the
  # minimal number of features selected
  best <- which.min(scores)
  lasso_coeffs <- coef(lasso, s = lambdas[best])
  
  if (fam == "multinomial") {
    # if the data has multiple responses, the coefficients are a matrix
    # that is returned as a list of columns. type.multinomial = "grouped"
    # forced features to be selected for all responses or for none, so we
    # can get the selected features also by only considering the first
    # column. we just replace lasso_coeffs by this column and proceed as usual
    lasso_coeffs <- lasso_coeffs[[1]]
  }
  
  # remove the intercept and the entries that are zero
  lasso_coeffs <- lasso_coeffs[-1,]
  lasso_coeffs <- lasso_coeffs[which(lasso_coeffs != 0)]
  
  # return the names of the selected features. previous code turned
  # coefficients into a vector, so use names() rather than rownames()
  return(names(lasso_coeffs))
}
opts_sel <- list(n_trials = 100, threshold = 0.35
                 , return_count = FALSE)
sel_features <- select_repeat(data.matrix(X), y, 
                              selector = select_lasso, 
                              options = opts_sel)
X_sel<-X[,sel_features]

model <- ropls::opls(X_sel, y,crossValI = 10, permI = 0,predI=2)

# the following block of code helps format the scores plot in Figure 4A
setEPS()
postscript("Figure-4A.eps")
options = list()
if (!("alpha" %in% names(options))) {
  options$alpha <- 1
}
if (!("size" %in% names(options))) {
  options$size <- 3.5
}
if (!("stroke" %in% names(options))) {
  options$stroke <- 0.5
}
# confidence levels for ellipses
if (!("level" %in% names(options))) {
  options$level <- 0.95
}
if (!("LV_ind" %in% names(options))) {
  options$LV_ind <- c(1,2)
} else if (ropls::getSummaryDF(model)$pre +
           ropls::getSummaryDF(model)$ort < max(options$LV_ind)) {
  stop("required LV exceed existing LVs")
}
# ----------------- END OPTIONS ----------------- #
# ----------------- GET SCORES ----------------- #
#
if (ropls::getSummaryDF(model)$ort > 0) {
  if (options$LV_ind[1] == 1) {
    df_scores <- data.frame(LV1 = ropls::getScoreMN(model),
                            LV2 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1],
                            y = y)
  } else {
    df_scores <- data.frame(LV1 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[1] - 1],
                            LV2 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1],
                            y = y)
  }
} else {
  df_scores <- data.frame(LV1 = ropls::getScoreMN(model)[,options$LV_ind[1]],
                          LV2 = ropls::getScoreMN(model)[,options$LV_ind[2]],
                          y = y)
}
# --------------------- END GET SCORES ------------------- #
# ---------------------- BEGIN PLOT ---------------------- #
plt_scores <-ggplot2::ggplot(df_scores,aes(LV1, LV2, color = as.factor(y)))  + geom_point(size =6)
plt_scores+
  ggplot2::geom_vline(xintercept = 0, size = 0.3) +
  ggplot2::geom_hline(yintercept = 0, size = 0.3) +
  ggplot2::labs(x = paste("scores on LV", options$LV_ind[1], "(",
                          toString(model@modelDF$R2X[options$LV_ind[1]] * 100), "%)", sep = ""),
                y = paste("scores on LV", options$LV_ind[2], "(",
                          toString(model@modelDF$R2X[options$LV_ind[2]] * 100), "%)", sep = ""),
                fill = y_name, color = y_name) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "right",
                 aspect.ratio = 1,
                 axis.text = ggplot2::element_text(color = "black"))
dev.off()

# Figure 4A Model Validation (Supplemental Fig 4A)

final<-matrix(0, 1, 18)
for (rep in 1:18){
  X_small<-X_sel[-rep,]
  y_small<-as.factor(y[-rep])
  model <- ropls::opls(X_small, y_small,crossValI = 10, permI = 0,predI=2)
  final[rep]<-predict_ropls(model, X_sel[rep,])
  
}
plotting<-data.frame(final[1,],y)
colnames(plotting)<-c("X","y")
k<-ggplot(plotting, aes(x =X, y = y))
k + geom_point(size = 6,colour="blue")+ geom_smooth(method=lm,,colour="blue")

## Figure 4B ##

# load the data - you can request our data by emailing the corresponding authors
full <- read.csv(file ='/Users/kpullen/breastmilk_full_Sx.csv',header= TRUE)
BM<-na.omit(full[full$sample.type=="BM",])
y<-BM$severity
y_name<-"maternal disease severity"
luminex_data <- BM[,c(7:110)]
luminex_data[] <- lapply(luminex_data, function(x) as.numeric(as.character(x)))

# quality control: remove luminex features with median below 10,000 MFI
luminex_data<-luminex_data[colMedians(data.matrix(luminex_data)) > 10000]
X <- data.frame(luminex_data)
X <- select(X, -contains("EBOV"))

#log transform data to make individual feature distributions closer to a normal distribution
X<-log10(X)
X <- data.frame(X,BM[,c(111:117,3)])

# normalize (z-score) 
X<-scale(X)
X<-data.frame(X)
X<-select(X,c(contains("SARS"),contains("SASR"),contains("Spike"),contains("spike"),contains("Nucleocapsid"),contains("RBD"),contains("NT50")))
X<-select(X,-c(contains("FcRn")))

# ElasticNet Feature Selection
model <-train(x = data.matrix(X), y = y,method="glmnet",trControl = trainControl("LOOCV"),tuneLength = 10)

opts_sel <- list(n_trials = 100, threshold = 0.35
                 , return_count = FALSE)
sel_features <- select_repeat(data.matrix(X), y, 
                              selector = select_lasso, 
                              options = opts_sel)
X_sel<-X[,sel_features]

model <- ropls::opls(X_sel, y,crossValI = 10, permI = 0,predI=2)

# the following block of code helps format the scores plot in Figure 4A
setEPS()
postscript("Figure-4B-Left.eps")
options = list()
if (!("alpha" %in% names(options))) {
  options$alpha <- 1
}
if (!("size" %in% names(options))) {
  options$size <- 3.5
}
if (!("stroke" %in% names(options))) {
  options$stroke <- 0.5
}
# confidence levels for ellipses
if (!("level" %in% names(options))) {
  options$level <- 0.95
}
if (!("LV_ind" %in% names(options))) {
  options$LV_ind <- c(1,2)
} else if (ropls::getSummaryDF(model)$pre +
           ropls::getSummaryDF(model)$ort < max(options$LV_ind)) {
  stop("required LV exceed existing LVs")
}
# ----------------- END OPTIONS ----------------- #
# ----------------- GET SCORES ----------------- #
#
if (ropls::getSummaryDF(model)$ort > 0) {
  if (options$LV_ind[1] == 1) {
    df_scores <- data.frame(LV1 = ropls::getScoreMN(model),
                            LV2 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1],
                            y = y)
  } else {
    df_scores <- data.frame(LV1 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[1] - 1],
                            LV2 = ropls::getScoreMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1],
                            y = y)
  }
} else {
  df_scores <- data.frame(LV1 = ropls::getScoreMN(model)[,options$LV_ind[1]],
                          LV2 = ropls::getScoreMN(model)[,options$LV_ind[2]],
                          y = y)
}
# --------------------- END GET SCORES ------------------- #
# ---------------------- BEGIN PLOT ---------------------- #
plt_scores <-ggplot2::ggplot(df_scores,aes(LV1, LV2, color = as.factor(y)))  + geom_point(size =6)
plt_scores+
  ggplot2::geom_vline(xintercept = 0, size = 0.3) +
  ggplot2::geom_hline(yintercept = 0, size = 0.3) +
  ggplot2::labs(x = paste("scores on LV", options$LV_ind[1], "(",
                          toString(model@modelDF$R2X[options$LV_ind[1]] * 100), "%)", sep = ""),
                y = paste("scores on LV", options$LV_ind[2], "(",
                          toString(model@modelDF$R2X[options$LV_ind[2]] * 100), "%)", sep = ""),
                fill = y_name, color = y_name) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "right",
                 aspect.ratio = 1,
                 axis.text = ggplot2::element_text(color = "black"))
dev.off()

loading<-as.data.frame(model@loadingMN)
loading[ "variable" ] <- rownames(loading)
loading$p2<-NULL
colnames(loading)<-c("LV1","Feature")
o<-order(abs(loading$LV1), decreasing = FALSE)
loading<-loading[o,]
loading$Feature <- factor(loading$Feature, levels = loading$Feature)
setEPS()
postscript("Figure-4B-Right.eps")
ggplot(loading,aes(x = Feature,y = LV1)) +geom_bar(stat = "identity",position = "dodge") +coord_flip()
dev.off()

# Figure 4B Model Validation (Supplemental Fig 4B)

final<-matrix(0, 1, 18)
for (rep in 1:18){
  X_small<-X_sel[-rep,]
  y_small<-(y[-rep])
  model <- ropls::opls(X_small, y_small,crossValI = 10, permI = 0,predI=2)
  final[rep]<-predict_ropls(model, X_sel[rep,])
  
}
plotting<-data.frame(final[1,],y)
colnames(plotting)<-c("X","y")
k<-ggplot(plotting, aes(x =X, y = y))
k + geom_point(size = 6,colour="blue")+ geom_smooth(method=lm,,colour="blue")


