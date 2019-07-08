library(data.table)
library(tidyverse)

df          <- fread('tot_nomis.csv')[,-1] # organized quantitative proteomics data
msi_toplist <- fread('msi_toplist.csv')[-1,-1] # lists generated after Limma analysis (https://bioconductor.org/packages/release/bioc/html/limma.html)


com     <- intersect(msi_toplist$V2,colnames(df)) 
com     <- c('sample','msi','gender',com)
df     <- df[, ..com]
df$idx <- c(1:68)
    
## normalization

df[is.na(df)] <- 0  # replace na = 0
df$msi.v      <- ifelse(df$msi=='Lo', 0,1)

minmx <- function(x){return((x- min(x)) /(max(x)-min(x)))  }


## split and normallize
set.seed(31876)
te_idx <- sample(1:68, 58, replace = FALSE)

train  <- df[-te_idx, minmx(df[te_idx,4:100])]
test   <- df[te_idx, minmx(df[-te_idx,4:100])]


# create a bootstrap data

set.seed(78562)
b <- 20 # make 100 copy
n <- 58
i <- sample(rep(te_idx, b))
g <- rep(1:b, each = n)
t <- data.table(idx=i, g=g)

train.b1 <-  merge(t, df, by= 'idx')
#train.b  <- train.b1[, minmx(train.b1[,-c('idx', 'g', 'sample', 'msi', 'msi.v',"mismatch", "gender", "msi","msi.v")])]

# PCA analysis #
dat <- train

df.pc <- prcomp(dat, center = TRUE, scale = FALSE)

summary(df.pc)

# plot
screeplot(df.pc, type = "l", npcs = 15, main = "Screeplot of the first 15 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(df.pc$sdev^2 / sum(df.pc$sdev^2))

plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.62, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)

# Scatter plot observations by components 1 and 2

plot(df.pc$x[, c(1, 3)], col = (df$msi.v + 1), 
     xlab = "PC1", ylab = "PC2")
legend(x="topright", pch=1, col = c("red", "black"), legend = c("Lo", "Hi"))

#### bi-plot #Bi-plot using covariance matrix: Looking at the descriptive statistics of "area_mean" and "area_worst"
table(df$msi.v[-te_idx]) # train (lo = 50, hi = 13)

cex.before <- par("cex")
par(cex = 0.7)
biplot(df.pc)

## Linear Discriminant Analysis (LDA) ##

# 1. extract first 10 pcs to build lda

tr.pc     <- data.table(df.pc$x[,1:10])
tr.pc$msi <- df[te_idx, 'msi.v'] # add target to the train.pc
#tr.pc$msi <- train.b1$msi.v

# 2. calculate the lda using MASS package

library(MASS)
tr.lda <- lda(msi~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +PC7 +PC8 
              + PC9 +PC10, data = tr.pc )
(tr.lda)

#### prediction

test.pc <- prcomp(test, center = TRUE, scale = FALSE)
te.pc     <- data.table(test.pc$x[,1:10])
te.pc$msi <- df[-te_idx, 'msi.v']     

te_pred <- predict(tr.lda, newdata = te.pc)
table(te_pred$class, te.pc$msi)

# Evaluate the model

library("ROCR")

te_post   <- as.data.frame(te_pred$posterior[,2])
pred      <- prediction(te_post, te.pc$msi)
roc.perf  <- performance(pred, measure = "tpr", x.measure = "fpr")

auc.train <- performance(pred, measure = "auc")
auc.train <- auc.train@y.values
plot(roc.perf)
abline(a=0, b= 1)
text(x = .40, y = .6,paste("AUC = ", round(auc.train[[1]],3), sep = ""))

