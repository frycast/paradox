#### Setup and libraries -----------------------------------------------------

# devtools::install_github("suren-rathnayake/EMMIX")
library(EMMIX)
source('sslfunctions.R')

#### --------------------------------------------------------------------
##
## SAR Data is available on the Zenodo archive
## https://zenodo.org/record/4008883
##
## The original pre-processed data:
## "https://zenodo.org/record/4008883/files/CC_sub_norm.tif?download=1"
## "https://zenodo.org/record/4008883/files/atlas.tif?download=1"
##
## Description of further processing:
##
## The 2x2 log-Cholesky transformed 1-lag autocovariance matrix 
## of each pixel (over time) in the sequence of SAR images was obtained,
## after down sampling using mean aggregation within
## 10x10 patches of pixels.
## 
## Each 2x2 log-Chol matrix was then 
## transformed to a 3-tuple, using equation (6) 
## from https://arxiv.org/pdf/2008.03454.pdf
##
## The above processing occurred within the script R/SAR_application.R in
## https://github.com/frycast/kmspd
## With output filname: cov_chol_d10_m1_cc.Rds
## Labels from GDE atlas: labsp_d10.Rds
##
#### --------------------------------------------------------------------

## This is the processed dataset, ready for fitting
chol <- readRDS("cov_chol_d10_m1_cc.Rds")
labs <- as.integer(readRDS("labsp_d10.Rds") > 0.5)

## Choosing the missing labels
m <- rbinom(n, size= 1, prob=0.1) # Random for testing only

# Take a small sample for testing purposes in case it's slow
n <- 1000
samp <- sample(1:nrow(chol), n)
X <- chol[samp,]                     # observations
Z <- makelabelmatrix(labs[samp]+1L)  # labels one-hot

g <- ncol(Z) # Number of clusters
p <- ncol(X) # Number of features

# get initial values from unsupervised fit
inits <- EMMIX(X, g=g, distr = 'mvn')

# list2par involves Cholesky decomposition of the g estimated Sigmas
par <- list2par(inits, g, p)

X1 <- X[m==0,] # n times p matrix of features for labeled observations
Z1 <- Z[m==0,]  # n times g indicator matrix of class labels for labeled observations
X2 <- X[m==1,] # n times p matrix of features for unlabeled observations

# run optimisation and store objective function values
syst <- system.time({
  opttrace <- capture.output(optfit <- stats::nlminb(
    start = par, objective = neg_profile_objective_function, 
    gradient=NULL, g=g, X1=X1, Z1=Z1, X2=X2, m=m, 
    control=list(
      trace=TRUE, eval.max=500, iter.max=500, rel.tol=1e-10)))
}); syst

# getting objective function values
objval <- sapply(opttrace, function(x) {
  as.numeric(strsplit(x, ':')[[1]][2])})
names(objval) <- 1:length(objval)

# updates seem to decrease objective at each iteration
plot(objval, type='o', pch=19)

# get parameter estimates in list form
parhat <- par2list(optfit$par, g, p)
parhat
