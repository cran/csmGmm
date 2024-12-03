## -----------------------------------------------------------------------------
library(csmGmm)
library(dplyr)

# number of SNPs and proportion in each case
J <- 40000
K <- 2
case0 <- 0.958 * J
case1 <- 0.02 * J
case2 <- 0.02 * J
case3 <- 0.002 * J
# effect size of association
effSize <- 4

# generate data
set.seed(0)
medDat <- rbind(cbind(rnorm(n=case0), rnorm(n=case0)),
                cbind(rnorm(n=case1, mean=effSize), rnorm(n=case1)),
                cbind(rnorm(n=case2), rnorm(n=case2, mean=effSize)),
                cbind(rnorm(n=case3, mean=effSize), rnorm(n=case3, mean=effSize)))

# intial starting values
maxMeans = matrix(data=c(8,8), nrow=2)
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=runif(n=4, min=0, max=min(maxMeans)), nrow=2, ncol=2), matrix(data=runif(n=4, min=0, max=min(maxMeans)), nrow=2, ncol=2), maxMeans)
initPiList <- list(c(0.82), c(0.02, 0.02),c(0.02, 0.02), c(0.1))
# fit the model
csmGmmOutput <- symm_fit_ind_EM(testStats = medDat, initMuList = initMuList, initPiList = initPiList,
                                checkpoint=FALSE)

# rejections at q=0.1
outputDF <- data.frame(Z1 = medDat[, 1], Z2 = medDat[, 2], origIdx = 1:J,
                         lfdrValue=csmGmmOutput$lfdrResults) %>%
    arrange(lfdrValue) %>%
    mutate(lfdrAvg = cummean(lfdrValue)) %>%
    mutate(Causal = ifelse(origIdx >= J - case3 + 1, 1, 0)) %>%
    mutate(Rej = ifelse(lfdrAvg < 0.1, 1, 0))

# number of false rejections
length(which(outputDF$Causal == 0 & outputDF$Rej == 1))
# number of true rejections
length(which(outputDF$Causal == 1 & outputDF$Rej == 1))


## -----------------------------------------------------------------------------

# number of SNPs and proportion in each case
J <- 40000
K <- 2
case0 <- 0.958 * J
case1 <- 0.02 * J
case2 <- 0.02 * J
case3 <- 0.002 * J
# effect size of association
effSize <- 4

# generate data
set.seed(0)
corMat <- matrix(data=c(1, 0.3, 0.3, 1), nrow=2)
pleioDat <- rbind(mvtnorm::rmvnorm(n=case0, sigma=corMat),
                mvtnorm::rmvnorm(n=case1, mean=c(effSize, 0), sigma=corMat),
                mvtnorm::rmvnorm(n=case2, mean=c(0, effSize), sigma=corMat),
                mvtnorm::rmvnorm(n=case3, mean=c(effSize, effSize), sigma=corMat))

# estimate the correlation from data
estCor <- cor(pleioDat)[1,2]
estCorMat <- matrix(data=c(1, estCor, estCor, 1), nrow=2)

# intial starting values
maxMeans = matrix(data=c(8,8), nrow=2)
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=runif(n=2, min=0, max=min(maxMeans)), nrow=2, ncol=1), matrix(data=runif(n=2, min=0, max=min(maxMeans)), nrow=2, ncol=1), maxMeans)
initPiList <- list(c(0.82), c(0.04),c(0.04), c(0.1))
# fit the model
c_csmGmm <- symm_fit_cor_EM(testStats = pleioDat, initMuList = initMuList, initPiList = initPiList,
                            corMat = estCorMat, checkpoint=FALSE)

# rejections at q=0.1
c_outputDF <- data.frame(Z1 = pleioDat[, 1], Z2 = pleioDat[, 2], origIdx = 1:J,
                         lfdrValue=c_csmGmm$lfdrResults) %>%
    arrange(lfdrValue) %>%
    mutate(lfdrAvg = cummean(lfdrValue)) %>%
    mutate(Causal = ifelse(origIdx >= J - case3 + 1, 1, 0)) %>%
    mutate(Rej = ifelse(lfdrAvg < 0.1, 1, 0))

# number of false rejections
length(which(c_outputDF$Causal == 0 & c_outputDF$Rej == 1))
# number of true rejections
length(which(c_outputDF$Causal == 1 & c_outputDF$Rej == 1))


## -----------------------------------------------------------------------------

# number of SNPs and proportion in each case
J <- 40000
K <- 2
case0 <- 0.958 * J
case1 <- 0.02 * J
case2 <- 0.02 * J
case3 <- 0.002 * J
# effect size of association
effSize <- 4

# generate data
set.seed(0)
repDat <- rbind(cbind(rnorm(n=case0), rnorm(n=case0)),
                cbind(rnorm(n=case1/2, mean=effSize), rnorm(n=case1/2)),
                cbind(rnorm(n=case1/2, mean=-effSize), rnorm(n=case1/2)),
                cbind(rnorm(n=case2/2), rnorm(n=case2/2, mean=effSize)),
                cbind(rnorm(n=case2/2), rnorm(n=case2/2, mean=-effSize)),
                cbind(rnorm(n=case3/4, mean=-effSize), rnorm(n=case3/4, mean=effSize)),
                cbind(rnorm(n=case3/4, mean=effSize), rnorm(n=case3/4, mean=-effSize)),
                cbind(rnorm(n=case3/4, mean=effSize), rnorm(n=case3/4, mean=effSize)),
                cbind(rnorm(n=case3/4, mean=-effSize), rnorm(n=case3/4, mean=-effSize)))

# intial starting values
maxMeans = matrix(data=c(8,8), nrow=2)
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=runif(n=4, min=0, max=min(maxMeans)), nrow=2, ncol=2), matrix(data=runif(n=4, min=0, max=min(maxMeans)), nrow=2, ncol=2), maxMeans)
initPiList <- list(c(0.82), c(0.02, 0.02),c(0.02, 0.02), c(0.1))
# fit the model
r_csmGmm <- symm_fit_ind_EM(testStats = repDat, initMuList = initMuList, initPiList = initPiList,
                                sameDirAlt=TRUE, checkpoint=FALSE)

# rejections at q=0.1
r_outputDF <- data.frame(Z1 = repDat[, 1], Z2 = repDat[, 2], origIdx = 1:J,
                         lfdrValue=r_csmGmm$lfdrResults) %>%
    arrange(lfdrValue) %>%
    mutate(lfdrAvg = cummean(lfdrValue)) %>%
    mutate(Causal = ifelse(origIdx >= J - case3/2 + 1, 1, 0)) %>%
    mutate(Rej = ifelse(lfdrAvg < 0.1, 1, 0))

# number of false rejections
length(which(r_outputDF$Causal == 0 & r_outputDF$Rej == 1))
# number of true rejections
length(which(r_outputDF$Causal == 1 & r_outputDF$Rej == 1))


