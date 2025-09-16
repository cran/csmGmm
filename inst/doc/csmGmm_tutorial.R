## ----Mediation, cache=TRUE, warning=FALSE, message=FALSE----------------------
library(csmGmm)
library(dplyr)

# define number of SNPs and proportion in each case
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

# define intial starting values
initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                   matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
# fit the model
csmGmmOutput <- symm_fit_ind_EM(testStats = medDat, initMuList = initMuList, initPiList = initPiList,
                                checkpoint=FALSE)

# summarize output and find which rows are rejected at q=0.1
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


## ----three pleiotropy, cache=TRUE---------------------------------------------
# define number of SNPs and proportion in each case
J <- 10000
K <- 3
# proportion of SNPs with no association
case0 <- 0.959 * J
# proportion of SNPs with association in only one column (divided by three columns)
oneCol <- (0.03 / 3) * J
# proportion of SNPs with associations in two columns (divided by three possible pairs)
twoCol <- (0.01 / 3) * J
# proportions of SNPs with associations in all three columns 
threeCol <- 0.001 * J

# effect size of association
effSize <- 4

# generate data
set.seed(0)
pleioDat3D <- rbind(cbind(rnorm(n=case0 + 1), rnorm(n=case0 + 1), 
                          rnorm(n=case0 + 1)), # +1 is for rounding purposes, to get to 100k sets
                cbind(rnorm(n=oneCol, mean=effSize), rnorm(n=oneCol), rnorm(n=oneCol)),
                cbind(rnorm(n=oneCol), rnorm(n=oneCol, mean=effSize), rnorm(n=oneCol)),
                cbind(rnorm(n=oneCol), rnorm(n=oneCol), rnorm(n=oneCol, mean=effSize)),
                cbind(rnorm(n=twoCol, mean=effSize), rnorm(n=twoCol, mean=effSize), rnorm(n=twoCol)),
                cbind(rnorm(n=twoCol, mean=effSize), rnorm(n=twoCol), rnorm(n=twoCol, mean=effSize)),
                cbind(rnorm(n=twoCol), rnorm(n=twoCol, mean=effSize), rnorm(n=twoCol, mean=effSize)),
                cbind(rnorm(n=threeCol, mean=effSize), rnorm(n=threeCol, mean=effSize), rnorm(n=threeCol, mean=effSize)))
# flip half of the effect sizes
flipRows <- sample(x=1:J, size=J/2, replace=FALSE)
pleioDat3D[flipRows, ] <- -1 * pleioDat3D[flipRows, ]

# define initial starting values
initPiList3D <- list(c(0.82))
for (i in 2:7) {initPiList3D[[i]] <- c(0.08 / 12, 0.08 / 12)}
initPiList3D[[8]] <- c(0.1)
initMuList3D <- list(matrix(data=0, nrow=3, ncol=1)) # package will add the appropriate 0s to initMuList
for (i in 2:7) {
  initMuList3D[[i]] <- cbind(rep(2, 3), rep(5, 3))
}
initMuList3D[[8]] <- matrix(data=c(8, 8, 8), nrow=3)

# fit the model
csmGmmOutput3D <- symm_fit_ind_EM(testStats = pleioDat3D, initMuList = initMuList3D, 
                                  initPiList = initPiList3D, eps = 10^(-4), checkpoint=FALSE)

# summarize output and find which rows are rejected at q=0.1
outputDF3D <- data.frame(Z1 = pleioDat3D[, 1], Z2 = pleioDat3D[, 2], Z3 = pleioDat3D[, 3], origIdx = 1:J,
                         lfdrValue=csmGmmOutput3D$lfdrResults) %>%
    arrange(lfdrValue) %>%
    mutate(lfdrAvg = cummean(lfdrValue)) %>%
    mutate(Causal = ifelse(origIdx >= J - threeCol + 1, 1, 0)) %>%
    mutate(Rej = ifelse(lfdrAvg < 0.1, 1, 0))

# number of false rejections
length(which(outputDF3D$Causal == 0 & outputDF3D$Rej == 1))
# number of true rejections
length(which(outputDF3D$Causal == 1 & outputDF3D$Rej == 1))


## ----four pleiotropy, cache=TRUE----------------------------------------------

# define number of SNPs and proportion in each case
J <- 1000
K <- 4

# generate data - no signals, all rows are null
set.seed(0)
fourDimZ <- cbind(rnorm(J), rnorm(J), rnorm(J), rnorm(J))

# define initial values for higher dimensions
initPiListHD <- list(c(0.82))
for (i in 2:(2^K)) {initPiListHD[[i]] <- 0.18 / (2^K - 1)}
initMuListHD <- list(matrix(data=rep(0, K), nrow=K, ncol=1))
for (i in 2:(2^K)) {
  initMuListHD[[i]] <- matrix(data=rep(3, K), nrow=K, ncol=1)
}

# fit the model (note: we usually change the convergence criterion for higher dimensions to increase speed)
csmGmmOutputHD <- symm_fit_ind_EM(testStats = fourDimZ, initMuList = initMuListHD, initPiList = initPiListHD, eps=10^(-1), checkpoint=F)

# summarize output and find which rows are rejected at q=0.1
outputDFHD <- data.frame(Z1 = fourDimZ[, 1], Z2 = fourDimZ[, 2], Z3 = fourDimZ[, 3], Z4 = fourDimZ[, 4],
                         origIdx = 1:J,
                         lfdrValue=csmGmmOutputHD$lfdrResults) %>%
    arrange(lfdrValue) %>%
    mutate(lfdrAvg = cummean(lfdrValue)) %>%
    mutate(Rej = ifelse(lfdrAvg < 0.1, 1, 0))

# no true causal SNPs, and no rejections
head(outputDFHD)


## ----c-csmGmm, cache=T--------------------------------------------------------

# define number of SNPs and proportion in each case
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

# define intial starting values
initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                   matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
# fit the model
c_csmGmm <- symm_fit_cor_EM(testStats = pleioDat, initMuList = initMuList, initPiList = initPiList,
                            corMat = estCorMat, checkpoint=FALSE)

# summarize output and find which rows are rejected at q=0.1
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


## ----r-csmGmm, cache=TRUE-----------------------------------------------------

# define number of SNPs and proportion in each case
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

# define intial starting values
initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                   matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
# fit the model
r_csmGmm <- symm_fit_ind_EM(testStats = repDat, initMuList = initMuList, initPiList = initPiList,
                                sameDirAlt=TRUE, checkpoint=FALSE)

# summarize output and find which rows are rejected at q=0.1
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


