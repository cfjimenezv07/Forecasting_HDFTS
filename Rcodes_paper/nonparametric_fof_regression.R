###########################################
# function-on-function regression with GCV
###########################################

ffunopare.knn.gcv <- function(RESPONSES, CURVES, PRED, ..., Knearest = NULL,
                              kind.of.kernel = "quadratic", semimetric = "deriv")
{
    if(is.vector(RESPONSES)) RESPONSES <- as.matrix(RESPONSES)
    if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
    testfordim <- sum(dim(CURVES)==dim(PRED))==2
    twodatasets <- T
    if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
    sm <- get(paste("semimetric.", semimetric, sep = ""))
    if(semimetric == "mplsr") stop("semimetric option mlpsr not allowed!")
    SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
    kernel <- get(kind.of.kernel)
    n1 <- ncol(SEMIMETRIC1)
    if(is.null(Knearest))
    {
        step <- ceiling(n1/100)
        if(step == 0) step <- 1
        Knearest <- seq(from = 5, to = n1 %/% 2, by = step)	
    }
    kmax <- max(Knearest)	
    p <- ncol(CURVES)
    RESPONSES.estimated <- matrix(0, nrow = n1, ncol = p)
    Bandwidth.opt <- 0
    HAT.RESP <- matrix(0, nrow = n1, ncol = length(Knearest))
    BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
    lKnearest <- length(Knearest)
    HAT.RESP <- array(0,c(n1,lKnearest,p))
    for(i in 1:n1) 
    {
        Norm.diff <- SEMIMETRIC1[, i]	
        Norm.order <- order(Norm.diff)	
        zz <- sort(Norm.diff)[2:(kmax + 2)]	
        BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
        z <- zz[ - (kmax + 1)]
        ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
        UMAT <- ZMAT/BANDWIDTH[i,  ]
        KNUM <- kernel(UMAT)
        KNUM[col(KNUM) > row(KNUM)] <- 0
        Kdenom <- apply(KNUM[Knearest,  ], 1, sum)
        WEIGHTS <- KNUM[Knearest,  ]/Kdenom
        Ind.curves <- Norm.order[2:(kmax + 1)]
        HAT.RESP[i,,] <- WEIGHTS %*% RESPONSES[Ind.curves,]
    }
    CRITARR <- array(0,c(n1,p,lKnearest))
    for(i in 1:n1)
    {
        CRITARR[i,,] <- (t(HAT.RESP[i,,]) - RESPONSES[i,])^2
    }
    Criterium <- apply(CRITARR, 3, sum)
    index.opt <- order(Criterium)[1]
    RESPONSES.estimated <- HAT.RESP[, index.opt,]
    knearest.opt <- Knearest[index.opt]
    Bandwidth.opt <- BANDWIDTH[, knearest.opt]
    Cv.estimated <- sum((RESPONSES.estimated - RESPONSES)^2)/(n1*p)
    if(twodatasets) 
    {
        SEMIMETRIC2 <- sm(CURVES, PRED, ...)
        Bandwidth2 <- 0
        n2 <- ncol(SEMIMETRIC2)
        for(k in 1:n2) 
        {
            Sm2k <- SEMIMETRIC2[, k]
            Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
        }
        KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
        KERNEL[KERNEL < 0] <- 0
        KERNEL[KERNEL > 1] <- 0
        Denom <- apply(KERNEL, 2, sum)
        NUM <- t(KERNEL) %*% RESPONSES
        RESPONSES.predicted <- NUM/Denom
        return(list(Estimated.values = RESPONSES.estimated, 
                    Predicted.values = RESPONSES.predicted, Bandwidths = 
                      Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                      Cv.estimated))
    }
    else 
    {
        return(list(Estimated.values = RESPONSES.estimated, Bandwidths
                    = Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                      Cv.estimated))
    }
}

############################
# Quadratic kernel function
############################

quadratic <- function(u)
{
    return(1 - (u)^2)
}

###########################
# semi-metric based on PCA
###########################

semimetric.pca <- function(DATA1, DATA2, q)
{
    if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
    if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
    testfordim <- sum(dim(DATA1)==dim(DATA2))==2
    twodatasets <- T
    if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
    qmax <- ncol(DATA1)
    if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
    n <- nrow(DATA1)
    COVARIANCE <- t(DATA1) %*% DATA1/n
    EIGENVECTORS <- eigen(COVARIANCE, sym = T)$vectors[, 1:q]
    COMPONENT1 <- DATA1 %*% EIGENVECTORS
    if(twodatasets) 
    {
        COMPONENT2 <- DATA2 %*% EIGENVECTORS
    }
    else 
    {
        COMPONENT2 <- COMPONENT1
    }
    SEMIMETRIC <- 0
    for(qq in 1:q)
        SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, qq], COMPONENT2[, qq], "-")^2
    return(sqrt(SEMIMETRIC))
}

