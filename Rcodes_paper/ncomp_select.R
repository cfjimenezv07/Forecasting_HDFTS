#################################
# select m_component and p_order
#################################

ChoosingVersionI <- function(FTS, K, dianismap, sampling_points, criterion = c("FFPE", "FPE", "AIC", "AICC", "BIC", "Shibata"))
{
    n <- length(FTS$coef[1,])
    pca <- pca.fd(FTS, nharm = K, centerfns = FALSE)
    scores <- pca$scores #h j grami antiprosopeui to \xi_(1,j),\xi_(2,j)...\xi_(m,j),\xi_{j,t}=<X_t,v_j>
    ###################################Step A####################################
    
    #####################calculate the ratio#########################
    pca1 <- pca.fd(FTS, nharm = K, centerfns = FALSE)
    scores1 <- pca1$scores
    ################################Step A#################################
    jomega <- jsizigis <- matrix(NA, K, n)
    for(i in 1:K)
    {
        jomega[i,] <- fft(t(scores1)[i,])
        jsizigis[i,] <- fft(t(scores1)[i,], inverse = TRUE)
    }
    periodograma <- list()
    for(vs in 1:length(jomega[1,]))
    {
        period <- (1/(2*pi*n)) * (jomega[,vs] %*% t(jsizigis[,vs]))
        periodograma[vs] <- list(period)
    }
    
    ##############################STEP B####################################
    metrogiaq <- function(j1, j2)
    {
        athrisma <- 0
        for(vs in 2:(floor(n/2)+1))
        {
            athrisma <- athrisma + (Mod(periodograma[[vs]][j1,j2]))^2 + (Mod(Conj(t(periodograma[[vs]]))[j1,j2]))^2
        }
        athrisma <- ((2*pi)/n) * athrisma
        result <- list(athrismametroug1 = athrisma)
    }
    
    pinakasq <- matrix(NA, K, K)
    for(v in 1:K)
    {
        for(t in 1:K)
        {
          pinakasq[v,t] <- metrogiaq(v, t)$athrismametroug1
        }
    }
    
    ##############################STEP C####################################
    t = sampling_points
    #h entoli fft ipologizi gia j=0,1,...,n-1(ara thesi 1 antistixi se j<-0)
    jomega <- matrix(NA, length(sampling_points), n)
    jomegaInverse <- matrix(NA, length(sampling_points), n)
    for(el in 1:(length(sampling_points)))
    {
        jomega[el,] <- fft(eval.fd(t, FTS)[el,])
        jomegaInverse[el,] <- fft(eval.fd(t, FTS)[el,], inverse = TRUE)
    }
    periodogramaStaXi <- list()
    #to periodograma gia arnitika j einai to anastrofo gia thetika j
    for(ts in 1:length(jomega[1,]))
    {
        periodStaXi <- (1/(2*pi*n)) * (jomega[,ts] %*% t(jomegaInverse[,ts]))
        periodStaXi <- mean((Mod(periodStaXi))^2)
        periodogramaStaXi[ts] <- list(periodStaXi)
    }
    #stin thesi 1 tou periodogramaStaXi(periodograma[[1]]) antistixi i timi j<-0
    
    paronomastis <- 0
    for(vs in 2:(floor(n/2)+1))
    {
        paronomastis <- paronomastis + 2 * periodogramaStaXi[[vs]]
    }
    paronomastis <- ((2*pi)/n)*paronomastis
    
    arxikom <- 1
    while((sum(apply(as.matrix(pinakasq[1:arxikom,1:arxikom]),1,sum))+(1/(2*pi))*sum((pca1$values[(arxikom+1):K])^2))/(paronomastis) < 0.95)
    {
        arxikom<-arxikom+1
    }
    telikom <- arxikom
    
    cn <- sqrt(n)/log10(n)
    arxikomperispomeni <- 1
    while((pca$values[1]/pca$values[arxikomperispomeni]) < cn)
    {
        arxikomperispomeni<-arxikomperispomeni+1
    }
    telikomperispomeni <- arxikomperispomeni-1
    
    ############################################
    teldianFrobenius <- rep(NA,length(dianismap))
    for(i in 1:length(dianismap))
    {
        p <- dianismap[i]
        ################################Step A-B#################################
        #Find the matrix \Sigma_e
        ################Step 1#################
        m = telikom
        #################Step 2#################
        pca1 <- pca.fd(FTS, nharm = m, centerfns = FALSE) #if centerfns = TRUE the scores are computed for the centered functions X_i-\bar{X_n}(vivlio sel34)
        scores1 <- pca1$scores   #h j grami antiprosopeui to \xi_(1,j),\xi_(2,j)...\xi_(m,j),\xi_{j,t}=<X_t,v_j>
        #################Step 3#################
        #https://en.wikipedia.org/wiki/Vector_autoregression
        #http://www-stat.wharton.upenn.edu/~steele/Courses/956/ResourceDetails/YWSourceFiles/YW-Eshel.pdf
        ektimitriaYW <- ar(scores1, method = "yule-walker", order.max = p, aic = FALSE, demean = FALSE)
        dianismae <- matrix(ektimitriaYW$resid, ncol = m)
        pinakasSigma <- t(dianismae[(p+1):n,]) %*% dianismae[(p+1):n,]
        
        if(criterion == "FFPE")
        {
            Frobenius<-(1/(n^2))*sum(apply(pinakasSigma^2,1,sum))
            teldianFrobenius[i]<-(((n+p*m+1)/(n-p*m-1))^2)*Frobenius
        }
        if(criterion == "FPE")
        {
            Frobenius<-(1/n)*pinakasSigma
            teldianFrobenius[i]<-(((n+p*m+1)/(n-p*m-1))^m)*(det(Frobenius))
        }
        if(criterion == "AIC")
        {
            Frobenius<-(1/n)*pinakasSigma
            teldianFrobenius[i]<-log(det(Frobenius))+(2*(m^2)*i)/n
        }
        if(criterion == "AICC")
        {
            Frobenius <- (1/n) * pinakasSigma
            teldianFrobenius[i] <- log(det(Frobenius)) + ((m^2)*i + n * m)/(n - i * m - m - 1)
        }
        if(criterion == "BIC")
        {
            Frobenius<-(1/n)*pinakasSigma
            teldianFrobenius[i]<-log(det(Frobenius))+(log(n)*(m^2)*i)/n
        }
        if(criterion == "Shibata")
        {
            Frobenius<-(1/n)*pinakasSigma
            teldianFrobenius[i]<-det(Frobenius)*(1+2*((m*i+1)/n))^m
        }
    }
    return(list(teldianFrobenius = teldianFrobenius, telikom = telikom,
                   telikomperispomeni = telikomperispomeni))
}
