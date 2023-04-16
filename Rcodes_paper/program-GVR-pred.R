#############
library(fda)
#############

############################################################################
#h sinartisi fma1 dimiourga functional moving average (1)
# x_t(ti)=\mu(\ti)+\epsilon_t(ti)+\theta*\epsilon_{t-1}(ti)
############################################################################


# Function simulates functional moving average(1) time serie
# Input:
#   n - number of observations
#   Theta1 - 1st moving average operator
#   sd - noise standard deviation
#   basis - fda basis
#   first - coefs of the first element
#   mesitimi- \mu(\ti)

K = L = 21
basis = create.fourier.basis(rangeval = c(0, 1), nbasis = L)
n = 300

#Arxika ipothetoume oti \mu(\ti)=0

fma1 = function(n, basis = NULL, Theta1 = NULL, Sigma = NULL, sampling_points = seq(0, 1, length.out = 21))
{
	  # if no basis than create 21 Fourier harmonics
	  if(is.null(basis))
		{
	      basis = create.fourier.basis(rangeval = c(0, 1), nbasis = length(sampling_points))
	  }

  	# get number of coefficients
	  D = basis$nbasis

  	# if no operator is provided, Theta1 und Theta2 are specified as 0.8*identity
	  if(is.null(Theta1))
		{
	      Theta1 = 0.8 * diag(D)
	  }

  	# if no Sigma is specified then use Sigma=1/seq(1:D)
	  if(is.null(Sigma))
		{
	      Sigma = 1/seq(1:D)
	  }
	  # build coefficients matrix (initially null)
	  coef = matrix(data = 0, nrow = D, ncol = n)

	  # define the MA1
    zlag0 = matrix(0,D,n)
	  for(i in 1:n)
	  {
	      zlag0[,i] = rnorm(n = D, mean = 0, sd = Sigma)
	  }
	  zlag1 = matrix(0,D,n)
	  for(i in 2:n)
	  {
	      zlag1[,i] = Theta1 %*% zlag0[,i-1]
	  }
    
	  coef = zlag1 + zlag0
 	
	  # return time series
	  Fdata <- fd(coef, basis = basis)
    meanfdata <- mean.fd(Fdata)

    fdvals <- eval.fd(sampling_points, meanfdata)
    result <- list(fdvals = fdvals, Fdata = Fdata)
}

############################################################################

#FTS<-fma1(n,basis)$Fdata

############################################################################
#
# Input
# FTS - H xronosira X_t
# m - number of functional principal components (=d)
# p - order of autoregressive time series
# K - number of functional principal components
############################################################################


ChoosingVersionI <- function(FTS, K, dianismap, sampling_points = seq(0, 1, length.out = 21),
                             criterion = c("FFPE", "FPE", "AIC", "AICC", "BIC", "Shibata"))
{
    n <- length(FTS$coef[1,])
    pca <- pca.fd(FTS, nharm = K, centerfns = FALSE)
    scores <- pca$scores 
    
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
    while((sum(apply(as.matrix(pinakasq[1:arxikom,1:arxikom]),1,sum))+(1/(2*pi))*sum((pca1$values[(arxikom+1):K])^2))/(paronomastis) < 0.85)
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

        ##...FFPE  OURS .....
        if(criterion == "FFPE")
        {
            Frobenius <- (1/(n^2)) * sum(apply(pinakasSigma^2, 1, sum))
            teldianFrobenius[i] <- (((n + p * m + 1)/(n - p * m - 1))^2)*Frobenius
        }
        if(criterion == "FPE")
        {
            Frobenius <- (1/n) * pinakasSigma
            teldianFrobenius[i] <- (((n+p*m+1)/(n-p*m-1))^m)*(det(Frobenius))
        }
        if(criterion == "AIC")
        {
            Frobenius <- (1/n)*pinakasSigma
            teldianFrobenius[i] <- log(det(Frobenius))+(2*(m^2)*i)/n
        }
        if(criterion == "AICC")
        {
            Frobenius <- (1/n) * pinakasSigma
            teldianFrobenius[i] <- log(det(Frobenius)) + ((m^2)*i + n * m)/(n - i * m - m - 1)
        }
        if(criterion == "BIC")
        {
            Frobenius <- (1/n) * pinakasSigma
            teldianFrobenius[i] <- log(det(Frobenius)) + (log(n) * (m^2) * i)/n
        }
        if(criterion == "Shibata")
        {
            Frobenius <- (1/n) * pinakasSigma
            teldianFrobenius[i] <- det(Frobenius) * (1 + 2 * ((m * i + 1)/n))^m
        }
    }
    result <- list(teldianFrobenius = teldianFrobenius, telikom = telikom,
                   telikomperispomeni = telikomperispomeni)
}


###########################################################
katametrisi <- matrix(0, 12, 8)
#porderselect<-1:12
#mdimselect<-1:10
dianismamperispomeni <- rep(0, 10)
dianismam85 <- rep(0, 10)
dianismaMAXm <- rep(0, 10)

#Apotelesmta
R<-1000
for(r in 1:R)
{
    print(r)
    dianismap<-1:12
    #dianismad<-1:8 #m=d
    FTS <- fma1(n,basis)$Fdata
    sinartisi <- ChoosingVersionI(FTS,K,dianismap)
    apotelesma <- Re(sinartisi$teldianFrobenius)
    telikom <- sinartisi$telikom
    telikomperispomeni <- sinartisi$telikomperispomeni
    megistom <- max(telikom,telikomperispomeni)
    dianismamperispomeni[telikomperispomeni] <- dianismamperispomeni[telikomperispomeni] + 1
    dianismam85[telikom] <- dianismam85[telikom]+1
    dianismaMAXm[megistom] <- dianismaMAXm[megistom]+1

    inds = which(apotelesma == min(apotelesma), arr.ind = TRUE)
    print(paste("p =",inds,"m =",telikom))
    katametrisi[inds,telikom] <- katametrisi[inds,telikom] + 1
}

rbind(dianismam85, dianismamperispomeni, dianismaMAXm)
porderselect <- rowSums(katametrisi)
mdimselect <- colSums(katametrisi)
print(katametrisi)
print(porderselect)
print(mdimselect)

