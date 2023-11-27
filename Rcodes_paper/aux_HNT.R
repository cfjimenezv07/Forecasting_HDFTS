# auxiliary functions for the HNT method

source('./MITS_class.R') 

init_data_new <-  function(data,p,Names){
  NAMES = Names
  N <- length(data)
  out = c('data','mixed_fts')
  mixed_fts <- mits$new()
  
  age_max = dim(data[[1]])[1]
  T = dim(data[[1]])[2]
  basis <- create.bspline.basis(c(0, 1), nbasis = p, norder = 4)
  args <- seq(0, 1, length=age_max)
  
  fd<-list()
  for(i in 1:N){
    fd[[i]] <- Data2fd(args, data[[i]], basis)
    mixed_fts$loadFunc(fd[[i]])
  }
  mixed_fts$actualize()
  names(fd) <- NAMES
  names(data) <- NAMES
  return(mget(out))
}

trim_data_new <- function(data_frame,T){
  out <- c('data')
  data <- list()
  names <- names(data_frame)
  for(i in 1:length(data_frame)){
    data[[names[i]]] <- data_frame[[names[i]]][,T:(T+n_year-h_max-1)]
    rownames(data[[names[i]]]) <- rownames(data_frame[[names[i]]])
    colnames(data[[names[i]]]) <- colnames(data_frame[[names[i]]])[T:(T+n_year-h_max-1)]
  }
  return(mget(out))
}


# Remove the zeros
remove_zeroes <- function(data_raw,n_prefectures,n_year,n_age) {
  N       <- n_prefectures
  T       <- n_year
  age_max <- n_age
  
  data_nozero <- list()
  for(i in 1:N) {
    data_nozero[[i]] <- matrix(data_raw[[i]],age_max,T)
  }
  for(i in 1:N) {
    for(j in 2:age_max) {#For j=1, the two zeroes have been removed before
      for(k in 1:T){ 
        if(data_nozero[[i]][j,k]==-Inf){
          data_nozero[[i]][j,k] <- data_nozero[[i]][j-1,k]
        }
      }
    }
  }
  return(data_nozero)
}
