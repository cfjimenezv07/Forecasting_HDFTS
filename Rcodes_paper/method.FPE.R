FPE.trace <- ftsa:::FPE.trace

method.FPE <- function (object, D = 21, var_type = "none", Pmax) 
{
    eigen_values = eigen(var(t(object$y)))$values
    vartot = sum(eigen_values)
    pca.FTS = ftsm(object, order = D)
    pca.FTS_scores = pca.FTS$coeff[, 2:(D + 1)]
    values = matrix(NA, D, Pmax)
    for (d in 1:D) {
      scores = pca.FTS_scores[, 1:d]
      var.explain = sum(eigen_values[1:d])
      for (p in 1:Pmax) {
        if (d == 1) {
          res = as.matrix(arima(scores, order = c(p, 0, 0))$residuals)
        }
        else {
          if (p == 0) {
            mean = t(matrix(rep(colMeans(scores), nrow(scores)), 
                            d))
            res = scores - mean
          }
          else {
            colnames(scores) <- as.character(seq(1:d))
            res = resid(vars::VAR(scores, p = p, type = var_type))
          }
        }
        values[d, p] = FPE.trace(res = res, p = p) + 
          vartot - var.explain
      }
    }
    
    
    
    rownames(values) = 1:D
    colnames(values) = 1:Pmax
    hat.p = (which.min(values))%/%D + 1
    hat.d = which.min(values)%%D
    if (hat.d == 0) 
      hat.d = hat.d + D
    out = c(hat.p, hat.d)
    return(out)
}
