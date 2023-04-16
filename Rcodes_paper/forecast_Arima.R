getxreg = forecast:::getxreg
arima.string = forecast:::arima.string
fitted.Arima = forecast:::fitted.Arima
future_msts = forecast:::future_msts
copy_msts = forecast:::copy_msts
residuals.Arima = forecast:::residuals.Arima

forecast_Arima <- function (object, h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 
                             10), level = c(80, 95), fan = FALSE, xreg = NULL, lambda = object$lambda, 
          bootstrap = FALSE, npaths = 5000, biasadj = NULL, ...) 
{
  all.args <- names(formals())
  user.args <- names(match.call())[-1L]
  check <- user.args %in% all.args
  if (!all(check)) {
    error.args <- user.args[!check]
    warning(sprintf("The non-existent %s arguments will be ignored.", 
                    error.args))
  }
  use.drift <- is.element("drift", names(object$coef))
  x <- object$x <- forecast::getResponse(object)
  usexreg <- (use.drift | is.element("xreg", names(object)))
  if (!is.null(xreg) && usexreg) {
    if (!is.numeric(xreg)) 
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg)
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- if (ncol(xreg) == 1) 
        "xreg"
      else paste("xreg", 1:ncol(xreg), sep = "")
    }
    origxreg <- xreg <- as.matrix(xreg)
    h <- nrow(xreg)
  }
  else {
    if (!is.null(xreg)) {
      warning("xreg not required by this model, ignoring the provided regressors")
      xreg <- NULL
    }
    origxreg <- NULL
  }
  if (fan) {
    level <- seq(51, 99, by = 3)
  }
  else {
    if (min(level) > 0 & max(level) < 1) {
      level <- 100 * level
    }
    else if (min(level) < 0 | max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }
  level <- sort(level)
  if (use.drift) {
    n <- length(x)
    if (!is.null(xreg)) {
      xreg <- `colnames<-`(cbind(drift = (1:h) + n, xreg), 
                           make.unique(c("drift", if (is.null(colnames(xreg)) && 
                                                      !is.null(xreg)) rep("", NCOL(xreg)) else colnames(xreg))))
    }
    else {
      xreg <- `colnames<-`(as.matrix((1:h) + n), "drift")
    }
  }
  if (!is.null(object$constant)) {
    if (object$constant) {
      pred <- list(pred = rep(x[1], h), se = rep(0, h))
    }
    else {
      stop("Strange value of object$constant")
    }
  }
  else if (usexreg) {
    if (is.null(xreg)) {
      stop("No regressors provided")
    }
    object$call$xreg <- getxreg(object)
    if (NCOL(xreg) != NCOL(object$call$xreg)) {
      stop("Number of regressors does not match fitted model")
    }
    if (!identical(colnames(xreg), colnames(object$call$xreg))) {
      warning("xreg contains different column names from the xreg used in training. Please check that the regressors are in the same order.")
    }
    pred <- predict(object, n.ahead = h, newxreg = xreg)
  }
  else {
    pred <- predict(object, n.ahead = h)
  }
  if (!is.null(x)) {
    tspx <- tsp(x)
    nx <- max(which(!is.na(x)))
    if (nx != length(x) | is.null(tsp(pred$pred)) | is.null(tsp(pred$se))) {
      tspx[2] <- time(x)[nx]
      start.f <- tspx[2] + 1/tspx[3]
      pred$pred <- ts(pred$pred, frequency = tspx[3], start = start.f)
      pred$se <- ts(pred$se, frequency = tspx[3], start = start.f)
    }
  }
  nint <- length(level)
  if (bootstrap) {
    sim <- matrix(NA, nrow = npaths, ncol = h)
    for (i in 1:npaths) sim[i, ] <- simulate(object, nsim = h, 
                                             bootstrap = TRUE, xreg = origxreg, lambda = lambda)
    lower <- apply(sim, 2, quantile, 0.5 - level/200, type = 8)
    upper <- apply(sim, 2, quantile, 0.5 + level/200, type = 8)
    if (nint > 1L) {
      lower <- t(lower)
      upper <- t(upper)
    }
    else {
      lower <- matrix(lower, ncol = 1)
      upper <- matrix(upper, ncol = 1)
    }
  }
  else {
    lower <- matrix(NA, ncol = nint, nrow = length(pred$pred))
    upper <- lower
    for (i in 1:nint) {
      qq <- qnorm(0.5 * (1 + level[i]/100))
      lower[, i] <- pred$pred - qq * pred$se
      upper[, i] <- pred$pred + qq * pred$se
    }
    if (!is.finite(max(upper))) {
      warning("Upper prediction intervals are not finite.")
    }
  }
  colnames(lower) <- colnames(upper) <- paste(level, "%", sep = "")
  lower <- ts(lower)
  upper <- ts(upper)
  tsp(lower) <- tsp(upper) <- tsp(pred$pred)
  method <- arima.string(object, padding = FALSE)
  seriesname <- if (!is.null(object$series)) {
    object$series
  }
  else if (!is.null(object$call$x)) {
    object$call$x
  }
  else {
    object$call$y
  }
  fits <- fitted.Arima(object)
  if (!is.null(lambda) & is.null(object$constant)) {
    pred$pred <- InvBoxCox(pred$pred, lambda, biasadj, pred$se^2)
    if (!bootstrap) {
      lower <- InvBoxCox(lower, lambda)
      upper <- InvBoxCox(upper, lambda)
    }
  }
  return(structure(list(method = method, model = object, level = level, 
                        mean = future_msts(x, pred$pred), sim = sim, lower = future_msts(x, 
                                                                              lower), upper = future_msts(x, upper), x = x, series = seriesname, 
                        fitted = copy_msts(x, fits), residuals = copy_msts(x, 
                                                                           residuals.Arima(object))), class = "forecast"))
}
