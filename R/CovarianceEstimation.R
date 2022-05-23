## Initial estimate of covariance matrix

CovEstStructural = function(covariates, treatments, outcomes, zeroMat) {
  
  startTimes = apply(treatments, 1, function(x) return(which(x == 1)[1]))
  lastTime = min(startTimes) - 1
  outcomes = outcomes[,1:lastTime]
  
  N = nrow(outcomes)
  TT = ncol(outcomes)
  
  outcomes = as.matrix(outcomes)
  
  
  ## Now find the fitted values. First the mean component
  fittedValues = matrix(0, N, TT)
  for (ii in 1 : nrow(outcomes)) {
    modData = data.frame(y=outcomes[ii,], time=1:TT)
    sFit = lm(y ~ splines::ns(time, df=5), data=modData)
    fittedValues[ii,] = sFit$fitted.values
  }
  
  ## then estimate the covariance matrix
  SSCP.E <- crossprod(t(as.matrix(outcomes)) - t(fittedValues))
  SigmaHat <- SSCP.E / (TT - 6)
  
  ## set rho to be very small so that there is no penalty on the magnitude
  ## of terms, and only the desired sparsity is fixed
  test = glasso::glasso(SigmaHat, rho=0.00000001, zero=zeroMat)
  
  return(test$w)
  
}
