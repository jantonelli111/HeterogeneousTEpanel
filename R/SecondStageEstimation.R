## Functions that take the posterior distribution of the unit
## and time specific treatment effects and summarize them into
## either marginal or heterogeneous treatment effects

DeltaQ = function(indTauPost, nTimesOut) {
  
  timeSpecificEst = apply(apply(indTauPost[,1:nTimesOut,], 2:3, mean, na.rm=TRUE), 1, mean)
  timeSpecificSD = apply(apply(indTauPost[,1:nTimesOut,], 2:3, mean, na.rm=TRUE), 1, sd)
  timeSpecificCIlower = apply(apply(indTauPost[,1:nTimesOut,], 2:3, mean, na.rm=TRUE), 1, quantile, .025)
  timeSpecificCIupper = apply(apply(indTauPost[,1:nTimesOut,], 2:3, mean, na.rm=TRUE), 1, quantile, .975)
  
  cumulativeEst = cumulativeSD = cumulativeCIlower = 
    cumulativeCIupper = rep(NA, nTimesOut)
  
  cumulativeEst[1] = timeSpecificEst[1]
  cumulativeSD[1] = timeSpecificSD[1]
  cumulativeCIlower[1] = timeSpecificCIlower[1]
  cumulativeCIupper[1] = timeSpecificCIupper[1]
  
  for (jjj in 2 : nTimesOut) {
    cumulativeEst[jjj] = mean(apply(indTauPost[,1:jjj,], c(1,3), sum, na.rm=TRUE))
    cumulativeSD[jjj] = sd(apply(apply(indTauPost[,1:jjj,],  c(1,3), sum, na.rm=TRUE), 2, mean))
    cumulativeCIlower[jjj] = quantile(apply(apply(indTauPost[,1:jjj,],  c(1,3), sum, na.rm=TRUE), 2, mean), .025)
    cumulativeCIupper[jjj] = quantile(apply(apply(indTauPost[,1:jjj,], c(1,3), sum, na.rm=TRUE), 2, mean), .975)
  }
  
  
  estimates = list(timeSpecific = timeSpecificEst, cumulative = cumulativeEst)
  SD = list(timeSpecific = timeSpecificSD, cumulative = cumulativeSD)
  CIlower = list(timeSpecific = timeSpecificCIlower, cumulative = cumulativeCIlower)
  CIupper = list(timeSpecific = timeSpecificCIupper, cumulative = cumulativeCIupper)
  
  timeSpecific = list(estimates = timeSpecificEst,
                      SD = timeSpecificSD,
                      CIlower = timeSpecificCIlower,
                      CIupper = timeSpecificCIupper)
  
  cumulative = list(estimates = cumulativeEst,
                    SD = cumulativeSD,
                    CIlower = cumulativeCIlower,
                    CIupper = cumulativeCIupper)
  
  l = list(timeSpecific=timeSpecific,
           cumulative=cumulative)
  
  return(l)
}

DeltaQsmooth = function(indTauPost, df=3) {
  
  nSubjects = dim(indTauPost)[1]
  nTimesOut = dim(indTauPost)[2]
  nKeep = dim(indTauPost)[3]
  
  deltaQ = matrix(NA, nKeep, nTimesOut)
  cumulativeDeltaQ = matrix(NA, nKeep, nTimesOut)
  for (ni in 1 : nKeep) {
    newdata = data.frame(id = rep(1:nSubjects, each=nTimesOut),
                         lag = rep(1:nTimesOut, nSubjects))
    
    newdata2 = newdata
    
    newdata2$tau = rep(NA, nrow(newdata2))
    for (i in 1 : nSubjects) {
      for (j in 1 : nTimesOut) {
        w = which(newdata2$id == i & newdata2$lag == j)
        newdata2$tau[w] = indTauPost[i,j,ni]
      }
    }
    
    ## remove the id variable
    newdata2$id = NULL
    
    newdataPred = newdata2
    newdataPred$tau = NULL
    
    modTau = lm(tau ~ . -lag + ns(lag, df=df), data=newdata2)
    fittedValues = predict(modTau, newdataPred)
    
    for (j in 1 : nTimesOut) {
      wj = which(newdata2$lag == j)
      deltaQ[ni,j] = mean(fittedValues[wj])
      cumulativeDeltaQ[ni,j] = sum(deltaQ[ni,1:j])
    }
  }
  
  timeSpecific = list(estimates = apply(deltaQ, 2, mean),
                      SD = apply(deltaQ, 2, sd),
                      CIlower = apply(deltaQ, 2, quantile, .025),
                      CIupper = apply(deltaQ, 2, quantile, .975))
  
  cumulative = list(estimates = apply(cumulativeDeltaQ, 2, mean),
                    SD = apply(cumulativeDeltaQ, 2, sd),
                    CIlower = apply(cumulativeDeltaQ, 2, quantile, .025),
                    CIupper = apply(cumulativeDeltaQ, 2, quantile, .975))
  
  l = list(timeSpecific=timeSpecific,
           cumulative=cumulative)
  
  return(l)
}

Linear2StageSmooth = function(indTauPost, newX, covariates, df=3) {
  nSubjects = dim(indTauPost)[1]
  nTimesOut = dim(indTauPost)[2]
  nKeep = dim(indTauPost)[3]
  
  nCovariates = dim(covariates)[2]
  
  covariatesNew = covariates
  for (j in 1 : (nTimesOut - 1)) {
    covariatesNew = rbind(covariatesNew, covariates)
  }
  
  diffXpost = matrix(NA, nKeep, nCovariates)
  betaPost = matrix(NA, nKeep, nCovariates)
  
  if (!is.null(newX)) {
    newXpred = matrix(NA, nKeep, nrow(newX))
  }
  
  for (ni in 1 : nKeep) {
    newdata = data.frame(id = rep(1:nSubjects, nTimesOut),
                         lag = rep(1:nTimesOut, each=nSubjects))
    
    newdata2 = cbind(newdata, covariatesNew)
    
    newdata2$tau = rep(NA, nrow(newdata2))
    for (i in 1 : nSubjects) {
      for (j in 1 : nTimesOut) {
        w = which(newdata2$id == i & newdata2$lag == j)
        newdata2$tau[w] = indTauPost[i,j,ni]
      }
    }
    
    ## remove the id variable
    newdata2$id = NULL
    
    newdataPred = newdata2
    newdataPred$tau = NULL
    
    modTau = lm(tau ~ . -lag + ns(lag, df=df), data=newdata2)
    fittedValues = predict(modTau, newdataPred)
    
    betaPost[ni,] = modTau$coefficients[2 : (nCovariates + 1)]
    
    if (!is.null(newX)) {
      newXpred[ni,] = predict(modTau, newX) 
    }
    
    for (jj in 1 : nCovariates) {
      newdataPred1 = newdataPred
      newdataPred0 = newdataPred
      
      newdataPred1[,jj+1] = quantile(covariates[,jj], .75)
      newdataPred0[,jj+1] = quantile(covariates[,jj], .25)
      
      diffXpost[ni,jj] = sum(predict(modTau, newdataPred1) - 
                               predict(modTau, newdataPred0)) / nSubjects 
    }
  }
  
  if (!is.null(newX)) {
    newXlist = list(estimates = apply(newXpred, 2, mean, na.rm=TRUE),
                    SD = apply(newXpred, 2, sd, na.rm=TRUE),
                    CIlower = apply(newXpred, 2, quantile, .025, na.rm=TRUE),
                    CIupper = apply(newXpred, 2, quantile, .975, na.rm=TRUE))
  } else {
    newXlist = NULL
  }
  
  diffX = apply(diffXpost, 2, mean, na.rm=TRUE)
  diffXsd = apply(diffXpost, 2, sd, na.rm=TRUE)
  diffXlower.ci = apply(diffXpost, 2, quantile, 0.025, na.rm=TRUE)
  diffXupper.ci = apply(diffXpost, 2, quantile, 0.975, na.rm=TRUE)
  
  diffXlist = list(estimates = diffX,
                   SD = diffXsd,
                   CIlower = diffXlower.ci,
                   CIupper = diffXupper.ci)
  
  betalist = list(estimates = apply(betaPost, 2, mean, na.rm=TRUE),
                  SD = apply(betaPost, 2, sd, na.rm=TRUE),
                  CIlower = apply(betaPost, 2, quantile, .025, na.rm=TRUE),
                  CIupper = apply(betaPost, 2, quantile, .975, na.rm=TRUE))
  
  l = list(beta = betalist, newX=newXlist)
  
  return(l)
}

Linear2Stage = function(indTauPost, newX, covariates) {
  nSubjects = dim(indTauPost)[1]
  nTimesOut = dim(indTauPost)[2]
  nKeep = dim(indTauPost)[3]
  
  nCovariates = dim(covariates)[2]
  
  covariatesNew = covariates
  for (j in 1 : (nTimesOut - 1)) {
    covariatesNew = rbind(covariatesNew, covariates)
  }
  
  diffXpost = matrix(NA, nKeep, nCovariates)
  betaPost = matrix(NA, nKeep, nCovariates)
  
  if (!is.null(newX)) {
    newXpred = matrix(NA, nKeep, nrow(newX))
  }
  
  for (ni in 1 : nKeep) {
    newdata = data.frame(id = rep(1:nSubjects, nTimesOut),
                         lag = rep(1:nTimesOut, each=nSubjects))
    
    newdata2 = cbind(newdata, covariatesNew)
    
    newdata2$tau = rep(NA, nrow(newdata2))
    for (i in 1 : nSubjects) {
      for (j in 1 : nTimesOut) {
        w = which(newdata2$id == i & newdata2$lag == j)
        newdata2$tau[w] = indTauPost[i,j,ni]
      }
    }
    
    ## remove the id variable
    newdata2$id = NULL
    
    newdataPred = newdata2
    newdataPred$tau = NULL
    
    modTau = lm(tau ~ . -lag + as.factor(lag), data=newdata2)
    fittedValues = predict(modTau, newdataPred)
    
    betaPost[ni,] = modTau$coefficients[2 : (nCovariates + 1)]
    
    if (!is.null(newX)) {
      newXpred[ni,] = predict(modTau, newX) 
    }
    
    for (jj in 1 : nCovariates) {
      newdataPred1 = newdataPred
      newdataPred0 = newdataPred
      
      newdataPred1[,jj+1] = quantile(covariates[,jj], .75)
      newdataPred0[,jj+1] = quantile(covariates[,jj], .25)
      
      diffXpost[ni,jj] = sum(predict(modTau, newdataPred1) - 
                               predict(modTau, newdataPred0)) / nSubjects 
    }
  }
  
  if (!is.null(newX)) {
    newXlist = list(estimates = apply(newXpred, 2, mean, na.rm=TRUE),
                    SD = apply(newXpred, 2, sd, na.rm=TRUE),
                    CIlower = apply(newXpred, 2, quantile, .025, na.rm=TRUE),
                    CIupper = apply(newXpred, 2, quantile, .975, na.rm=TRUE))
  } else {
    newXlist = NULL
  }
  
  diffX = apply(diffXpost, 2, mean, na.rm=TRUE)
  diffXsd = apply(diffXpost, 2, sd, na.rm=TRUE)
  diffXlower.ci = apply(diffXpost, 2, quantile, 0.025, na.rm=TRUE)
  diffXupper.ci = apply(diffXpost, 2, quantile, 0.975, na.rm=TRUE)
  
  diffXlist = list(estimates = diffX,
                   SD = diffXsd,
                   CIlower = diffXlower.ci,
                   CIupper = diffXupper.ci)
  
  betalist = list(estimates = apply(betaPost, 2, mean, na.rm=TRUE),
                  SD = apply(betaPost, 2, sd, na.rm=TRUE),
                  CIlower = apply(betaPost, 2, quantile, .025, na.rm=TRUE),
                  CIupper = apply(betaPost, 2, quantile, .975, na.rm=TRUE))
  
  l = list(beta = betalist, newX=newXlist)
  
  return(l)
}

