## Estimate multivariate time series model and find posterior predictive
## distribution of unknown potential outcomes

JointEstStructural = function(covariates, treatments, outcomes,
                              cov_sigma,
                              nTimesOut = 10,
                              nScans = 1500, 
                              nBurn=500, 
                              thin=4,
                              priorA = 0.33,
                              priorB = 0.0000033) {
  
  nKeep = floor((nScans - nBurn) / thin)
  keep = (1:nKeep)*thin + nBurn
  
  outcomes = as.matrix(outcomes)
  
  ## center them to have sd 1 and first point equal to zero
  for (ii in 1 : nrow(outcomes)) {
    outcomes[ii,] = outcomes[ii,] - outcomes[ii,1]
  }
  
  chSig = chol(cov_sigma)
  sigmaInv = solve(cov_sigma)
  
  startTimes = apply(treatments, 1, function(x) return(which(x == 1)[1]))
  
  N = nrow(outcomes)
  TT = ncol(outcomes)
  
  ## Store MCMC samples
  muPost = deltaPost = array(NA, dim=c(nKeep,N, TT))
  sigmaUPost = sigmaDPost = array(NA, dim=c(nKeep, N))
  
  muPostTemp = deltaPostTemp = array(NA, dim=c(N, TT))
  sigmaUPostTemp = sigmaDPostTemp = rep(NA, N)
  
  
  deltaPriorMeanZero = rep(0, N)
  
  tempVarU = tempVarD = rep(NA, N)
  
  for (i in 1 : nrow(outcomes)) {
    tempMod = lm(outcomes[i,1:(startTimes[i]-1)] ~ ns(1:(startTimes[i]-1), 6))
    
    tempVarD[i] = mean(diff(diff(tempMod$fitted.values))^2)
    tempVec = c()
    for (tt in 2 : (startTimes[i]-2)) {
      d1 = tempMod$fitted.values[tt] - tempMod$fitted.values[tt-1]
      d2 = tempMod$fitted.values[tt+1] - tempMod$fitted.values[tt]
      
      tempVec = c(tempVec, (d2-d1)^2)
    }
    tempVarU[i] = mean(tempVec)
    
    muPostTemp[i,1:(startTimes[i]-1)] = tempMod$fitted.values
    deltaPriorMeanZero[i] = muPostTemp[i,2] - muPostTemp[i,1]
    
    deltaPostTemp[i,1:(startTimes[i]-1)] = 
      c(diff(muPostTemp[i,1:(startTimes[i]-1)]), 
        muPostTemp[i,startTimes[i]-1] - muPostTemp[i,startTimes[i]-2])
  }
  
  sigmaUPostTemp[1:N] = tempVarU
  sigmaDPostTemp[1:N] = tempVarD
  
  for (ni in 2 : nScans) {
    
    Du = diag(sigmaUPostTemp)
    Dd = diag(sigmaDPostTemp)
    DuInv = diag(1/sigmaUPostTemp)
    DdInv = diag(1/sigmaDPostTemp)
    
    if (ni %% 100 == 0) print(paste(ni, "MCMC scans have finished"))
    
    ## sample first time point
    tempVar = solve(sigmaInv + 2*DuInv)
    mu2star = muPostTemp[,2] - deltaPostTemp[,1]
    tempMu = tempVar %*% t(t(outcomes[,1]) %*% sigmaInv + t(mu2star) %*% DuInv)
    
    muPostTemp[,1] = tempMu + chol(tempVar) %*% rnorm(N)
    
    tempVar = solve(DuInv + 2*DdInv)
    mu2star = muPostTemp[,2] - muPostTemp[,1]
    tempMu = tempVar %*% t(t(deltaPostTemp[,2]) %*% DdInv + 
                             t(mu2star) %*% DuInv +
                             t(deltaPriorMeanZero) %*% DdInv)
    
    deltaPostTemp[,1] = tempMu + chol(tempVar) %*% rnorm(N)
    
    ## Sample middle time points
    for (tt in 2 : (TT - 1)) {
      
      wPrevious = which(startTimes <= tt)
      wNext = which(startTimes == (tt + 1))
      wLater = which(startTimes > (tt + 1))
      
      if (length(wPrevious) == 0) {
        
        ## Update mu
        if (length(wNext) != 0) {
          DuInvZero = DuInv
          DuInvZero[wNext,wNext] = 0
          tempVar = solve(sigmaInv + DuInv + DuInvZero)
          
          
          mu2star1 = muPostTemp[,tt-1] + deltaPostTemp[,tt-1]
          mu2star2 = rep(0, N)
          
          if (length(wNext) != N) {
            mu2star2[-wNext] = muPostTemp[-wNext,tt+1] - 
              deltaPostTemp[-wNext,tt]          
          }
          
          tempMu = tempVar %*% t(t(outcomes[,tt]) %*% sigmaInv + 
                                   t(mu2star1) %*% DuInv + t(mu2star2) %*% DuInvZero)
          
          muPostTemp[,tt] = tempMu + chol(tempVar) %*% rnorm(N)
        } else {
          tempVar = solve(sigmaInv + 2*DuInv)
          
          mu2star1 = muPostTemp[,tt-1] + deltaPostTemp[,tt-1]
          mu2star2 = muPostTemp[,tt+1] - deltaPostTemp[,tt]
          tempMu = tempVar %*% t(t(outcomes[,tt]) %*% sigmaInv + 
                                   t(mu2star1) %*% DuInv + t(mu2star2) %*% DuInv)
          
          muPostTemp[,tt] = (tempMu + chol(tempVar) %*% rnorm(N))
        }
        
        
        ## Update delta
        if (length(wNext) != 0) {
          deltaPostTemp[wNext,tt] = 
            rnorm(length(wNext), mean=deltaPostTemp[wNext,tt-1],
                  sd=sqrt(sigmaDPostTemp[wNext]))
        }
        
        if (length(wLater) != 0) {
          tempVar = solve(DuInv[wLater,wLater] + 2*DdInv[wLater,wLater])
          mu2star1 = muPostTemp[wLater,tt+1] - muPostTemp[wLater,tt]
          tempMu = tempVar %*% t(t(deltaPostTemp[wLater,tt+1]) %*% DdInv[wLater,wLater] + 
                                   t(deltaPostTemp[wLater,tt-1]) %*% DdInv[wLater,wLater] +
                                   t(mu2star1) %*% DuInv[wLater,wLater])
          
          deltaPostTemp[wLater,tt] = tempMu + chol(tempVar) %*% rnorm(length(wLater)) 
        }
        
      } else if (length(wPrevious) < N & length(wPrevious) > 0) {
        
        
        
        
        
        ## Update Mu
        if (length(wNext) != 0) {
          wAll = wNext
          if (length(wLater) != 0) wAll = c(wNext, wLater)
          
          
          if (length(wNext) == 1 & length(wLater) == 0) {
            tempVar = 1/((1/cov_sigma[wAll, wAll]) + 
                           DuInv[wAll, wAll])
            
            mu2star1 = muPostTemp[wAll,tt-1] + deltaPostTemp[wAll,tt-1]
            
            tempMu = tempVar * ((outcomes[wAll,tt] / cov_sigma[wAll, wAll]) + 
                                  (mu2star1 * DuInv[wAll, wAll]))
            
            muPostTemp[wAll,tt] = tempMu + chol(tempVar) %*% rnorm(length(wAll))
            
          } else {
            DuInvZero = DuInv[wAll, wAll]
            DuInvZero[1:length(wNext),1:length(wNext)] = 0
            
            tempVar = solve(solve(cov_sigma[wAll, wAll]) + 
                              DuInv[wAll, wAll] + DuInvZero)
            
            mu2star1 = muPostTemp[wAll,tt-1] + deltaPostTemp[wAll,tt-1]
            
            mu2star2 = rep(0, length(wAll))
            
            if (length(wLater) != 0) {
              mu2star2[(length(wNext)+1) : length(wAll)] = 
                muPostTemp[wLater,tt+1] - deltaPostTemp[wLater,tt]         
            }
            
            tempMu = tempVar %*% t(t(outcomes[wAll,tt]) %*% solve(cov_sigma[wAll, wAll]) + 
                                     t(mu2star1) %*% DuInv[wAll, wAll] +
                                     t(mu2star2) %*% DuInvZero)
            
            muPostTemp[wAll,tt] = tempMu + chol(tempVar) %*% rnorm(length(wAll)) 
          }
          
        } else {
          wAll = wLater
          tempVar = solve(solve(cov_sigma[wAll, wAll]) + 2*DuInv[wAll, wAll])
          
          mu2star1 = muPostTemp[wAll,tt-1] + deltaPostTemp[wAll,tt-1]
          mu2star2 = muPostTemp[wAll,tt+1] - deltaPostTemp[wAll,tt]
          tempMu = tempVar %*% t(t(outcomes[wAll,tt]) %*% solve(cov_sigma[wAll, wAll]) + 
                                   t(mu2star1) %*% DuInv[wAll, wAll] + 
                                   t(mu2star2) %*% DuInv[wAll, wAll])
          
          muPostTemp[wAll,tt] = (tempMu + chol(tempVar) %*% rnorm(length(wAll)))
        }
        
        
        ## Update delta
        if (length(wNext) != 0) {
          deltaPostTemp[wNext,tt] = rnorm(length(wNext), 
                                          mean=deltaPostTemp[wNext,tt-1],
                                          sd=sqrt(sigmaDPostTemp[wNext]))
        }
        
        if (length(wLater) != 0) {
          tempVar = solve(DuInv[wLater,wLater] + 2*DdInv[wLater,wLater])
          mu2star1 = muPostTemp[wLater,tt+1] - muPostTemp[wLater,tt]
          tempMu = tempVar %*% t(t(deltaPostTemp[wLater,tt+1]) %*% DdInv[wLater,wLater] + 
                                   t(deltaPostTemp[wLater,tt-1]) %*% DdInv[wLater,wLater] +
                                   t(mu2star1) %*% DuInv[wLater,wLater])
          
          deltaPostTemp[wLater,tt] = tempMu + chol(tempVar) %*% rnorm(length(wLater)) 
        }
      }
      
    }
    
    ## Update variances
    for (ii in 1 : N) {
      aStar = priorA + (startTimes[ii]-1)/2
      bStar = priorB + (muPostTemp[ii,1]^2)/2 + 
        sum((diff(muPostTemp[ii,1:(startTimes[ii]-1)]) - 
               deltaPostTemp[ii,1:(startTimes[ii]-2)])^2)/2
      sigmaUPostTemp[ii] = 1/rgamma(1, aStar, bStar)
      
      aStar = priorA + (startTimes[ii]-1)/2
      bStar = priorB + (deltaPostTemp[ii,1]^2)/2 + 
        sum(diff(deltaPostTemp[ii,1:(startTimes[ii]-1)])^2)/2
      sigmaDPostTemp[ii] = 1/rgamma(1, aStar, bStar)
      
    }
    
    if (ni %in% keep) {
      tempNI = which(keep == ni)
      sigmaDPost[tempNI,] = sigmaDPostTemp
      sigmaUPost[tempNI,] = sigmaUPostTemp
      deltaPost[tempNI,,] = deltaPostTemp
      muPost[tempNI,,] = muPostTemp
    }
    
  }
  
  
  ## Forecasting
  
  ## First update the remaining mu and delta parameters for future time points
  for (ni in 1 : length(keep)) {
    for (tt in 2 : TT) {
      ## which indices need updating?
      wMiss = which(is.na(muPost[ni,,tt]) == TRUE)
      
      if (length(wMiss) == 1) {
        deltaPost[ni,wMiss,tt] = deltaPost[ni,wMiss,tt-1] + 
          rnorm(1, sd=sqrt(sigmaDPost[ni, wMiss]))
        
        muPost[ni,wMiss,tt] = deltaPost[ni,wMiss,tt-1] + 
          muPost[ni,wMiss,tt-1] +
          rnorm(1, sd=sqrt(sigmaUPost[ni, wMiss]))
      } else if (length(wMiss) >= 1) {
        Dd = diag(sigmaDPost[ni, wMiss])
        Du = diag(sigmaUPost[ni, wMiss])
        deltaPost[ni,wMiss,tt] = deltaPost[ni,wMiss,tt-1] + 
          chol(Dd)%*%rnorm(length(wMiss))
        
        muPost[ni,wMiss,tt] = deltaPost[ni,wMiss,tt-1] + 
          muPost[ni,wMiss,tt-1] +
          chol(Du)%*%rnorm(length(wMiss)) 
      }
    }
  }
  
  forecasts = array(NA, dim=c(length(keep), N, TT))
  
  for (ni in 1 : length(keep)) {
    for (tt in 1 : TT) {
      
      w = which(startTimes > tt)
      wMiss = which(startTimes <= tt)
      
      if (length(wMiss) == 0) {
        forecasts[ni,,tt] = outcomes[,tt]
      } else if (length(wMiss) > 0 & length(wMiss) < N) {
        mMean = muPost[ni,,tt]
        mVar = cov_sigma
        
        cMean = mMean[-w] + cov_sigma[-w,w] %*% solve(cov_sigma[w, w]) %*% 
          (outcomes[w,tt] - mMean[w])
        cVar = cov_sigma[-w,-w] -  cov_sigma[-w,w] %*% 
          solve(cov_sigma[w, w]) %*% cov_sigma[w, -w]
        
        cVar = (cVar + t(cVar)) / 2
        
        draw = cMean + chol(cVar) %*% rnorm(length(cMean))
        forecasts[ni,w,tt] = outcomes[w,tt]
        forecasts[ni,-w,tt] = draw
      } else {
        mMean = muPost[ni,,tt]
        mVar = cov_sigma
        
        draw = mMean + chol(mVar) %*% rnorm(length(mMean))
        forecasts[ni,,tt] = draw
      }
    }
  }
  
  ## Now get posterior draws of the treatment effect
  indTauPost = array(NA, dim=c(N, nTimesOut, length(keep)))
  
  for (ni in 1 : length(keep)) {
    for (ii in 1 : N) {
      tempStart = startTimes[ii]
      for (jj in 1 : nTimesOut) {
        if (tempStart + jj - 1 <= TT) {
          indTauPost[ii,jj,ni] = outcomes[ii,tempStart + jj - 1] - 
            forecasts[ni,ii,tempStart + jj -1] 
        }
      }
    }
  }
  
  l = list(indTauPost=indTauPost)
  
  return(l)
}