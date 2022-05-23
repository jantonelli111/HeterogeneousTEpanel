#' Estimate causal effects of a policy with staggered adoption of treatment
#'
#' @param outcomes            An N by T matrix of outcomes for each unit over T time periods
#' @param treatments          An N by T matrix of binary indicators indicating treatment
#'                            status for each unit at each point in time
#' @param covariates          An N by Q data frame of covariate values for each unit
#' @param nScans              The number of MCMC scans to run
#' @param nBurn               The number of MCMC scans that will be dropped as a burn-in
#' @param thin                This number represents how many iterations between each scan
#'                            that is kept
#' @param smoothEffects       Whether or not to smooth treatment effects over time
#' @param newX                An M by Q matrix of covariate values to estimate the 
#'                            treatment effect at. The default is NULL.
#' @param zeroMat             An optional matrix of indices of entries of the inverse
#'                            covariance matrix of the time series to be set to zero. 
#'                            If NULL, then no elements will be forced to zero. You do
#'                            not need to specify both the (i,j) and (j,i) elements
#'                            as the solution is necessarily symmetric. For more details
#'                            see the helpfile for glasso().
#' @param priorA              First hyperprior for inverse gamma prior distribution
#'                            on variance components for the time series state parameters.
#'                            We don't recommend changing this parameter, but generally
#'                            larger values lead to smoother time series forecasts.
#' @param priorB              Second hyperprior for inverse gamma prior distribution
#'                            on variance components for the time series state parameters.
#'                            We don't recommend changing this parameter, but generally
#'                            smaller values lead to smoother time series forecasts.                           
#'          
#'
#' @importFrom glasso glasso
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom splines ns
#' @return Description
#'
#' @export


HeterogeneousTEpanel = function(outcomes,
                                treatments,
                                covariates=NULL,
                                nTimesOut = 10,
                                nScans = 5000, 
                                nBurn=2000, 
                                thin=6,
                                zeroMat = NULL,
                                smoothEffects = TRUE,
                                newX=NULL,
                                priorA = 0.33,
                                priorB = 0.0000033) {
  
  ## Standardize data and then un-standardize effect estimates after
  sdSave = rep(NA, nrow(outcomes))
  for (ii in 1 : nrow(outcomes)) {
    sdSave[ii] = sd(as.numeric(outcomes[ii,]), na.rm=TRUE)
    outcomes[ii,] = (outcomes[ii,] - mean(as.numeric(outcomes[ii,]), na.rm=TRUE))/
      sd(as.numeric(outcomes[ii,]), na.rm=TRUE)
  }
  
  ## Estimate the residual covariance matrix
  covEst = CovEstStructural(covariates = covariates, treatments=treatments,
                            outcomes=outcomes, zeroMat=zeroMat)
  
  ## Get posterior distribution of individual treatment effects
  JointAutoSplines = JointEstStructural(covariates=covariates, 
                                        treatments=treatments,
                                        outcomes=outcomes,
                                        cov_sigma=covEst,
                                        nTimesOut = nTimesOut,
                                        nScans = nScans, 
                                        nBurn=nBurn, 
                                        thin=thin,
                                        priorA = priorA,
                                        priorB = priorB)
  
  indEffectsPosterior = JointAutoSplines$indTauPost
  
  ## Un-standardize results
  for (ii in 1 : nrow(outcomes)) {
    indEffectsPosterior[ii,,] = indEffectsPosterior[ii,,] * sdSave[ii]
  }
  
  ## Obtain marginal effects
  if (smoothEffects == TRUE) {
    ATEresults = DeltaQsmooth(indTauPost = indEffectsPosterior,
                              df = 3)
  } else {
    ATEresults = DeltaQ(indTauPost = indEffectsPosterior, 
                        nTimesOut = nTimesOut)
  }
  
  ## Obtain heterogeneous treatment effects
  if (is.null(covariates)) {
    HEresults = NULL
  } else {
    if (smoothEffects == TRUE) {
      HEresults = Linear2StageSmooth(indTauPost=indEffectsPosterior,
                                     covariates=covariates,
                                     newX=newX,
                                     df=3)
    } else {
      HEresults = Linear2Stage(indTauPost=indEffectsPosterior,
                               covariates=covariates,
                               newX=newX)
    } 
  }
  
  
  l = list(Marginal = ATEresults,
           Heterogeneous = HEresults,
           indEffectsPosterior = indEffectsPosterior) 
  
  return(l)
}
