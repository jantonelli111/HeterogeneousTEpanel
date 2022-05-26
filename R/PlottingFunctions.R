#' Plot marginal effects of a policy with staggered adoption of treatment
#'
#' @param modelFit            The result of a call to HeterogeneousTEpanel()                          
#'          
#'
#' @return Plot with estimates and 95% credible intervals for marginal estimands 
#'         at each time point considered.
#'
#' @export


PlotMarginalEffects = function(modelFit) {
  nTimesOut = length(modelFit$Marginal$timeSpecific$estimates)
  
  plot(1 : nTimesOut, modelFit$Marginal$timeSpecific$estimates, 
       pch=18, cex=1.3, 
       ylim=range(c(modelFit$Marginal$timeSpecific$CIlower,
                    modelFit$Marginal$timeSpecific$CIupper)),
       ylab="Treatment effects",
       xlab="Time periods after treatment")
  segments(x0=1:nTimesOut, x1=1:nTimesOut,
           y0=modelFit$Marginal$timeSpecific$CIlower,
           y1=modelFit$Marginal$timeSpecific$CIupper, lwd=1.5)
  abline(h=0, lty=2, col="darkgrey")
  
}


#' Plot heterogeneous effects of a policy with staggered adoption of treatment
#'
#' @param modelFit            The result of a call to HeterogeneousTEpanel()                          
#'          
#'
#' @return Plots estimates and 95% credible intervals for the effect of each
#'         covariate considered.
#'
#' @export


PlotHeterogeneousEffects = function(modelFit) {
  Q = length(modelFit$Heterogeneous$beta$estimates)
  
  plot(1 : Q, modelFit$Heterogeneous$beta$estimates, 
       pch=18, cex=1.3, 
       ylim=range(c(modelFit$Heterogeneous$beta$CIlower,
                    modelFit$Heterogeneous$beta$CIupper)),
       ylab="",
       xlab="Covariate number")
  segments(x0=1:Q, x1=1:Q,
           y0=modelFit$Heterogeneous$beta$CIlower,
           y1=modelFit$Heterogeneous$beta$CIupper, lwd=1.5)
  abline(h=0, lty=2, col="darkgrey")
  
}
  
  
  