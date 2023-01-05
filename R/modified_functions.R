summary.ergm.david <- function (object, ..., 
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE)
{
  # Warn if the object was produced by an earlier version of ergm.
  myver <- packageVersion("ergm")
  objver <- NVL(object$ergm_version, as.package_version("3.9.4")) # 3.9.4 was the last version that didn't store the version information.
  nextver <- as.package_version(paste(objver$major, objver$minor+1, sep="."))
  if(objver < paste(myver$major, myver$minor, sep=".")){
    warning(paste0("This object was fit with ", sQuote("ergm"), " version ", objver, " or earlier. Summarizing it with version ", nextver, " or later may return incorrect results or fail."))
  }
  
  if("digits" %in% names(list(...))) warn("summary.ergm() no longer takes a digits = argument.")
  control <- object$control
  pseudolikelihood <- object$estimate=="MPLE"
  independence <- NVL(object$MPLE_is_MLE, is.dyad.independent(object))
  coef <- coef(object)
  
  ans <- list(formula=object$formula,
              call=object$call,
              correlation=correlation,
              offset = object$offset,
              drop = NVL(object$drop, rep(0,length(object$offset))),
              estimable = NVL(object$estimable, rep(TRUE,length(object$offset))),
              covariance=covariance,
              pseudolikelihood=pseudolikelihood,
              independence=independence,
              estimate=object$estimate,
              estimate.desc=object$estimate.desc,
              control=object$control)
  
  asycov <- vcov(object, sources=if(total.variation) "all" else "model")
  asyse <- sqrt(diag(asycov))
  # Convert to % error  
  est.se <- sqrt(diag(vcov(object, sources="estimation")))
  mod.se <- sqrt(diag(vcov(object, sources="model")))
  tot.se <- sqrt(diag(vcov(object, sources="all")))
  est.pct <- rep(NA,length(est.se))
  if(any(!is.na(est.se))){
    # We want (sqrt(V.model + V.MCMC)-sqrt(V.model))/sqrt(V.model + V.MCMC) * 100%,
    est.pct[!is.na(est.se)] <- ifelse(est.se[!is.na(est.se)]>0, round(100*(tot.se[!is.na(est.se)]-mod.se[!is.na(est.se)])/tot.se[!is.na(est.se)]), 0)
  }
  
  zval <- coef / asyse
  pval <- 2 * pnorm(q=abs(zval), lower.tail=FALSE)
  or <- exp(coef)
  lci <- exp(coef - 1.96 * asyse)
  uci <- exp(coef + 1.96 * asyse)
  
  count <- 1
  coefmat <- cbind(
    `Estimate` = coef,
    `Std. Error` = asyse,
    `MCMC %` = est.pct,
    `Wald` = zval,
    `OR` = round(or, 3),
    `Lower` = round(lci, 3),
    `Upper` = round(uci, 3),
    `Pr(>|z|)` = pval)
  
  rownames(coefmat) <- param_names(object)
  
  devtext <- "Deviance:"
  if (object$estimate!="MPLE" || !independence || object$reference != as.formula(~Bernoulli)) {
    if (pseudolikelihood) {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect.\n"
    } 
    else if(object$estimate == "MLE" && any(is.na(est.se) & !ans$offset & !ans$drop==0 & !ans$estimable) && 
            (!independence || control$force.main) ) {
      ans$message <- "\nWarning:  The standard errors are suspect due to possible poor convergence.\n"
    }
  } else {
    ans$message <- "\nFor this model, the pseudolikelihood is the same as the likelihood.\n"
  }
  mle.lik<-try(logLik(object,...), silent=TRUE)
  
  if(inherits(mle.lik,"try-error")) ans$objname<-deparse(substitute(object))
  else if(!is.na(mle.lik)){
    # Only evaluate the null likelihood if the MLE likelihood is defined.
    null.lik<-try(logLikNull(object,...), silent=TRUE)
    ans$null.lik.0 <- is.na(null.lik)
    
    df <- length(coef)
    dyads<- sum(as.rlebdm(object$constrained, object$constrained.obs, which="informative"))
    rdf <- dyads - df
    ans$devtable <- matrix(c(if(is.na(null.lik)) 0 else -2*null.lik, -2*mle.lik,
                             c(dyads, rdf)), 2,2, dimnames=list(c("Null","Residual"),
                                                                c("Resid. Dev", "Resid. Df")))
    ans$devtext <- devtext
    
    ans$aic <- AIC(object)
    ans$bic <- BIC(object)
    ans$mle.lik <- ERRVL(mle.lik, NA)
    ans$null.lik <- ERRVL(null.lik, NA)
  }else ans$devtable <- NA
  
  ans$coefficients <- coefmat
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}