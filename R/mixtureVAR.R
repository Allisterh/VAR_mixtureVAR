#' Estimate mixtureVARs
#'
#' @param Y timeseries
#' @param P list with lag length per model for the mixture
#' @param EXO list of length length(P) of exogenous timeseries
#' @param data data frame containg the additional variables
#' @param maxiter maximum number of iteration for the iterative fitting procedure
#' @param debug if TRUE some debug output is provided
#' @param eps positive convergence tolerance
#' @param nreps number of inital models to be estimated, for lowering the problem of possible convergence to local maxima
#' @param Theta_start optional starting parameters
#' @param multinommod Model formula for latent model.
#'
#' @return list
#'
#' @export
mixtureVAR <- function(YFormula,
                   LatentFormula = NULL,
                   LatentLag = NULL,
                   EXOFormula = NULL,
                   data,
                   P,
                   maxiter = 1000,
                   debug = FALSE,
                   eps = 1e-8,
                   nreps = 5,
                   Theta_start = NULL) {
  ################
  #### Preparations and stating values ####
  ################
  p     <- max(unlist(P))
  PDiff <- lapply(P, function(pp) p - pp)
  if(plyr::is.formula(YFormula[[1]])){
    Y <- lapply(YFormula,function(form)t(model.matrix(form,data=data)))
  } else Y <- YFormula
  X <- mapply(function(p,y) t(cbind(1, embed(t(y), p + 1)[,-seq_len(nrow(y))])),P,Y,SIMPLIFY = FALSE)
  if (!is.null(EXOFormula[[1]])) {
    if(plyr::is.formula(EXOFormula)){
      EXO <- lapply(EXOFormula,function(form) t(model.frame(form,data=data)))
    } else EXO <- EXOFormula
    X   <- mapply(function(x, exo) {
      xout <- rbind(x, exo[, (ncol(exo) - ncol(x) + 1):ncol(exo)])
      xout
    }, X, EXO, SIMPLIFY = FALSE)
    EXO.orig <- EXO
    EXO     <- mapply(function(exo, p) exo[, -c(1:(p))], EXO.orig, P, SIMPLIFY = FALSE)
  }
  K     <- length(P) # number of Models
  returnlist <- list()
  Y.orig<-Y
  Y     <- mapply(function(p,y) y[, -c(1:(p))],P,Y,SIMPLIFY = FALSE)

  rep <- 1
  if (nreps > 1)
    pb <- txtProgressBar(min = 0, max = nreps, style = 3)
  while (rep <= nreps) {
    diditwork <- try({
      Alpha <- lapply(rdirichlet(1, rep(1 / K, K)), identity)
      Sigma <- lapply(Y, function(y) diag(nrow(y)))
      if (nreps > 1 & !is.null(Theta_start)) {
        stop("using Theta_start for more than one nreps is not allowed yet")
      } else if (nreps == 1 & !is.null(Theta_start)) {
        Theta <- Theta_start
        Theta.names <- lapply(Theta, dimnames)
      } else {
        if (!is.null(EXOFormula)) {
          Theta <- mapply(VARlist2, Y, P, EXO, SIMPLIFY = FALSE)
        } else {
          Theta <- mapply(VARlist2, Y, P, SIMPLIFY = FALSE)
        }
        Theta.names <- lapply(Theta, dimnames)
        Theta <- mapply("+", Theta,
                        lapply(Theta, function(bla)
                          matrix(
                            rnorm(length(as.vector(bla)), sd = 0.001),
                            nrow = nrow(bla),
                            ncol = ncol(bla)
                          )), SIMPLIFY = FALSE)
      }
      E     <- mapply(res, Y, X, Theta, PDiff, SIMPLIFY = FALSE)

      liklhood    <- mapply(likl, Alpha, Sigma, E, PDiff)
      Tau         <- do.call(c, apply(apply(liklhood, 1, function(x) x / sum(x)), 1, list))
      Tau_alt     <- lapply(Tau, function(xx) 5 * xx) # ist die Summe==1 bei den Taus noch gegeben?? Ist das relevant?
      Theta_alt   <- lapply(Theta, function(xx) 5 * xx)
      Sigma_alt   <- lapply(Sigma, function(xx) 5 * xx)

      if (debug == TRUE) {
        Sigma_hist <- list()
        Theta_hist <- list()
      }

      i <- 0
      while (max(c(
        unlist(Theta_alt) - unlist(Theta),
        unlist(Sigma_alt) - unlist(Sigma)
      )) > eps & i <= maxiter) {
        i <- i + 1
        print(i)
        ################
        #### E-Step ####
        ################

        liklhood    <- mapply(likl, Alpha, Sigma, E, PDiff)
        if(!is.null(LatentFormula) & !is.null(LatentLag)) {
          Tau.tmp    <- apply(liklhood, 1, function(x) x / sum(x))
          datatmp    <- cbind(t(Tau.tmp),data[-c(1:p),])
          colnames(datatmp)[1:length(Y)] <- paste0("Y",1:length(Y))
          laggedtaus <- sapply(LatentLag, function(llag){
            tmplag <- matrix(NA,nrow=nrow(datatmp),ncol=nrow(Tau.tmp)-1,
                             dimnames = list(colnames(Tau.tmp),
                                             paste0("tau",1:(nrow(Tau.tmp)-1),"_lag",llag)))
            tmplag[(llag+1):nrow(datatmp),] <- Tau.tmp[-nrow(Tau.tmp),1:(ncol(Tau.tmp)-llag)]
            tmplag
          })
          laggedtausnames <- unlist(lapply(LatentLag, function(llag){
            paste0("tau",1:(nrow(Tau.tmp)-1),"_lag",llag)
          }))
          colnames(laggedtaus)<-laggedtausnames
          datatmp <- data.frame(datatmp,laggedtaus)
          maxLatentLag <- max(which(apply(datatmp,1,function(tmpna)any(is.na(tmpna)))))
          datatmp <- datatmp[-c(1:maxLatentLag),]

          LatentFormulaLag <- update(LatentFormula,paste(c(" ~ .",colnames(laggedtaus)),sep="",collapse="+"))
          LatentMod <-
            multinom(update(LatentFormulaLag,reformulate(".",parse(text=paste0("cbind(",paste0("Y",1:length(Y),collapse=","),")"))[[1]])),
                     data=datatmp, MaxNWts = 100000, trace = FALSE)
          Tau <- predict(LatentMod, type = "probs")
          Tau <- rbind(t(Tau.tmp[,1:maxLatentLag]),Tau)
          Tau <- do.call(c, apply(Tau, 2, list))
        } else if(!is.null(LatentFormula)) {
          Tau.tmp    <- apply(liklhood, 1, function(x) x / sum(x))
          datatmp <- cbind(t(Tau.tmp),data[-c(1:p),])
          colnames(datatmp)[1:length(Y)] <- paste0("Y",1:length(Y))
          LatentMod <-
            multinom(update(LatentFormula,reformulate(".",parse(text=paste0("cbind(",paste0("Y",1:length(Y),collapse=","),")"))[[1]])), data=datatmp,
                     MaxNWts = 100000,
                     trace = FALSE)
          Tau <- predict(LatentMod, type = "probs")
          Tau <- do.call(c, apply(Tau, 2, list))
        } else  {
          Tau <- do.call(c, apply(apply(liklhood, 1, function(x)x / sum(x)), 1, list))
        }
        ################
        #### M-Step ####
        ################

        Tau_alt     <- Tau
        Theta_alt   <- Theta
        Sigma_alt   <- Sigma

        Alpha       <- lapply(Tau, mean)
        Theta       <- mapply(thetastep2, Tau, X, Y, PDiff, SIMPLIFY = FALSE)
        E           <- mapply(res, Y, X, Theta, PDiff, SIMPLIFY = FALSE)
        Sigma       <- mapply(function(tau, e) 1 / sum(tau) * e %*% (tau * t(e)), Tau, E, SIMPLIFY = FALSE)
        if (debug == TRUE) {
          if (i %% 10 == 0)
            cat(sum(apply( sapply(Tau, "[") * liklhood, 2, sum )), "\n")
          Sigma_hist[[i]] <- Sigma
          Theta_hist[[i]] <- Theta
        }
      }
      if (!i < maxiter)
        warning("Iteration limit reached. Try to increase 'maxiter'.")

      Theta <- mapply(function(theta, theta.names) {
        dimnames(theta) <- theta.names
        theta
      }, Theta, Theta.names, SIMPLIFY = FALSE)

      ### Hier finales submodel berechnen
      if(!is.null(LatentFormula) & !is.null(LatentLag)) {
        Tau.tmp    <- apply(liklhood, 1, function(x) x / sum(x))
        datatmp    <- cbind(t(Tau.tmp),data[-c(1:p),])
        colnames(datatmp)[1:length(Y)] <- paste0("Y",1:length(Y))
        laggedtaus <- sapply(LatentLag, function(llag){
          tmplag <- matrix(NA,nrow=nrow(datatmp),ncol=nrow(Tau.tmp)-1,
                           dimnames = list(colnames(Tau.tmp),
                                           paste0("tau",1:(nrow(Tau.tmp)-1),"_lag",llag)))
          tmplag[(llag+1):nrow(datatmp),] <- Tau.tmp[-nrow(Tau.tmp),1:(ncol(Tau.tmp)-llag)]
          tmplag
        })
        laggedtausnames <- unlist(lapply(LatentLag, function(llag){
          paste0("tau",1:(nrow(Tau.tmp)-1),"_lag",llag)
        }))
        colnames(laggedtaus)<-laggedtausnames
        datatmp <- data.frame(datatmp,laggedtaus)
        maxLatentLag <- max(which(apply(datatmp,1,function(tmpna)any(is.na(tmpna)))))
        datatmp <- datatmp[-c(1:maxLatentLag),]

        LatentFormulaLag <- update(LatentFormula,paste(c(" ~ .",colnames(laggedtaus)),sep="",collapse="+"))
        LatentMod <-
          multinom(update(LatentFormulaLag,reformulate(".",parse(text=paste0("cbind(",paste0("Y",1:length(Y),collapse=","),")"))[[1]])),
                   data=datatmp, MaxNWts = 100000, trace = FALSE)
        Tau <- predict(LatentMod, type = "probs")
        Tau <- rbind(t(Tau.tmp[,1:maxLatentLag]),Tau)
        Tau <- do.call(c, apply(Tau, 2, list))
      } else if(!is.null(LatentFormula)) {
        Tau.tmp    <- apply(liklhood, 1, function(x) x / sum(x))
        datatmp <- cbind(t(Tau.tmp),data[-c(1:p),])
        colnames(datatmp)[1:length(Y)] <- paste0("Y",1:length(Y))
        LatentMod <-
          multinom(update(LatentFormula,reformulate(".",parse(text=paste0("cbind(",paste0("Y",1:length(Y),collapse=","),")"))[[1]])), data=datatmp,
                   MaxNWts = 100000,
                   trace = FALSE)
        Tau <- predict(LatentMod, type = "probs")
        Tau <- do.call(c, apply(Tau, 2, list))
      } else  {
        Tau <- do.call(c, apply(apply(liklhood, 1, function(x)x / sum(x)), 1, list))
      }

      if (debug == TRUE)
        list(
          Theta = Theta,
          Sigma = Sigma,
          Tau = Tau,
          Tau.tmp = t(Tau.tmp),
          Alpha = Alpha,
          e = E,
          iter = i,
          Sigma_hist = Sigma_hist,
          Theta_hist = Theta_hist,
          P = unlist(P),
          data=data,
          datatmp = datatmp,
          Y = Y,
          X = X,
          EXO.orig = if(!is.null(EXOFormula)) EXO.orig else NULL,
          EXOFormula = if (!is.null(EXOFormula)) EXOFormula else NULL,
          Y.orig = Y.orig,
          YFormula = YFormula,
          LatentFormula = if (!is.null(LatentFormula)) LatentFormula else NULL,
          LatentLag = if (!is.null(LatentLag)) LatentLag else NULL,
          LatentFormulaLag = if (!is.null(LatentLag)) LatentFormulaLag else NULL,
          LatentMod = if (!is.null(LatentFormula)) LatentMod else NULL,
          maxiter = maxiter
        )
      else
        list(
          Theta = Theta,
          Sigma = Sigma,
          Tau = Tau,
          Tau.tmp = t(Tau.tmp),
          Alpha = Alpha,
          e = E,
          iter = i,
          data=data,
          datatmp = datatmp,
          P = unlist(P),
          Y = Y,
          X = X,
          EXO.orig = if(!is.null(EXOFormula)) EXO.orig else NULL,
          EXOFormula = if (!is.null(EXOFormula)) EXOFormula else NULL,
          Y.orig = Y.orig,
          YFormula = YFormula,
          LatentFormula = if (!is.null(LatentFormula)) LatentFormula else NULL,
          LatentLag = if (!is.null(LatentLag)) LatentLag else NULL,
          LatentFormulaLag = if (!is.null(LatentLag)) LatentFormulaLag else NULL,
          LatentMod = if (!is.null(LatentFormula)) LatentMod else NULL,
          maxiter = maxiter
        )
    }, silent = TRUE)
    if (!inherits(diditwork, "try-error") | nreps == 1) {
      returnlist[[rep]] <- diditwork
      if (nreps > 1)
        setTxtProgressBar(pb, rep)
      rep <- rep + 1
    } else {
      #cat("error\n")
      rep <- rep
    }
  }
  if (nreps == 1) {
    ret <- returnlist[[1]]
    if (!inherits(ret, "try-error"))
      class(ret) <- "mixvar"
    return(ret)
  }
  else {
    liks <-
      sapply(returnlist, function(mod)
        sum(apply(
          sapply(mod$Tau, "[") * mapply(likl, mod$Alpha, mod$Sigma, mod$e, PDiff) ,
          2,
          sum
        )))
    returnlist <-
      list(
        finalmodel = returnlist[[which.max(liks)]],
        allmodels = returnlist,
        likelihoods = liks
      )
    class(returnlist) <- "mixvar"
    class(returnlist$finalmodel) <- "mixvar"
    return(returnlist)
  }

}
