#' #' Make VAR Object from msvar
#' #'
#' #' @param msvar
#' #' @param y
#' #' @param nstep
#' #'
#' #' @return list
#' #'
#' #' @export
#' makeVARobjfrommsvar <- function(msvar,y,nstep){
#'   outputlist <- list()
#'   for(h in 1:msvar$h){
#'     m     <- msvar$m
#'     p     <- msvar$p
#'     ncoef <- (m * p) + 1
#'     ndum  <- m + 1
#'     n     <- nrow(msvar$fp) + p
#'     capT  <- n - p + ndum
#'     Ts    <- n - p
#'     Bh    <- msvar$hreg$Bk[,,h]
#'     u     <- msvar$hreg$e[,,h]
#'     Sh    <- crossprod(u)/Ts # Ist das korrekt oder doch lieber msvar$hreg$Sigmak[,,h]
#'     Sh1   <- crossprod(u)/capT
#'     intercept     <- Bh[(m * p + 1), ]
#'     ar.coefs      <- t(Bh[1:(m * p), ])
#'     dim(ar.coefs) <- c(m, m, p)
#'     ar.coefs      <- aperm(ar.coefs, c(2, 1, 3))
#'     exog.coefs    <- NA
#'     num.exog      <- 0
#'     X       <- msvar$init.model$X[(p+2):(n+1),]
#'     XX      <- crossprod(X)
#'     Y       <- msvar$init.model$Y[(p+2):(n+1),]
#'     pfit    <- list(capT = capT, ncoef = ncoef, num.exog = num.exog, Sh1 = Sh1)
#'     output  <- list(intercept = intercept, ar.coefs = ar.coefs,
#'                     Bhat = Bh, vcv = Sh, exog.coefs = exog.coefs, residuals = u,
#'                     mean.S = Sh, hstar = XX, X = X, Y = Y, y = y, pfit = pfit)
#'     class(output)           <- c("VAR")
#'     attr(output, "eqnames") <- colnames(y)
#'     outputlist[[h]]         <- MSBVAR::irf(output,nstep)
#'   }
#'   outputlist
#' }
#'
#' #' Make VAR Object from LSTAR
#' #'
#' #' @param msvar
#' #' @param y
#' #' @param nstep
#' #'
#' #' @return list
#' #'
#' #' @export
#' makeVARobjfromLSTAR <- function(msvar,y,nstep){
#'   outputlist <- list()
#'   for(h in 1:2){
#'     m     <- msvar$k
#'     p     <- msvar$lag
#'     ncoef <- (m * p) + 1
#'     ndum  <- m + 1
#'     n     <- msvar$T
#'     capT  <- n - p + ndum
#'     Ts    <- n - p
#'     Bh    <- t(msvar$coefficients[,((h-1)*ncoef + 1:ncoef)])
#'     u     <- msvar$residuals[round(msvar$model.specific$regime[-c(1:p)])==h,]
#'     Sh    <- crossprod(u)/Ts
#'     Sh1   <- crossprod(u)/capT
#'     intercept     <- Bh[1, ]
#'     ar.coefs      <- t(Bh[-1, ])
#'     dim(ar.coefs) <- c(m, m, p)
#'     ar.coefs      <- aperm(ar.coefs, c(2, 1, 3))
#'     exog.coefs    <- NA
#'     num.exog      <- 0
#'     X       <- msvar$model[-(1:p),m + ((h-1)*ncoef + 1:ncoef)]
#'     XX      <- crossprod(X)
#'     Y       <- msvar$model[-(1:p),1:m]
#'     pfit    <- list(capT = capT, ncoef = ncoef, num.exog = num.exog, Sh1 = Sh1)
#'     output  <- list(intercept = intercept, ar.coefs = ar.coefs,
#'                     Bhat = Bh, vcv = Sh, exog.coefs = exog.coefs, residuals = u,
#'                     mean.S = Sh, hstar = XX, X = X, Y = Y, y = msvar$model[,1:m], pfit = pfit)
#'     class(output)           <- c("VAR")
#'     attr(output, "eqnames") <- colnames(msvar$model[,1:m])
#'     outputlist[[h]]         <- MSBVAR::irf(output,nstep)
#'   }
#'   outputlist
#' }
#'
#' #' Make VAR Object from STAR (Version 2)
#' #'
#' #' @param msvar
#' #' @param y
#' #' @param nstep
#' #'
#' #' @return list
#' #'
#' #' @export
#' makeVARobjfromLSTARVersion2 <- function(msvar,y,nstep){
#'   outputlist <- list()
#'   for(h in 1:2){
#'     m     <- msvar$k
#'     p     <- msvar$lag
#'     ncoef <- (m * p) + 1
#'     ndum  <- m + 1
#'     n     <- msvar$T
#'     capT  <- n - p + ndum
#'     Ts    <- n - p
#'     Bh    <- t(msvar$coefficients[,((h-1)*ncoef + 1:ncoef)])
#'     X     <- msvar$model[-(1:p),m + ((h-1)*ncoef + 1:ncoef)]
#'     XX    <- crossprod(X)
#'     Y     <- msvar$model[-(1:p),1:m]
#'     u     <- Y-X%*%Bh
#'     Sh    <- crossprod(u)/Ts
#'     Sh1   <- crossprod(u)/capT
#'     intercept     <- Bh[1, ]
#'     ar.coefs      <- t(Bh[-1, ])
#'     dim(ar.coefs) <- c(m, m, p)
#'     ar.coefs      <- aperm(ar.coefs, c(2, 1, 3))
#'     exog.coefs    <- NA
#'     num.exog      <- 0
#'     pfit    <- list(capT = capT, ncoef = ncoef, num.exog = num.exog, Sh1 = Sh1)
#'     output  <- list(intercept = intercept, ar.coefs = ar.coefs,
#'                     Bhat = Bh, vcv = Sh, exog.coefs = exog.coefs, residuals = u,
#'                     mean.S = Sh, hstar = XX, X = X, Y = Y, y = msvar$model[,1:m], pfit = pfit)
#'     class(output)           <- c("VAR")
#'     attr(output, "eqnames") <- colnames(msvar$model[,1:m])
#'     outputlist[[h]]         <- MSBVAR::irf(output,nstep)
#'   }
#'   outputlist
#' }
#'
#' #' resid LVSTAR (Version 1)
#' #'
#' #' @param msvar
#' #' @param y
#' #' @param nstep
#' #'
#' #' @return list
#' #'
#' #' @export
#' residLVSTARVersion1 <- function(msvar,y,nstep){
#'   outputlist <- list()
#'   for(h in 1:2){
#'     p     <- msvar$lag
#'     u     <- msvar$residuals[round(msvar$model.specific$regime[-c(1:p)])==h,]
#'     outputlist[[h]] <- u
#'   }
#'   outputlist
#' }
#'
#' #' resid LVSTAR (Version 2)
#' #'
#' #' @param msvar
#' #'
#' #' @return list
#' #'
#' #' @export
#' residLVSTARVersion2 <-  function(msvar){
#'   outputlist <- list()
#'   for(h in 1:2){
#'     m     <- msvar$k
#'     p     <- msvar$lag
#'     ncoef <- (m * p) + 1
#'     ndum  <- m + 1
#'     n     <- msvar$T
#'     capT  <- n - p + ndum
#'     Ts    <- n - p
#'     Bh    <- t(msvar$coefficients[,((h-1)*ncoef + 1:ncoef)])
#'     X     <- msvar$model[-(1:p),m + ((h-1)*ncoef + 1:ncoef)]
#'     XX    <- crossprod(X)
#'     Y     <- msvar$model[-(1:p),1:m]
#'     u     <- Y-X%*%Bh
#'     outputlist[[h]] <- u
#'   }
#'   outputlist
#' }
