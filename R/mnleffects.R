#' mnlAveEffPlot2
#'
#' @param obj mixvar
#' @param varname
#' @param ColNames
#' @param R
#' @param nvals
#' @param plot
#'
#' @export
mnlAveEffPlot2 <- function (obj, varname, ColNames=NULL, R = 1500, nvals = 25, plot = TRUE, ...)
{
  data <- datatmp<- obj$datatmp
  Y    <- 1:length(obj$YFormula)
  y    <- obj$Tau.tmp
  obj  <- obj$LatentMod


  vars <- all.vars(formula(obj))[-1]
  if (any(!(vars %in% names(data)))) {
    vars <- vars[-which(!vars %in% names(data))]
  }
  rn <- vars
  var.classes <- sapply(vars, function(x) class(data[[x]]))
  TMPdata<-list()
  TMPdata$Y<-.GlobalEnv$Y
  TMPdata$datatmp<-.GlobalEnv$datatmp
  Y<<-Y
  datatmp<<-datatmp
  b <- mvrnorm(R, c(t(coef(obj))), nnet:::vcov.multinom(obj))
  d0 <- list()
  if (is.numeric(data[[varname]])) {
    s <- seq(min(data[[varname]], na.rm = T), max(data[[varname]],
                                                  na.rm = T), length = nvals)
    for (i in 1:length(s)) {
      d0[[i]] <- data
      d0[[i]][[varname]] <- s[i]
    }
  }
  if (!is.numeric(data[[varname]])) {
    s <- obj$xlevels[[varname]]
    for (j in 1:length(s)) {
      d0[[j]] <- data
      d0[[j]][[varname]] <- factor(j, levels = 1:length(s),
                                   labels = s)
    }
  }
  Xmats <- lapply(d0, function(x) model.matrix(formula(obj),
                                               data = x))
  # y <- model.response(model.frame(obj))
  if(is.matrix(y)){
    if(is.null(ColNames)) ylev <- colnames(y) else ylev <- ColNames
  } else ylev <- levels(y)

  b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
  xb <- lapply(Xmats, function(x) lapply(1:nrow(b), function(z) cbind(1,
                                                                      exp(x %*% t(matrix(c(t(b[z, ])), ncol = ncol(coef(obj)),
                                                                                         byrow = T))))))
  probs <- lapply(xb, function(x) lapply(x, function(z) z/rowSums(z)))
  out.ci2 <- lapply(probs, function(x) sapply(x, colMeans))
  out.ci <- lapply(out.ci2, apply, 1, quantile, c(0.5, 0.025,0.975,0.16,0.84))
  out.ciddendum <- lapply(out.ci2, apply, 1, function(x)c(min(x),max(x)))
  tmp <- data.frame(mean = do.call("c", lapply(out.ci, function(x) x[1, ])),
                    lower1 = do.call("c", lapply(out.ci, function(x) x[2,])),
                    upper1 = do.call("c", lapply(out.ci, function(x) x[3,])),
                    lower2 = do.call("c", lapply(out.ci, function(x) x[4,])),
                    upper2 = do.call("c", lapply(out.ci, function(x) x[5,])),
                    y = rep(ylev, length(out.ci)))
  if (is.numeric(data[[varname]])) {
    tmp$s <- rep(s, each = length(ylev))
  }
  else {
    tmp$s <- factor(rep(1:length(s), each = length(ylev)),
                    labels = s)
  }
  tmp$min <- do.call("c", lapply(out.ciddendum, function(x) x[1, ]))
  tmp$max <- do.call("c", lapply(out.ciddendum, function(x) x[2, ]))
  SummaryOutput <<- tmp
  tmp$min <- NULL
  tmp$max <- NULL
  cat("Variable >> SummaryOutput << created\n")
  # if (plot) {
  #   if (is.factor(data[[varname]])) {
  #     pl <- xyplot(mean ~ s | y, data = tmp, xlab = "",
  #                  ylab = "", ..., lower = tmp$lower,
  #                  upper = tmp$upper, prepanel = prepanel.ci,
  #                  scales = list(x = list(at = 1:length(s), labels = s)),
  #                  panel = function(x, y, lower, upper, subscripts) {
  #                       panel.points(x, y, col = "black", lty = 1,pch = 16)
  #                       panel.segments(x, lower[subscripts], x, upper[subscripts],lty = 1, col = "blue")
  #                  }
  #     )
  #   }
  #   if (is.numeric(data[[varname]])) {
  #     pl <- xyplot(mean ~ s | y, data = tmp, xlab = "",
  #                  ylab = "Predicted Value", ..., lower = tmp$lower,
  #                  upper = tmp$upper, prepanel = prepanel.ci, panel = panel.ci,
  #                  zl = FALSE)
  #   }
  #   return(pl)
  # }
  # else {
  #   return(tmp)
  # }
  mmcolors <- c("black","darkgrey")
  mmlty    <- c(1,1)
  mmpch    <- c(1,4)
  relation="same"
  PLOT <- xyplot( mean ~ s| y , data=tmp, xlab="", ylab="",
                  upper1 = tmp$upper1, lower1 = tmp$lower1,upper2 = tmp$upper2, lower2 = tmp$lower2,
                  prepanel = mixtureVAR:::prepanel.default.xyplot2,
                  panel = function(x, y,upper1, lower1,upper2, lower2, subscripts, ...){
                    upper1 <- upper1[subscripts]
                    lower1 <- lower1[subscripts]
                    upper2 <- upper2[subscripts]
                    lower2 <- lower2[subscripts]
                    panel.polygon(c(x, rev(x)), c(upper1, rev(lower1)), fill=.5, col="lightgrey", border=FALSE,...)
                    panel.polygon(c(x, rev(x)), c(upper2, rev(lower2)), fill=.5, col="grey", border=FALSE,...)
                    panel.xyplot(x, y, type='b', cex=0.2, lty=1, col=mmcolors[1], pch = mmpch[1],...)
                  },
                  strip=strip.custom(  ),
                  strip.left=FALSE)
  PLOT <- update(PLOT, scales = list(y=list(relation=relation,alternating=3,at=c(0,.2,.4,.6,.8,1), at.lab=c(0,.2,.4,.6,.8,1)),
                                     x=list(alternating=1, tck=c(1,0))),ylim=c(0,1),
                 par.settings=list(strip.background = list( col="transparent" ),
                                   strip.border     = list( col = NA )))
  #useOuterStrips(PLOT,strip.left=FALSE)
  .GlobalEnv$Y <- TMPdata$Y
  .GlobalEnv$datatmp <- TMPdata$datatmp

  PLOT
}
