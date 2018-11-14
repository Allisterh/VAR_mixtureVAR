#' data plot
#'
#' @param x
#' @param iname
#'
#' @return list
#'
#' @export
dataplot <- function(x, iname) {
  impulses <- x$irf[[iname]]
  range <- range(impulses)
  upper <- NULL
  lower <- NULL
  x$boot <- FALSE
  if (x$boot) {
    upper <- x$Upper[[iname]]
    lower <- x$Lower[[iname]]
    range <- range(cbind(impulses, upper, lower))
  }
  text1 <- paste("Orthogonal Impulse Response from",
                 iname, sep = " ")
  text2 <- ""
  result <- list(impulses = impulses, upper = upper, lower = lower,
                 range = range, text1 = text1, text2 = text2)
  return(result)
}

plot.multiple <- function(dp, nc = nc, main=NULL, ylim = NULL, mar.multi=NULL, oma.multi = NULL,
                          ylab = NULL, lty = NULL, lwd=NULL, col=NULL, adj.mtext=NULL, padj.mtext = NULL,
                          col.mtext = NULL, sub=NULL, ...) {
  x <- dp$impulses
  y <- dp$upper
  z <- dp$lower
  ifelse(is.null(main), main <- dp$text1, main <- main)
  ifelse(is.null(sub), sub <- dp$text2, sub <- sub)
  ifelse(is.null(ylim), ylim <- dp$range, ylim <- ylim)
  range <- range(c(x, y, z))
  nvr <- ncol(x)
  if (missing(nc)) {
    nc <- ifelse(nvr > 4, 2, 1)
  }
  nr <- ceiling(nvr/nc)
  par(mfrow = c(nr, nc), mar = mar.multi, oma = oma.multi)
  if (nr > 1) {
    for (i in 1:(nvr - nc)) {
      ifelse(is.null(ylab), ylabel <- colnames(x)[i],
             ylabel <- ylab)
      xy <- xy.coords(x[, i])
      plot(xy, axes = FALSE, type = "l", ylab = ylabel,
           ylim = ylim, col = col[1], lty = lty[1], lwd = lwd[1],
           ...)
      axis(2, at = pretty(range)[-1])
      abline(h = 0, col = "red")
      if (!is.null(y))
        lines(y[, i], col = col[3], lty = lty[3], lwd = lwd[3])
      if (!is.null(z))
        lines(z[, i], col = col[3], lty = lty[3], lwd = lwd[3])
      box()
    }
    for (j in (nvr - nc + 1):nvr) {
      ifelse(is.null(ylab), ylabel <- colnames(x)[j],
             ylabel <- ylab)
      xy <- xy.coords(x[, j])
      plot(xy, axes = FALSE, type = "l", ylab = ylabel,
           ylim = ylim, col = col[1], lty = lty[1], lwd = lwd[1],
           ...)
      axis(2, at = pretty(range)[-1])
      axis(1, at = 1:(nrow(x)), labels = c(0:(nrow(x) -
                                                1)))
      box()
      abline(h = 0, col = "red")
      if (!is.null(y))
        lines(y[, j], col = col[3], lty = lty[3], lwd = lwd[3])
      if (!is.null(z))
        lines(z[, j], col = col[3], lty = lty[3], lwd = lwd[3])
    }
    mtext(main, 3, line = 2, outer = TRUE, adj = adj.mtext,
          padj = padj.mtext, col = col.mtext, ...)
    mtext(sub, 1, line = 4, outer = TRUE, adj = adj.mtext,
          padj = padj.mtext, col = col.mtext, ...)
  }
  else {
    for (j in 1:nvr) {
      ifelse(is.null(ylab), ylabel <- colnames(x)[j],
             ylabel <- ylab)
      xy <- xy.coords(x[, j])
      plot(xy, type = "l", ylab = ylabel, ylim = ylim,
           col = col[1], lty = lty[1], lwd = lwd[1], ...)
      if (!is.null(y))
        lines(y[, j], col = col[3], lty = lty[3], lwd = lwd[3])
      if (!is.null(z))
        lines(z[, j], col = col[3], lty = lty[3], lwd = lwd[3])
      abline(h = 0, col = "red")
    }
    mtext(main, 3, line = 2, outer = TRUE, adj = adj.mtext,
          padj = padj.mtext, col = col.mtext, ...)
    mtext(sub, 1, line = 4, outer = TRUE, adj = adj.mtext,
          padj = padj.mtext, col = col.mtext, ...)
  }
}


#' plot mixtureVAR
#'
#' @param x irf.mixvar Object
#'
#' @export
plot.mixvar <- function (x, plot.type = "multiple", names = NULL,
                         main = NULL, sub = NULL, lty = NULL, lwd = NULL, col = NULL,
                         ylim = NULL, ylab = NULL, xlab = NULL, nc, mar.multi = c(0, 4, 0, 4), oma.multi = c(6, 4, 6, 4), adj.mtext = NA,
                         padj.mtext = NA, col.mtext = NA, path=NULL,...)
{
  op <- graphics::par(no.readonly = TRUE)
  on.exit(par(op))
  plot.type <- match.arg(plot.type)
  inames <- x$impulse
  rnames <- x$response
  if (is.null(names)) {
    names <- inames
  }
  else {
    names <- as.character(names)
    if (!(all(names %in% inames))) {
      warning("\nInvalid variable name(s) supplied, using first variable.\n")
      inames <- inames[1]
    }
    else {
      inames <- names
    }
  }
  nvi <- length(inames)
  nvr <- length(rnames)
  ifelse(is.null(lty), lty <- c(1, 1, 2, 2), lty <- rep(lty,
                                                        4)[1:4])
  ifelse(is.null(lwd), lwd <- c(1, 1, 1, 1), lwd <- rep(lwd,
                                                        4)[1:4])
  ifelse(is.null(col), col <- c("black", "gray", "red", "red"),
         col <- rep(col, 4)[1:4])

  if (plot.type == "multiple") {
    for(k in 1:length(x$irf)){
      if(is.null(path)) pdf(file=paste0("regime",k,".pdf")) else pdf(file=paste0(path,"/regime",k,".pdf"))
      xtmp<-x
      xtmp$irf <- x$irf[[k]]
      for (i in 1:nvi) {
        dp <- dataplot(xtmp, iname = inames[i])
        plot.multiple(dp, nc = nc, main=main, ylim= ylim, mar.multi = mar.multi, oma.multi = oma.multi,
                      ylab = ylab, lty = lty, lwd = lwd, col= col, adj.mtext = adj.mtext, padj.mtext = padj.mtext,
                      col.mtext = col.mtext, ...)
        if (nvi > 1) par(ask = TRUE)
      }
      dev.off()
    }
  }
}
#' Generate plots for mixtureVAR
#'
#' @param mod the estimated model by mixtureVAR
#' @param type type of plot desired
#' @param start
#' @param frequency
#'
#' @export
plots <- function(mod,type="residuals",start=1,frequency=1){
  K<-length(mod$Alpha)
  if(type=="tau"){
    plotdat<-ts(sapply(mod$Tau,"["),start=start,frequency=frequency)
    ts.plot(plotdat,col=1:K)
  } else if(type=="residuals") {
    par(mfrow=c(K,1))
    ylim <- range(unlist(mod$e))
    for(i in 1:K){
      n <- ncol(mod$e[[i]])
      plot(y=mod$e[[i]][1,],x=1:n,col=1,type="l",lty=1,
           ylab=paste("e of mod",i), xlab = "time", ylim=ylim)
      for(j in 1:nrow(mod$e[[i]]))
        lines(y=mod$e[[i]][j,],x=1:n,col=j,lty=j)
    }
  }
}

#' Generate plots for mixtureVAR
#'
#' @param mod the estimated model by mixtureVAR
#' @param type type of plot desired
#' @param start
#' @param frequency
#'
#' @export
plots2 <- function(mod,type="residuals",start=1,frequency=1){
  K<-length(mod$Alpha)
  if(type=="tau"){
    plotdat<-ts(sapply(mod$Tau,"["),start=start,frequency=frequency)
    ts.plot(plotdat,col=1:K)
  } else if(type=="residuals") {
    par(mfrow=c(K,1))
    ylim <- range(unlist(mod$e))
    for(i in 1:K){
      n <- ncol(mod$e[[i]])
      plot(y=mod$e[[i]][1,]*mod$Tau[[i]],x=1:n,col=1,type="l",lty=1,
           ylab=paste("e of mod",i), xlab = "time", ylim=ylim)
      for(j in 1:nrow(mod$e[[i]]))
        lines(y=mod$e[[i]][j,]*mod$Tau[[i]],x=1:n,col=j,lty=j)
    }
  }
}
