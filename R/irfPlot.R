#' irfplot1
#'
#' @param DAT
#' @param grey
#' @param response
#' @param relation
#'
#' @export
irfplot1 <- function(DAT,grey=TRUE,response="all",relation="free"){
  LAT <- expand.grid(dimnames(DAT)[-2])
  LAT$lower1 <- as.vector(DAT[,1,])
  LAT$lower2 <- as.vector(DAT[,2,])
  LAT$mean   <- as.vector(DAT[,3,])
  LAT$upper2 <- as.vector(DAT[,4,])
  LAT$upper1 <- as.vector(DAT[,5,])
  if(response[1]=="all") response <- unique(LAT$variable)
  if(grey){
    mmcolors <- c("black","darkgrey")
    mmlty    <- c(1,1)
    mmpch    <- c(1,4)
  } else {
    mmcolors <- c("blue","black")
    mmlty    <- c(1,1)
    mmpch    <- c(1,4)
  }
  PLOT <- xyplot( mean ~ time| variable , data=LAT, subset = variable%in%response,
                  upper1 = LAT$upper1, lower1 = LAT$lower1,upper2 = LAT$upper2, lower2 = LAT$lower2, prepanel = mixtureVAR:::prepanel.default.xyplot2,
                  panel = function(x, y,upper1, lower1,upper2, lower2, subscripts, ...){
                    upper1 <- upper1[subscripts]
                    lower1 <- lower1[subscripts]
                    upper2 <- upper2[subscripts]
                    lower2 <- lower2[subscripts]
                    panel.polygon(c(x, rev(x)), c(upper1, rev(lower1)), fill=.5, col="lightgrey", border=FALSE,...)
                    panel.polygon(c(x, rev(x)), c(upper2, rev(lower2)), fill=.5, col="grey", border=FALSE,...)
                    panel.abline(h=0,col="red")
                    panel.xyplot(x, y, type='b', cex=0.2, lty=1, col=mmcolors[1], pch = mmpch[1],...)
                  },ylab="",xlab="",main="",
                  #main=paste("The impulse variable is",impulse,"with bands [",paste(range(bands),collapse = ";"),"]"),
                  #key=list(space="top", columns = 1,lines=list(col=mmcolors[1], type="b",  pch=mmpch[1], lty=mmlty[1], lwd=2),text=list(c("bootirf mean"))),
                  strip=strip.custom(  strip.names       = c(F,F),strip.levels=c(T,F)),
                  strip.left=FALSE)
  #library(latticeExtra)
  n.ahead <- dim(DAT)[1]
  at.seq <- seq(0,n.ahead,by = 3)+1
  at.lab <- at.seq-1
  at.lab[!((1:length(at.lab))%%2)] <- ""
  at.lab[1] <- ""
  PLOT <- update(PLOT, scales = list(y=list(relation=relation,alternating=3,rot=0),
                                     x=list(alternating=1,at=at.seq,labels=at.lab, rot = 0, tck=c(1,0))),
                 par.settings=list(strip.background = list( col="transparent" ),
                                   strip.border     = list( col = NA )))
  #if(length(response)==1 & response!="all") useOuterStrips(PLOT,strip.left=FALSE) else useOuterStrips(PLOT)
  PLOT
}
