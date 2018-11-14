### State-Dependent Transmission of Monetary Policy in the Euro Area
### by Jan Pablo Burgard, Matthias Neuenkirch, and Matthias Nöckel

### For questions concerning the R package, please contact Jan Pablo Burgard (burgardj@uni-trier.de).
### For questions concerning the article, please contact Matthias Neuenkirch (neuenkirch@uni-trier.de).

### Execute script step by step
### All other files are auxiliary files

mVAR <- function() {

  ### Seting working directory ###
  library(mixtureVAR)
  rm(list=ls())
  setwd("D:/Users/neuenkirch/Desktop/Dropbox/Forschung/State-Dependent Transmission/Package JMCB/mixtureVAR")

  ### Loading data ###
  Dat <<- read.csv("data/ECB.csv", sep=";", header=T)

  ### Detrending and creating lagged variables ###
  TREND <<- 1:204

  det_ip <<- lm(Dat$IP ~ TREND)
  summary(det_ip)
  Dat$IP <- residuals(det_ip)*100

  det_infl <<- lm(Dat$INFL ~ TREND)
  summary(det_infl)
  Dat$INFL <- residuals(det_infl)

  det_mrr <<- lm(Dat$MRR_SHADOW ~ TREND)
  summary(det_mrr)
  Dat$MRR_SHADOW <- residuals(det_mrr)

  det_vstoxx <<- lm(Dat$VSTOXX ~ TREND)
  summary(det_vstoxx)
  Dat$VSTOXX <- residuals(det_vstoxx)

  Dat$IP_lag <-NA
  Dat$IP_lag[(2):nrow(Dat)]  <- Dat$IP[1:(nrow(Dat)-1)]

  Dat$INFL_lag <-NA
  Dat$INFL_lag[(2):nrow(Dat)]  <- Dat$INFL[1:(nrow(Dat)-1)]

  Dat$MRR_SHADOW_lag <-NA
  Dat$MRR_SHADOW_lag[(2):nrow(Dat)]  <- Dat$MRR_SHADOW[1:(nrow(Dat)-1)]

  Dat$VSTOXX_lag <-NA
  Dat$VSTOXX_lag[(2):nrow(Dat)]  <- Dat$VSTOXX[1:(nrow(Dat)-1)]

  rm(det_infl, det_ip, det_mrr, det_vstoxx, TREND, pos = ".GlobalEnv")

  ### Time series plots (Figure A1)
  xticks <<- seq(1999,2016,1)

  plot(ts(Dat[,2],start=c(1999,1),frequency=12),
       main = "Log(Industrial Production)",lwd=2,xlab='',ylab='',xaxt = 'n')
  axis(1, at = xticks, labels = xticks)

  plot(ts(Dat[,3],start=c(1999,1),frequency=12),
       main = "Inflation",lwd=2,xlab='',ylab='',xaxt = 'n')
  axis(1, at = xticks, labels = xticks)

  plot(ts(Dat[,4],start=c(1999,1),frequency=12),
       main = "Interest Rate",lwd=2,xlab='',ylab='',xaxt = 'n')
  axis(1, at = xticks, labels = xticks)

  plot(ts(Dat[,5],start=c(1999,1),frequency=12),
       main = "VSTOXX",lwd=2,xlab='',ylab='',xaxt = 'n')
  axis(1, at = xticks, labels = xticks)

  ### Estimate MixtureVAR with four lags in both states and lagged latent variable in submodel ###
  YFormula4      <<- list(formula( ~ 0 + IP + INFL + MRR_SHADOW + VSTOXX),
                         formula( ~ 0 + IP + INFL + MRR_SHADOW + VSTOXX))
  LatentFormula4 <<- formula( ~ IP_lag + INFL_lag + MRR_SHADOW_lag + VSTOXX_lag,  LatentLag = c(1))

  set.seed(1)
  mod.mixres4 <<- mixtureVAR::mixtureVAR(YFormula=YFormula4, LatentFormula = LatentFormula4,
                  data = Dat, P = list(4,4), nreps = 2, maxiter=500, eps=1e-04)

  ### Residual analysis (including Figure A2) ###
  bestmod4 <<- mod.mixres4$finalmodel
  we4 <<- mapply(function(e,w)t(e)*w,bestmod4$e,bestmod4$Tau,SIMPLIFY = FALSE)
  wes4 <<- Reduce("+",we4)
  wes.ts4 <<-ts(wes4,start=c(1999,5),frequency=12)

  Box.test(wes.ts4[,1], lag = 6, type = "Ljung-Box", fitdf = 0)
  Box.test(wes.ts4[,2], lag = 6, type = "Ljung-Box", fitdf = 0)
  Box.test(wes.ts4[,3], lag = 6, type = "Ljung-Box", fitdf = 0)
  Box.test(wes.ts4[,4], lag = 6, type = "Ljung-Box", fitdf = 0)

  yticks_inf    <<- seq(-0.75,0.75,0.25)
  yticks_vstoxx <<- seq(-20,20,5)

  plot(ts(wes.ts4[,1],start=c(1999,3),frequency=12),
       main = "Log(Industrial Production): Residuals",
       lwd=2,xlab='',ylab='',xaxt = 'n',ylim=c(-2.5,2.5))
  abline(h = 0, col="red")
  axis(1, at = xticks, labels = xticks)

  plot(ts(wes.ts4[,2],start=c(1999,3),frequency=12),
       main = "Inflation: Residuals",
       lwd=2,xlab='',ylab='',yaxt = 'n',xaxt = 'n',ylim=c(-0.75,0.75))
  abline(h = 0, col="red")
  axis(1, at = xticks, labels = xticks)
  axis(2, at = yticks_inf, labels = yticks_inf)

  plot(ts(wes.ts4[,3],start=c(1999,3),frequency=12),
       main = "Interest Rate: Residuals",
       lwd=2,xlab='',ylab='',xaxt = 'n', ylim=c(-1,1))
  abline(h = 0, col="red")
  axis(1, at = xticks, labels = xticks)

  plot(ts(wes.ts4[,4],start=c(1999,3),frequency=12),
       main = "VSTOXX: Residuals",
       lwd=2,xlab='',ylab='',yaxt = 'n',xaxt = 'n',ylim=c(-20,20))
  abline(h = 0, col="red")
  axis(1, at = xticks, labels = xticks)
  axis(2, at = yticks_vstoxx, labels = yticks_vstoxx)

  ### State weights and correlations (Figure 1) ###
  ### States 1 and 2 might be flipped compared to paper (see also below) ###
  ### "Stronger" state represents "normal" times ###
  ### Switching state labels in figures might be necessary: stronger state is State 1 in paper ###
  States <- unlist(bestmod4$Tau)
  dim(States) <- c(200, 2)
  colnames(States) <- c("State 1","State 2")

  plot(ts(States[,1],start=c(1999,5),frequency=12),
       main = "Weight of State 1 (Normal)",
       lwd=2,xlab='',ylab='',xaxt = 'n',ylim=c(0,1))
  axis(1, at = xticks, labels = xticks)

  plot(ts(States[,2],start=c(1999,5),frequency=12),
       main = "Weight of State 2 (Crisis)",
       lwd=2,xlab='',ylab='',xaxt = 'n',ylim=c(0,1))
  axis(1, at = xticks, labels = xticks)

  cor(Dat[5:204,6],States[,2]) # IP_lag
  cor(Dat[5:204,7],States[,2]) # INFL_lag
  cor(Dat[5:204,8],States[,2]) # MRR_SHADOW_lag
  cor(Dat[5:204,9],States[,2]) # VSTOXX_lag

  ### Again, check state ordering: "stronger" state represents "normal" times ###
  bestmod4$Alpha
  ### Switching state labels in figures might be necessary: stronger state is State 1 in paper ###

  LatentFormula <<- LatentFormula4

  ### Predicted probabilities of logit submodel (Figure 2)
  mixtureVAR::mnlAveEffPlot2(bestmod4, "IP_lag"         ,
                             ColNames = c("Lagged Industr. Prod. (State 1: Normal)",
                             "Lagged Industr. Prod. (State 2: Crisis)"),R = 1500, nvals = 25, plot = TRUE)
  SummaryOutput$mean[1]  ### Predicted probability of State 1 at minimum of variable
  SummaryOutput$mean[49] ### Predicted probability of State 1 at maximum of variable

  mixtureVAR::mnlAveEffPlot2(bestmod4, "INFL_lag"       ,
                             ColNames = c("Lagged Inflation (State 1: Normal)",
                             "Lagged Inflation (State 2: Crisis)"),R = 1500, nvals = 25, plot = TRUE)
  SummaryOutput$mean[1]
  SummaryOutput$mean[49]

  mixtureVAR::mnlAveEffPlot2(bestmod4, "MRR_SHADOW_lag" ,
                             ColNames = c("Lagged Interest Rate (State 1: Normal)",
                             "Lagged Interest Rate (State 2: Crisis)"),R = 1500, nvals = 25, plot = TRUE)
  SummaryOutput$mean[1]
  SummaryOutput$mean[49]

  mixtureVAR::mnlAveEffPlot2(bestmod4, "VSTOXX_lag"     ,
                             ColNames = c("Lagged VSTOXX (State 1: Normal)",
                             "Lagged VSTOXX (State 2: Crisis)"),R = 1500, nvals = 25, plot = TRUE)
  SummaryOutput$mean[1]
  SummaryOutput$mean[49]

  ### Impulse responses for baseline ordering (Figure 3) ###
  ### Again, check state ordering: "stronger" state represents "normal" times ###
  bestmod4$Alpha
  ### Switching state labels in figures might be necessary: stronger state is State 1 in paper ###
  ### For alternative ordering (Figure A3) interchange interest rate indicator and VSTOXX
  bestmod4$maxiter <- 500
  set.seed(1)
  BOOTAUD2  <<- mixtureVAR::Boot_mixtureVAR(bestmod4,nboot = 10, eps = 1e-4)

  mixtureVAR::plotBootIrf(BOOTAUD2,bestmod4,impulse="MRR_SHADOW",response="IP",
              n.ahead=48,type="mean",relation="same",grey="TRUE",
              modnames =c("Industr. Prod. (State 1: Normal)", "Industr. Prod. (State 2: Crisis)"))

  ### Exporting summary table for all impulse responses ###
  writexl::write_xlsx(SummaryOutput, paste("Summary.xlsx", sep=""))

  mixtureVAR::plotBootIrf(BOOTAUD2,bestmod4,impulse="MRR_SHADOW",response="INFL",
              n.ahead=48,type="mean",relation="same",grey="TRUE",
              modnames =c("Inflation (State 1: Normal)", "Inflation (State 2: Crisis)"))

  mixtureVAR::plotBootIrf(BOOTAUD2,bestmod4,impulse="MRR_SHADOW",response="MRR_SHADOW",
              n.ahead=48,type="mean",relation="same",grey="TRUE",
              modnames =c("Interest Rate (State 1: Normal)", "Interest Rate (State 2: Crisis)"))

  mixtureVAR::plotBootIrf(BOOTAUD2,bestmod4,impulse="MRR_SHADOW",response="VSTOXX",
              n.ahead=48,type="mean",relation="same",grey="TRUE",
              modnames =c("VSTOXX (State 1: Normal)", "VSTOXX (State 2: Crisis)"))
}

LVARWCI <- function() {
  ### Linear VAR with Cholesky Identification (Figure 4) ###
  Dat <<- Dat[,2:5]
  VARMOD1 <<- vars::VAR(Dat,lag=4)

  IRF1_IP95     <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="IP",         ci = 0.95, runs = 500)
  IRF1_IP68     <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="IP",         ci = 0.68, runs = 500)
  IRF1_INFL95   <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="INFL",       ci = 0.95, runs = 500)
  IRF1_INFL68   <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="INFL",       ci = 0.68, runs = 500)
  IRF1_IR95     <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="MRR_SHADOW", ci = 0.95, runs = 500)
  IRF1_IR68     <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="MRR_SHADOW", ci = 0.68, runs = 500)
  IRF1_VSTOXX95 <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="VSTOXX",     ci = 0.95, runs = 500)
  IRF1_VSTOXX68 <<- irf(VARMOD1,impulse="MRR_SHADOW",n.ahead = 48, response="VSTOXX",     ci = 0.68, runs = 500)

  ### Adjust shock to weighted shock mixture VAR ###
  SHOCK        <<- IRF1_IR95$irf$MRR_SHADOW[1]
  SHOCK_FACTOR <<- 0.406032542202903/SHOCK

  IRF1_IP     <<- cbind(IRF1_IP95$Lower$MRR_SHADOW     * SHOCK_FACTOR,
                        IRF1_IP68$Lower$MRR_SHADOW     * SHOCK_FACTOR,
                        IRF1_IP95$irf$MRR_SHADOW       * SHOCK_FACTOR,
                        IRF1_IP68$Upper$MRR_SHADOW     * SHOCK_FACTOR,
                        IRF1_IP95$Upper$MRR_SHADOW     * SHOCK_FACTOR)
  IRF1_INFL   <<- cbind(IRF1_INFL95$Lower$MRR_SHADOW   * SHOCK_FACTOR,
                        IRF1_INFL68$Lower$MRR_SHADOW   * SHOCK_FACTOR,
                        IRF1_INFL95$irf$MRR_SHADOW     * SHOCK_FACTOR,
                        IRF1_INFL68$Upper$MRR_SHADOW   * SHOCK_FACTOR,
                        IRF1_INFL95$Upper$MRR_SHADOW   * SHOCK_FACTOR)
  IRF1_IR     <<- cbind(IRF1_IR95$Lower$MRR_SHADOW     * SHOCK_FACTOR,
                        IRF1_IR68$Lower$MRR_SHADOW     * SHOCK_FACTOR,
                        IRF1_IR95$irf$MRR_SHADOW       * SHOCK_FACTOR,
                        IRF1_IR68$Upper$MRR_SHADOW     * SHOCK_FACTOR,
                        IRF1_IR95$Upper$MRR_SHADOW     * SHOCK_FACTOR)
  IRF1_VSTOXX <<- cbind(IRF1_VSTOXX95$Lower$MRR_SHADOW * SHOCK_FACTOR,
                        IRF1_VSTOXX68$Lower$MRR_SHADOW * SHOCK_FACTOR,
                        IRF1_VSTOXX95$irf$MRR_SHADOW   * SHOCK_FACTOR,
                        IRF1_VSTOXX68$Upper$MRR_SHADOW * SHOCK_FACTOR,
                        IRF1_VSTOXX95$Upper$MRR_SHADOW * SHOCK_FACTOR)

  IRF1 <<- abind::abind('Industrial Production'=IRF1_IP,'Inflation'=IRF1_INFL,along=3)
  dimnames(IRF1) <- list(time=1:49,value=c("lower1","lower2","mean","upper2","upper1"),
                          variable=dimnames(IRF1)[[3]])

  IRF2 <<- abind::abind('Interest Rate'=IRF1_IR,'VSTOXX'=IRF1_VSTOXX,along=3)
  dimnames(IRF2) <- list(time=1:49,value=c("lower1","lower2","mean","upper2","upper1"),
                          variable=dimnames(IRF2)[[3]])

  mixtureVAR::irfplot1(IRF1,grey=TRUE,response = "Industrial Production")
  mixtureVAR::irfplot1(IRF1,grey=TRUE,response = "Inflation")
  mixtureVAR::irfplot1(IRF2,grey=TRUE,response = "Interest Rate")
  mixtureVAR::irfplot1(IRF2,grey=TRUE,response = "VSTOXX")


  ### Comparison to other non-linear VARs (Figure 5) ###
  regimes_nonlinear <<- read.csv("data/regimes_nonlinear.csv", sep=";", header=T)

  ### Again, check state ordering: "stronger" state represents "normal" times ###
  ### Hence, relying on State 2 in figures might be necessary ###
  bestmod4$Alpha
  state1_mixtureVAR <- ts(States[,1],              start = c(1999,5),frequency=12)
  regime1_msvar     <- ts(regimes_nonlinear$MSVAR, start = c(1999,5),frequency=12)
  regime1_lstvar    <- ts(regimes_nonlinear$LSTVAR,start = c(1999,5),frequency=12)

  ts.plot(ts(regime1_msvar,start = c(1999,5),frequency=12),state1_mixtureVAR,
          main="MSVAR vs Logit-MVAR", gpars=list(xlab="", ylab="",yaxt = 'n',xaxt = 'n',
          lwd=2, ylim=c(0,1),lty=c("solid","dashed")))
  axis(1, at = xticks, labels = xticks)
  axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1))

  ts.plot(ts(regime1_lstvar,start = c(1999,5),frequency=12),state1_mixtureVAR,
          main="LSTVAR vs Logit-MVAR", gpars=list(xlab="", ylab="",yaxt = 'n',xaxt = 'n',
          lwd=2, ylim=c(0,1), lty=c("solid","dashed")))
  axis(1, at = xticks, labels = xticks)
  axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1))
}
