#' impulse response functions
#'
#' compute impulse response function for mixtureVARs
#'
#' @param x an estimates mixturVAR model
#' @param impulse the names of the impulse variables
#' @param response the names of the response variables
#' @param n.ahead the number of timesteps going ahead
#' @param ortho TRUE/FALSE Should the residuals be orthogonalized
#' @param cumulative TRUE/FALSE Should a cumulative irf be used
#' @param ci nominal coverage rate of the interval
#' @param seed starting seed
#'
#' @return list
#'
#' @export
irf.mixvar = function (x, impulse = NULL, response = NULL, n.ahead = 10, ortho = TRUE,
                       cumulative = FALSE, ci = 0.95, seed = NULL,
                       ...)
{

  y.names <- rownames(x$Theta[[1]]) ### Sichergehen, dass Datensatz columnwise gegeben ist?
  if (is.null(impulse)) {
    impulse <- y.names
  }
  else {
    impulse <- as.vector(as.character(impulse))
    if (any(!(impulse %in% y.names))) {
      stop("\nPlease provide variables names in impulse\nthat are in the set of endogenous variables.\n")
    }
    impulse <- subset(y.names, subset = y.names %in% impulse)
  }
  if (is.null(response)) {
    response <- y.names
  }
  else {
    response <- as.vector(as.character(response))
    if (any(!(response %in% y.names))) {
      stop("\nPlease provide variables names in response\nthat are in the set of endogenous variables.\n")
    }
    response <- subset(y.names, subset = y.names %in% response)
  }
  irs <- .irf.mixvar(x = x, impulse = impulse, response = response,
                     y.names = y.names, n.ahead = n.ahead, ortho = ortho,
                     cumulative = cumulative)
  Lower <- NULL
  Upper <- NULL
  result <- list(irf = irs, Lower = Lower, Upper = Upper, response = response,
                 impulse = impulse, ortho = ortho, cumulative = cumulative,
                 ci = ci, model = class(x))
  class(result) <- "varirf"
  return(result)
}

