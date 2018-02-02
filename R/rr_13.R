#' Air quality data, rural backgroung PM10 in Germany, daily averages 2005-2006
#'
#' A subset of the \code{STFDF} \code{r5to10} used in the vignette of the package
#' \code{\link{gstat}}: "Introduction to Spatio-Temporal Variography".
#' This STFDF consists of PM10 daily data, measured from the 1st January 2005 to
#' the 31st December 2006 in 13 rural background stations, located in the central
#' Germany area.
#'
#' @docType data
#'
#' @usage data(rr_13)
#'
#' @format Object of class \code{\link[spacetime]{STFDF}} with dimensions
#' (s, t, attr): (13, 730, 1).
#'
#' @keywords dataset
#'
#' @note The primary dataset on air quality data has been compiled for R by
#' Benedict Graeler in \code{spacetime}
#'
#' @seealso \code{\link[spacetime]{air}}
#'
#' @references
#' \url{http://acm.eionet.europa.eu/databases/airbase}
#'
#' Graeler, B., Pebesma, E., Heuvelink, G., 2016. Spatio-Temporal Interpolation
#' using gstat. The R Journal \bold{8(1)} 204--218
#'
#' @examples
#' ## Not run:
#' # rr_13 is obtained by running the following command lines:
#' library(sp)
#' library(spacetime)
#' library(gstat)
#' rm(list = ls())
#' data(air)
#' ls()
#' if (!exists("rural")) rural = STFDF(stations, dates, data.frame(PM10 =
#' as.vector(air)))
#' rr = rural[,"2005::2010"]
#' unsel = which(apply(as(rr, "xts"), 2, function(x) all(is.na(x))))
#' r5to10 = rr[-unsel,]
#' rr_13 <- r5to10[c("DEHE046","DESN049","DETH026","DENW063","DETH061","DEBY047",
#' "DENW065","DEUB029","DENW068","DENI019","DEHE051","DERP016","DENI051"),
#' "2005::2006"]
#' ## End (Not run)
"rr_13"
