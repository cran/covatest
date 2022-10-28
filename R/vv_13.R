#' Precomputed variogram for PM10 in data set air
#'
#' Precomputed variogram for PM10 in a subset of the air quality data set
#' \code{\link[spacetime]{air}}
#'
#' @docType data
#'
#' @usage data(vv_13)
#'
#' @format Object of class StVariogram
#'
#' @keywords dataset StVariogram
#'
#' @examples
#' # --start define the STFDF rr_13-- #
#' library(sp)
#' library(spacetime)
#'
#'
#' data(air)
#'
#' ls()
#'
#' if (!exists("rural")) rural = STFDF(stations, dates, data.frame(PM10 =
#' as.vector(air)))
#'
#' rr <- rural[,"2005::2010"]
#'
#' unsel = which(apply(as(rr, "xts"), 2, function(x) all(is.na(x))))
#'
#' r5to10 = rr[-unsel,]
#'
#' rr_13 <- r5to10[c("DEHE046","DESN049","DETH026","DENW063","DETH061","DEBY047",
#' "DENW065","DEUB029","DENW068","DENI019","DEHE051","DERP016","DENI051"),
#' "2005::2006"]
#' # --end define the STFDF rr_13-- #
#'
#' ## Not run
#' ## To estimate the spatio-temporal variogram, paste and copy the following lines
#' ## (without the symbol '#') in the console:
#' #
#' ## vv_13 is obtained by running the function variogramST of the package gstat,
#' ## as follows
#' #
#' # vv_13 <- gstat::variogramST(PM10~1, rr_13, width=60, cutoff = 220, tlags=0:15)
#' ## End (Not run)
"vv_13"
