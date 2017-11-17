#' Precomputed variogram for PM10 in data set rr_13
#'
#' Precomputed variogram for PM10 in data set rr_13.
#'
#' @docType data
#'
#' @usage data(vv_13)
#'
#' @format Object of class StVariogram
#'
#' @keywords dataset, StVariogram
#'
#' @examples
#' ## To run the example, paste and copy the following lines
#' # (without the symbol '#') in the console:
#' #
#' # rr_13 is obtained by running the following command lines:
#' #
#' # library(spacetime)
#' # library(gstat)
#' # rm(list = ls())
#' # data(rr_13)
#' # vv_13 = variogram(PM10~1, rr_13, width=60, cutoff = 220, tlags=0:15)
#' ## End (Not run)
"vv_13"
