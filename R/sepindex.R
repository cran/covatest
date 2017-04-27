#' Class "sepindex"
#'
#' A class for the non-separability index (r) for different spatial
#' and temporal lags:
#' \deqn{r(h, u, \Theta)= \rho(h, u;\Theta)/ [\rho(h,0;\Theta)\rho(0,u;\Theta)]}
#' with \eqn{\rho(h, u;\Theta)>0}; \eqn{\rho(h,0;\Theta)>0} and \eqn{\rho(0,u;\Theta)>0}.
#' On the basis of this index, the type of non-separability of the covariance
#' function can be analyzed.
#'
#' @slot sep.index.ratio the empirical non-separability index ratio
#' @slot cov.tm the purely temporal sample covariance function
#' @slot cov.sp the purely spatial sample covariance function
#' @slot cov.st the spatio-temporal sample covariance function
#'
#' @rdname sepindex-class
#' @exportClass sepindex
setClass("sepindex", slots = c(sep.index.ratio = "matrix",
                               cov.st = "matrix",
                               cov.tm = "matrix",
                               cov.sp = "matrix"))

#' @param nt integer, the number of temporal lags in \code{vario_st}
#' @param ns integer, the number of spatial lags in \code{vario_st}
#' @param vario_st spatio-temporal sample variogram, output from
#' \code{\link[gstat]{variogramST}}
#' @param globalSill numeric, the value of the sample variance
#'
#' @references De Iaco, S., Posa, D., 2013, Positive and negative non-separability
#' for space-time covariance models. Journal of Statistical Planning and
#' Inference, \bold{143} 378--391.
#'
#' Pebesma, E., 2004, Multivariable geostatistics in S: the gstat
#' package. Computers & Geosciences, \bold{30} 683--691.
#'
#' Rodriguez, A., Diggle, P.J., 2010, A class of convolution-based
#' models for spatio-temporal processes with non-separable covariance structure.
#' Scandinavian Journal of Statistics, \bold{37(4)} 553--567.
#'
#' @seealso \code{\link[gstat]{variogramST}}
#' @importFrom methods new
#' @importFrom graphics boxplot
#' @usage sepindex(nt, ns, vario_st, globalSill)
#' @examples
#' library(covatest)
#' data(rr_13)
#' data(vv_13)
#' #compute the globalSill
#' C00_13 <- var(rr_13[, ,"PM10"]@data[[1]], na.rm = TRUE)
#' nonsepind <- sepindex(nt = 16, ns = 4, vario_st = vv_13, globalSill = C00_13)
#'
#' @rdname sepindex-class
#' @export
sepindex <- function(nt, ns, vario_st, globalSill) {
  # Set 0 to NA data in vario_st
  vario_st$dist[1] <- 0
  vario_st$gamma[1] <- 0

  # Compute ST covariance#
  cov_st <- vario_st
  cov_st$gamma <- globalSill - cov_st$gamma
  cov_st <- as.matrix(cbind(cov_st$avgDist, cov_st$gamma, cov_st$timelag))
  colnames(cov_st) <- c("avgDist", "C_st", "timelag")

  # cov_st<<-cov_st to save the object in the global environment

  # compute marginal temporal covariances cov_tm
  n <- nt * ns
  cov_tm <- matrix(0, ncol = 3, nrow = nt)
  j <- 1
  for (i in 1:nt) {
    cov_tm[i, ] <- cov_st[j, ]
    j <- j + ns

  }

  cov_tm <- cbind(rep(cov_tm[, 3], each = ns), rep(cov_tm[, 2], each = ns))
  colnames(cov_tm) <- c("ht", "Ct")

  # cov_tm<<-cov_tm to save the object in the global environment

  # compute marginal spatial covariances cov_sp
  cov_sp <- cov_st[1:ns, ]
  cov_sp <- cbind(rep(cov_sp[, 1], times = nt), rep(cov_sp[, 2], times = nt))
  colnames(cov_sp) <- c("hs", "C_sp")

  # cov_sp<<-cov_sp to save the object in the global environment

  # Intermediate computation for the non-separability index cov_sep
  cov_sep <- cbind(cov_sp[, 1], cov_tm[, 1], (cov_sp[, 2] * cov_tm[,
                                                              2]) / globalSill)
  colnames(cov_sep) <- c("ht", "hs", "C'_st")

  # cov_sep<<-cov_sep to save the object in the global environment

  # Compute the non-separability index
  sep_index_ratio <- cbind(cov_sp[, 1], cov_tm[, 1], cov_st[, 2] / cov_sep[,
                                                                         3])
  sep_index_ratio <- round(sep_index_ratio, digits = 3)
  colnames(sep_index_ratio) <- c("hs", "ht", "SepIndex")
  x <- sep_index_ratio[, 3]
  y <- x[x < 0]
  z <- length(y) / n * 100
  if (z > 2) {
    warning(round(z, digits = 3), "%  of non admissible values for the non-separability index (<0) have been removed")
  }
  sep_index_ratio_pos <- sep_index_ratio[sep_index_ratio[, 3] >= 0, ]

  new("sepindex", sep.index.ratio = sep_index_ratio_pos, cov.st = cov_st,
     cov.tm = cov_tm, cov.sp = cov_sp)

}

#' @param x object of class \code{sepindex}
#' @param ... any arguments that will be passed to the panel plotting functions
#' @rdname sepindex-class
#' @aliases sepindex-class
#' @aliases sepindex-method
#' @aliases plot
setMethod("boxplot", signature = c(x = "sepindex"),
          function(x, ...) {
            boxplot(x@sep.index.ratio[, 3] ~ x@sep.index.ratio[, 1],
                    data = x@sep.index.ratio,
                    ...)
            boxplot(x@sep.index.ratio[, 3] ~ x@sep.index.ratio[, 2],
                    data = x@sep.index.ratio,
                    ...)
          }

)
