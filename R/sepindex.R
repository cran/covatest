#' Class "sepindex"
#'
#' A class for the non-separability index (r) for different spatial
#' and temporal lags:
#' \deqn{r(h, u, \Theta)= \rho(h, u;\Theta)/ [\rho(h,0;\Theta)\rho(0,u;\Theta)]}
#' with \eqn{\rho(h, u;\Theta)>0}; \eqn{\rho(h,0;\Theta)>0} and \eqn{\rho(0,u;\Theta)>0}.
#' On the basis of this index, the type of non-separability of the covariance
#' function can be analyzed.
#'
#' @slot sep.index.ratio the empirical non-separability index ratio and the
#' corresponding spatio-temporal lags
#' @slot cov.st the spatio-temporal sample covariance function and the
#' corresponding spatio-temporal lags
#' @slot cov.tm the purely temporal sample covariance function and the
#' corresponding temporal lags
#' @slot cov.sp the purely spatial sample covariance function and the
#' corresponding spatial lags
#'
#' @rdname sepindex-class
#' @exportClass sepindex
setClass("sepindex", slots = c(sep.index.ratio = "matrix",
                               cov.st = "matrix",
                               cov.tm = "matrix",
                               cov.sp = "matrix"))

#' @param vario_st spatio-temporal sample variogram, output from
#' \code{\link[gstat]{variogramST}}
#' @param nt integer, the number of temporal lags in \code{vario_st}
#' @param ns integer, the number of spatial lags in \code{vario_st}
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
#' @import stats
#' @import methods
#' @importFrom methods new
#' @importFrom graphics boxplot
#'
#' @examples
#' # --start define the STFDF rr_13-- #
#' library(sp)
#' library(spacetime)
#' library(gstat)
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
#' # --end define the STFDF rr_13-- #
#'
#' #compute the Global Sill
#' C00_13<-var(rr_13[,,"PM10"]@data[[1]], na.rm = TRUE)
#'
#' #estimate the spatio-temporal variogram
#' data(vv_13)
#' nonsep.index<-sepindex(vario_st = vv_13, nt = 16, ns = 4, globalSill = C00_13)
#'
#' ##methods for sepindex
#'
#' #1. show
#' nonsep.index
#'
#' #2. summary
#' summary(nonsep.index)
#'
#' #3. boxplot
#' boxplot(nonsep.index, ylab="Non-separability ratio")
#'
#' #4. [ extract
#' nonsep.index[1:8, ] #selection of the first 8 rows
#' nonsep.index[1:8, 1:2] #selection of the first 2 columns
#' @rdname sepindex-class
#' @export
sepindex <- function(vario_st, nt, ns, globalSill) {

  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}

  ### SOME CHECKS ON THE ARGUMENTS ###

  if (is.scalar(nt) == FALSE || is.scalar(ns) == FALSE  || is.scalar(globalSill) == FALSE) {
    stop("Some of the arguments are not numeric. Stop running.")
  }


  if (!inherits(vario_st, "StVariogram") && !inherits(vario_st, "data.frame")){
    stop("vario_st argument has to be of class StVariogram or data.frame.")
  }

  if(nt != as.integer(nt) || ns != as.integer(ns)){
    nt <- as.integer(nt)
    ns <- as.integer(ns)
    warning("nt and ns arguments are forced to be integer numbers.")
  }

  # Set 0 to NA data in vario_st
  vario_st$dist[1] <- 0
  vario_st$gamma[1] <- 0

  # Compute ST covariance#
  cov_st <- vario_st
  cov_st$gamma <- globalSill - cov_st$gamma
  cov_st2 <- as.matrix(cbind(cov_st$avgDist, cov_st$gamma, cov_st$timelag))
  cov_st <- cov_st2
  cov_st[which(cov_st < 0, arr.ind = TRUE, useNames = FALSE)] <- NA
  colnames(cov_st) <- c("avgDist", "C_st", "timelag")

  # compute marginal temporal covariances cov_tm
  n <- nt * ns
  cov_tm <- matrix(0, ncol = 3, nrow = nt)
  j <- 1
  for (i in 1:nt) {
    cov_tm[i, ] <- cov_st[j, ]
    j <- j + ns

  }
  cov_tm2 <- cov_tm[,2]
  cov_tm3 <- cov_tm[,3]
  cov_tm3 <- cbind(cov_tm3, cov_tm2)

  cov_tm <- cbind(rep(cov_tm[, 3], each = ns), rep(cov_tm[, 2], each = ns))
  colnames(cov_tm) <- c("ht", "Ct")
  colnames(cov_tm3) <- c("ht", "Ct")

  # compute marginal spatial covariances cov_sp
  cov_sp <- cov_st[1:ns, ]
  cov_sp2 <- cov_st[1:ns, -3]
  cov_sp <- cbind(rep(cov_sp[, 1], times = nt), rep(cov_sp[, 2], times = nt))
  colnames(cov_sp) <- c("hs", "C_sp")

  # Intermediate computation for the non-separability index cov_sep
  cov_sep <- cbind(cov_sp[, 1], cov_tm[, 1], (cov_sp[, 2] * cov_tm[,
                                                                   2]) / globalSill)
  colnames(cov_sep) <- c("ht", "hs", "C'_st")

  # Compute the non-separability index
  sep_index_ratio <- cbind(cov_sp[, 1], cov_tm[, 1], cov_st[, 2] / cov_sep[,
                                                                           3])
  sep_index_ratio <- round(sep_index_ratio, digits = 3)
  colnames(sep_index_ratio) <- c("hs", "ht", "SepIndex")
  x <- sep_index_ratio[, 3]
  y <- x[is.na(x)]
  z <- length(y) / n * 100
  if (z > 2) {
    warning(round(z, digits = 3), "%  of non admissible values for the non-separability index (<0) have been removed")
  }

  new("sepindex", sep.index.ratio = sep_index_ratio, cov.st = cov_st,
      cov.tm = cov_tm3, cov.sp = cov_sp2)

}
#' @param x object of class \code{sepindex} for methods \code{boxplot} and \code{extract}
#' @param ... any arguments that will be passed to the panel plotting functions
#' @rdname sepindex-class
#' @aliases sepindex-method
#' @aliases boxplot
#' @export
setMethod("boxplot", signature = c(x = "sepindex"),
          function(x, ...) {
            sep.index.ratio.pos <- x@sep.index.ratio[!is.na(x@sep.index.ratio[, 3]), ]

            boxplot(sep.index.ratio.pos[, 3] ~ sep.index.ratio.pos[, 1],
                    data = sep.index.ratio.pos, xlab="Spatial lags", ...)
            boxplot(sep.index.ratio.pos[, 3] ~ sep.index.ratio.pos[, 2],
                    data = sep.index.ratio.pos, xlab="Temporal lags", ...)
          }

)
#' @param object object of class \code{sepindex} for methods \code{show} and \code{summary}
#' @rdname sepindex-class
#' @aliases sepindex-class
#' @aliases sepindex-method
#' @aliases show
#' @export
setMethod(f="show", signature="sepindex", definition=function(object) {
  cat("An object of class sepindex, with", "\n")
  cat("nt=", nrow(object@cov.tm), "\n")
  cat("ns=", nrow(object@cov.sp), "\n")
  globalSill <- object@cov.st[1,2]
  cat("globalSill=", globalSill, "\n")
  cat("\n")
  cat("Slot 'sep.index.ratio':")
  cat("\n")
  print(object@sep.index.ratio)
  cat("\n")
  cat("Slot 'cov.st':")
  cat("\n")
  print(object@cov.st)
  cat("\n")
  cat("Slot 'cov.tm':")
  cat("\n")
  print(object@cov.tm)
  cat("\n")
  cat("Slot 'cov.sp':")
  cat("\n")
  print(object@cov.sp)
}
)
#' @param i index specifing elements to extract. Each row includes data for specific
#' spatio-temporal lags
#' @param j index specifing elements to extract. Set \code{1} for spatial lags (hs),
#' \code{2} for temporal lags (ht) and \code{3} for the non-separability index
#' (SepIndex)
#' @rdname sepindex-class
#' @aliases sepindex-class
#' @aliases sepindex-method
#' @aliases select
#' @export
setMethod(f="[", signature="sepindex", definition=function(x, i, j) {
  x@sep.index.ratio[i,j]
}
)

#' @rdname sepindex-class
#' @aliases sepindex-class
#' @aliases sepindex-method
#' @aliases summary
#' @export
setMethod(f = "summary", signature = "sepindex",
          definition = function(object) {

            x <- object@sep.index.ratio[,3]
            if(all(is.na(x)) == FALSE){
              ns <- nrow(object@cov.sp)
              nt <- nrow(object@cov.tm)

              min1 <- min(x, na.rm = TRUE)
              Q1 <- quantile(x, probs = 0.25, na.rm = TRUE)
              median1 <- median(x, na.rm = TRUE)
              mean1 <- mean(x, na.rm = TRUE)
              Q3 <- quantile(x, probs = 0.75, na.rm = TRUE)
              max1 <- max(x, na.rm = TRUE)
              sd1 <- sd(x, na.rm = TRUE)
              cat("\n")
              cat("Global summary", "\n")
              cat("\n")
              cat("Min. =" , min1 , "\n")
              cat("1st Qu. =" , Q1 , "\n")
              cat("Median =" , median1 , "\n")
              cat("Mean =" , mean1, "\n")
              cat("3rd Qu. =" , Q3 , "\n")
              cat("Max. =" , max1 , "\n")
              cat("St. Dev. =" , sd1, "\n")


              cat("\n")
              cat("Spatial summary", "\n")
              cat("\n")
              x <- matrix(data = 0, nrow = nt, ncol = 1)
              y <- matrix(data = 0, nrow = ns, ncol = 8)
              colnames(y) <- c(" Min. ", "  Q1  ", "Median", " Mean ", "  Q3  ", " Max. ", "  #>1  ", "  #<1  ")
              for(i in 1:ns){
                ratio.sp <- 0
                npos <- 0
                nneg <- 0

                for(j in 1:nt){
                  if(!is.na(object@sep.index.ratio[i + (j - 1) * ns, 3]) == TRUE){
                    if(object@sep.index.ratio[i + (j - 1) * ns, 3] > 1){
                      npos <- npos + 1
                    }
                    if(object@sep.index.ratio[i + (j - 1) * ns, 3] < 1){
                      nneg <- nneg + 1
                    }
                    x[j, ] <- object@sep.index.ratio[i + (j - 1) * ns, 3]
                  }
                }
                y[i, 1] <- round(min(x, na.rm = TRUE), digits=2)
                y[i, 2] <- round(quantile(x, probs = 0.25, na.rm = TRUE), digits=2)
                y[i, 3] <- round(median(x, na.rm = TRUE), digits=2)
                y[i, 4] <- round(mean(x, na.rm = TRUE), digits=2)
                y[i, 5] <- round(quantile(x, probs = 0.75, na.rm = TRUE), digits=2)
                y[i, 6] <- round(max(x, na.rm = TRUE), digits=2)
                y[i, 7] <- round(npos, digits=2)
                y[i, 8] <- round(nneg, digits=2)
              }
              print(y)



              cat("\n")
              cat("Temporal summary", "\n")
              cat("\n")
              x <- matrix(data = 0, nrow = ns, ncol = 1)
              y <- matrix(data = 0, nrow = nt, ncol = 8)
              colnames(y) <- c(" Min. ", "  Q1  ", "Median", " Mean ", "  Q3  ", " Max. ", "  #>1  ", "  #<1  ")
              for(i in 1:nt){
                ratio.sp <- 0
                npos <- 0
                nneg <- 0
                for(j in 1:ns){
                  if(!is.na(object@sep.index.ratio[j + (i - 1) * ns, 3]) == TRUE){
                    if(object@sep.index.ratio[j + (i - 1) * ns, 3] > 1){
                      npos <- npos + 1
                    }
                    if(object@sep.index.ratio[j + (i - 1) * ns, 3] < 1){
                      nneg <- nneg + 1
                    }
                    x[j, ] <- object@sep.index.ratio[j + (i - 1) * ns, 3]
                  }
                }
                y[i, 1] <- round(min(x), digits=2)
                y[i, 2] <- round(quantile(x, probs = 0.25), digits=2)
                y[i, 3] <- round(median(x), digits=2)
                y[i, 4] <- round(mean(x), digits=2)
                y[i, 5] <- round(quantile(x, probs = 0.75), digits=2)
                y[i, 6] <- round(max(x), digits=2)
                y[i, 7] <- round(npos, digits=2)
                y[i, 8] <- round(nneg, digits=2)
              }
              print(y)
            }
            else{
              print("All values are not available, hence the summary has not benn computed.")
            }


          })

