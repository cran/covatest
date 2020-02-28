#' Class "blocks"
#'
#' A class for overlapped blocks of the time series associated with the given
#' spatial points specified in \code{couples}. Thus, it is necessary to
#' execute \linkS4class{couples} first and then \linkS4class{blocks}
#'
#' @slot mat.block matrix of dimension (\emph{lb} x \emph{overall
#' number of blocks}); the columns of this matrix are associated with the
#' different blocks, of length equal to \code{lb}, that can be extracted
#' from the time series related to the selected spatial points
#' defined in \code{stpairs@sel.staz}
#' @slot array.block array of dimension (\emph{lb} x \emph{number of
#' blocks for each selected spatial points} x \emph{number of spatial points}).
#' In each table of this array, the overlapped blocks for each spatial location
#' are available
#' @slot sel.staz numeric or character; contains the ID codes of the selected
#' spatial points
#'
#' @rdname blocks-class
#' @exportClass blocks
setClass("blocks", slots = c(mat.block = "matrix",
                                       array.block = "array",
                                       sel.staz = "ANY"))

#' @param lb integer, length of each block. The number of terms in each block
#' must be greater than 5 and smaller than the quarter part of the length of
#' each time series
#' @param ls integer, number of overlapped data between two consecutive blocks.
#' The number of overlapped terms between two consecutive blocks must in the
#' interval [0, lb/2]
#' @param matdata STFDF/STSDF or \code{data.frame}; which contains the
#' coordinates of the spatial points, the identification code of the spatial
#' points, the indentification code of the temporal points and the values of
#' the variable, typically output from \code{read.STdata}
#' @param pardata1 integer, it represents the column in which the spatial ID is
#' stored (if the spatio-temporal data set is given as data.frame) or the number
#' of variables in the STFDF/STSDF (if the data are given as a STFDF/STSDF)
#' @param pardata2 integer, it represents the column in which the values of the
#' variable are stored (if the spatio-temporal data set is given as data.frame)
#' or the slot in which the values of the variable of interest are stored
#' (if the data are given as a STFDF/STSDF). Note that for STFDF/STSDF the
#' argument is set, by default, equal to 1 if the number of variables is equal
#' to 1
#' @param stpairs object of class \code{couples}, containing the spatial
#' points and the corresponding temporal lags to be analyzed
#'
#' @details
#' A message informs the user of the number of blocks extracted
#'
#'
#' @note
#' \itemize{
#' \item "Error in matdata[, clvr]: subscript out of bounds" appears if \code{pardata2}
#' does not exist in the argument \code{matdata}
#'
#' \item If "Error in matdata[, clvr]" appears, no data for some of the
#' spatial points, specified in \code{stpairs}, are available. The user has to
#' go back to \code{couples} and revise the vector of the selected spatial points
#' (\code{sel.staz} and \code{sp.couples.in} arguments)
#'
#' \item A stop occurs if more than 75\% of consecutive data are missing in the time
#' series, since a large number of missing values do not guarantee the reliability
#' of the tests
#'
#' \item A stop occurs if the length of the time series for each
#' spatial points is less than 29
#'
#' \item A message appears if the length of the time series for each
#' spatial point is greater than 29 and less than 89, since the length of the
#' time series is low and may not guarantee the reliability of the tests
#'
#' \item A stop occurs if more than 80\% of consecutive data are missing in one
#' of the blocks, since the estimation of the covariance matrix is not reliable,
#' when a large number of missing values occur
#'
#' \item If, in the last block of each selected spatial point, more than 15\%
#' of data are missing a warning message appears, since the estimation of the
#' covariance matrix, when a large number of missing values occurs, is not
#' reliable
#'
#' \item A warning message appears if the number of blocks, computed by fixing
#' \code{lb} and \code{ls}, is less than 5. It is convenient that the number of
#' blocks is close to the number of spatio-temporal comparisons defined in
#' \code{couples}. This avoids singolarity in computing test statistics
#'}
#'
#' @references
#' Li, B., Genton, M.G., Sherman, M., 2007, A nonparametric assessment
#' of properties of spacetime covariance functions.
#' Journal of the American Statistical Association, \bold{102} 736--744.
#'
#' De Iaco, S., Palma, M., Posa, D., 2016. A general procedure for selecting a
#' class of fully symmetric space-time covariance functions.
#' Environmentrics, \bold{27(4)} 212--224.
#'
#' Cappello, C., De Iaco, S., Posa, D., 2018, Testing the type of
#' non-separability and some classes of covariance models for space-time data.
#' Stochastic Environmental Research and Risk Assessment,
#' \bold{32} 17--35
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
#'
#' sel.staz.sym <- c("DERP016", "DENW065", "DEHE051", "DETH026", "DENW063", "DENI019",
#' "DENW068", "DEHE046", "DEUB029", "DEBY047", "DETH061", "DESN049")
#'
#' sp.couples.in.sym <- matrix(data = c("DERP016", "DENW065", "DEHE051", "DETH026",
#' "DENW063", "DENI019", "DENW068", "DEHE046", "DEUB029", "DEBY047", "DETH061", "DESN049"),
#' ncol = 2, byrow = TRUE)
#'
#' t.couples.in.sym <- c(1, 2)
#'
#' couples.sym <- couples(sel.staz = sel.staz.sym, sp.couples.in = sp.couples.in.sym,
#' t.couples.in = t.couples.in.sym, typetest = "sym", typecode = character())
#'
#' block.sym <- blocks(lb = 40, ls = 10, matdata = rr_13, pardata1 = 1, pardata2 = 1,
#' stpairs = couples.sym)
#'
#' ### methods for blocks
#' #1. show
#' block.sym
#'
#' #2. [  extract
#' block.sym[1,] #select the 1st block of each spatial location
#' block.sym[,1] #select all blocks of the 1st spatial location
#' block.sym[1:2, 1:3] #select the first two blocks of the first 3 spatial locations
#'
#' #3. summary
#' summary(block.sym, 1:2, 1:3) #to obtain the summary associated to the first
#' #two blocks of the first 3 spatial locations
#'
#' summary(block.sym, 0, 1) #to obtain the summary associated to all blocks of
#' #the 1st spatial location
#'
#' #4. boxplot
#' boxplot(block.sym, 1:5, 1:2) #boxplots of the first 5 blocks of associated to
#' #the first 2 spatial locations
#'
#' boxplot(block.sym, 0 ,1) #boxplots of all blocks of associated to the 1st
#' #spatial location

#'
#' @seealso \code{\link{couples}}
#' @seealso \code{\link{read.STdata}}
#' @rdname blocks-class
#' @export
blocks <- function(lb, ls, matdata, pardata1, pardata2, stpairs) {
  # blocks creates a matrix which cointains all blocks of the same length
  # that can be extracted from each time series

  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}

  ### SOME CHECKS ON THE ARGUMENTS ###

  if (is.scalar(lb) == FALSE || is.scalar(ls) == FALSE  || is.scalar(pardata1) == FALSE ||
      is.scalar(pardata2) == FALSE) {
    message("Start error message. Some of the arguments are not numeric.")
    stop("End error message. Stop running.")
  }


  if(lb != as.integer(lb) || ls != as.integer(ls) || pardata1 != as.integer(pardata1) ||
     pardata2 != as.integer(pardata2)){
    lb <- as.integer(lb)
    ls <- as.integer(ls)
    pardata1 <- as.integer(pardata1)
    pardata2 <- as.integer(pardata2)
    message("Warning message: the arguments expected to be integer are forced to be integer numbers.")
  }


  if (!inherits(stpairs, "couples")){
    message("Start error message. stpairs argument has to be of class couples.")
    stop("End error message. Stop running.")
  }


  selstaz <- stpairs@sel.staz
  # ==================================#
  # == Some checks on class of data ==#
  # ==================================#

  if (class(matdata) == "matrix" || class(matdata) == "data.frame") {
    iclsp.id <- as.integer(pardata1)
    iclvr <- as.integer(pardata2)
  }
  flag <- 0
  info.na <- NA
  info.nna29 <- NA
  info.nna89 <- NA
  if (is.vector(selstaz) == TRUE && length(selstaz) >= 2) {
    matrix.names.matblock <- vector(mode = typeof(selstaz),
                                    length = length(selstaz))
    for (i in 1:length(selstaz)) {
      ### data in GSLIB format###
      if (class(matdata) == "matrix" || class(matdata) == "data.frame") {
        if (i == 1) {
          selstaz.names <- matdata[, iclsp.id]
        selstaz.inter <- intersect(selstaz.names, selstaz)
        if (length(selstaz.inter) != length(selstaz)) {
          message("Start error message. No data for some of the selected spatial points. Please go back to the function 'couples' and revise the vector of the selected spatial points.")
          stop("End error message. Stop running.")
        }
        }

        if (is.numeric(matdata[, iclvr]) == FALSE){
          message("Start error message. Check the column in which the values of the variable are stored. Data must be numeric.")
          stop("End error message. Stop running.")}

        datistaz <- matdata[matdata[, iclsp.id] == selstaz[i],
                            iclvr]
      } else {

        ### data in gstat format###
        if (class(matdata) == "STFDF") {
          if (i == 1) {
            selstaz.names <- row.names(matdata@sp)
            selstaz.inter <- intersect(selstaz.names, selstaz)
            if (length(selstaz.inter) != length(selstaz)) {
              message("Start error message. No data for some of the selected spatial points.Please go back to the function 'couples' and revise the vector of the selected spatial points.")
              stop("End error message. Stop running.")

            }
            nvr <- as.integer(pardata1)
            iclvr <- as.integer(pardata2)
            if (nvr == 1) {
              iclvr <- 1
            }
          }

          datistaz <- matrix(matdata[selstaz[i], ], ncol = (1 + nvr))[, iclvr]

        } else {


          ### data in gstat format###
          if (class(matdata) == "STSDF") {
            matdata <- as(matdata, "STFDF")
            if (i == 1) {
              selstaz.names <- row.names(matdata@sp)
              selstaz.inter <- intersect(selstaz.names, selstaz)
              if (length(selstaz.inter) != length(selstaz))
                message("Start error message. No data for some of the selected spatial points.Please go back to the function 'couples' and revise the vector of the selected spatial points.")
                stop("End error message. Stop running.")

              nvr <- as.integer(pardata1)
              iclvr <- as.integer(pardata2)

              if (nvr == 1) {
                iclvr <- 1
              }
            }
            datistaz <- matrix(matdata[selstaz[i], ], ncol = (1 + nvr))[, iclvr]

          } else {
            message("Start error message. The class of data must be matrix (gslib format), data.frame, STFDF or STSDF.")
            stop("End error message. Stop running.")
          }
        }
      }

      if (i == 1) {
        lt <- length(datistaz)
        # ======================================#
        # ===== Some checks on lb, ls, lt  =====#
        # ======================================#
        nflag <- 1
        while (nflag == 1) {

          is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) abs(x -
                                                                round(x)) < tol

          if (is.wholenumber(lb) == FALSE || is.wholenumber(ls) == FALSE) {
            message("Start error message. Check the first two parameters (lb,ls). Data must be integer.")
            stop("End error message. Stop running.")
          }

          if (ls < 0 || ls > (lb / 2)) {
            message("Start error message. The number of overlapped terms between two consecutive blocks must be in [0,lb/2].")
            stop("End error message. Stop running.")
          }

          if (lb <= 5) {
            message("Start error message. The number of terms in each block must be greater than 5.")
            stop("End error message. Stop running.")
          }

          if (lt <= 29) {
            message("Start error message. The length of the time series (equal to ", lt, ") for each spatial point must be greater than 29.")
            stop("End error message. Stop running.")
          }
          if (lt <= 89 && lt > 29) {

            message("*****************************************************************************")
            message("* The length of the time series (equal to ",lt, ") for each spatial point   *")
            message("* is low and may not guarantee the reliability of the tests.                *")
            message("* See the manual for more details.                                          *")
            message("*****************************************************************************")

           }


          if (lb > (lt / 4)) {
            message("Start error message. The number of terms in each block must be less than a quarter part of the length of each time series.")
            stop("End error message. Stop running.")
          }
          nflag <- 0
        }
        #=======================================#
        #== compute the number of blocks (nb) ==#
        #=======================================#

        ld <- lb - ls

        nb <- as.integer(((length(datistaz)) - 1) / ld)

        if (nb <= 5) {
          message("Warning message: the number of blocks is not greater than 5.")
        }
        if (nb <= 0) {
          message("Start error message. The number of blocks is not positive.")
          stop("End error message. Stop running.")
        }


        message("*****************************************************************************")
        message("* The number of blocks computed, on the basis of the arguments, is ", nb, " *")
        message("* The number of blocks has to be consistent with the number of contrasts    *")
        message("* See the manual for more details.                                          *")
        message("*****************************************************************************")

        arrayblock <- array(data = NA, dim = c(lb, nb, length(selstaz)))
        matblock <- matrix(0, nrow = lb, ncol = nb)
      }


      #====================================================#
      #== compute the number of consecutive missing values #
      #====================================================#


      count.na <- matrix(0, nrow = lt, ncol = 1)
      count.cons.na <- 0
      count.nna <- matrix(0, nrow = lt, ncol = 1)
      count.cons.nna <- 0
      for (ii in 1:lt) {
        if(is.na(datistaz[ii]) == TRUE){
          count.cons.na <- count.cons.na + 1
          count.na[ii, 1] <- count.cons.na
        }else{count.cons.na<- 0}
        if(is.na(datistaz[ii]) == FALSE){
          count.cons.nna <- count.cons.nna + 1
          count.nna[ii, 1] <- count.cons.nna
        }else{count.cons.nna<- 0}
      }
      max.count.na <-   max(count.na[,])/lt
      max.count.nna <-   max(count.nna[,])
      if(max.count.nna > 29){
      if(max.count.nna <= 89){
        if(is.na(info.nna89[1]) == TRUE){
          info.nna89 <- selstaz[i]
        }else{
          info.nna89 <- rbind(info.nna89, selstaz[i])
        }
              }
      if(max.count.na < 0.75){
      bmat1 <- matrix(0, nrow = lb, ncol = nb)
      n1 <- 1
      n2 <- lb
      lused <- lb + (nb - 2) * ld
      if ((lb - (ls + (lt - lused))) > (0.15 * lb) && flag == 0) {
        flag <- 1
        message("Warning message: the last block is not complete (missing values are greater than 15% of the fixed length) ")
      }
      for (j in 1:nb) {
        bmat1[, j] <- c(datistaz[n1:n2])
        n1 <- n1 + ld
        n2 <- n2 + ld
      }
        if (i == 1){
        matblock <- bmat1
      }else{
      if (i <= length(selstaz) - 1 & i > 1) {
        matblock <- cbind(matblock, bmat1)
      }
      }
      arrayblock[, , i] <- bmat1
      matrix.names.matblock[i] <- paste("SpatialPoint", selstaz[i], sep = "_")

      }else{
        if(is.na(info.na[1]) == TRUE){
          info.na <- selstaz[i]
        }else{
          info.na <- rbind(info.na, selstaz[i])
        }
      }

      }else{
        if(is.na(info.nna29[1]) == TRUE){
          info.nna29 <- selstaz[i]
        }else{
          info.nna29 <- rbind(info.nna29, selstaz[i])
        }
    }

      #==========================================================#
      #== END computing the number of consecutive missing values #
      #==========================================================#

    }
    #==========================================================#
    #== END cicle over selstaz #                               #
    #==========================================================#
    if(is.na(info.na[1]) == FALSE){
      message("Start error message. Too many consecutive NAs (greater than 75%). The following spatial points are non-admissible:")
      for (i in 1:length(info.na)){
        message((info.na[i]))
      }
      message("Please exclude/change the non-admissible spatial points from the selection.")
      stop("End error message. Stop running.")
    }

    if(is.na(info.nna29[1]) == FALSE){
      message("Start error message. The following spatial points are non-admissible. The number of valid consecutive values must be greater than 29.")
      for (i in 1:length(info.nna29)){
        message((info.nna29[i]))
      }
      message("Please exclude/change the non-admissible spatial points from the selection.")
      stop("End error message. Stop running.")
    }
    if(is.na(info.nna89[1]) == FALSE){
      message("Warning message: the number of valid consecutive values of the following spatial points is low (<=89) and may not guarantee the reliability of the tests.")
      for (i in 1:length(info.nna89)){
        message((info.nna89[i]))
      }
    }

    matblock <- cbind(matblock, bmat1)
    arrayblock[, , length(selstaz)] <- bmat1

    array.block <- array(arrayblock, dim = c(lb, nb, length(selstaz)),
                         dimnames = list(NULL, NULL, matrix.names.matblock))


    #==============================================================#
    #== compute the number of consecutive missing values per block #
    #==============================================================#

    info.block.na <- NA
    count.na<- matrix(0, nrow = lb, ncol = 1)
    count.cons.na <- 0
    for (ii in 1:length(selstaz)) {
      kk <- 0
      while (kk < nb) {
        kk <- kk + 1
        for (jj in 1:lb) {
          if(is.na(array.block[jj, kk, ii]) == TRUE){
            count.cons.na <- count.cons.na + 1
            count.na[jj, 1] <- count.cons.na
          }else{count.cons.na <- 0}
        }
        max.count.na <- max(count.na[, 1])/lb
        if(max.count.na >= 0.80){
          kk <- nb+1
          if(is.na(info.block.na[1]) == TRUE){
            info.block.na <- selstaz[ii]
          } else{
            info.block.na <- rbind(info.block.na, selstaz[ii])
          }
          #message(selstaz[i],'is a non-admissible spatial point. Too many consecutive NAs in the blocks.')
        }
      }
    }


    if(is.na(info.block.na[1]) == FALSE){
      message("Start error message. Too many consecutive NAs per block (greater than 80%). The following spatial points are non-admissible:")
      for (ii in 1:length(info.block.na)) {
         message(info.block.na[ii])
      }
      message("Please exclude/change the non-admissible spatial points from the selection.")
      stop("End error message. Stop running.")
    }

    #====================================================================#
    #== END computing the number of consecutive missing values per block #
    #====================================================================#


    new("blocks", mat.block = matblock, array.block = array.block, sel.staz = selstaz)


  } else {
    message("Start error message. The number of spatial points selected in function 'couples' must be a vector with at least two components.")
    stop("End error message. Stop running.")
  }


}
#' @include sepindex.R couples.R
NULL
#' @param x object of class \code{blocks} for methods \code{boxplot} and \code{extract}
#' @param i index specifing the block to be selected. If \code{i=0} all blocks are
#' selected automatically (option available only for \code{boxplot} and \code{summary}
#' methods)
#' @param j index specifing the spatial point to be selected. If \code{j=0} all
#' spatial points are selected automatically (option available only for \code{boxplot}
#' and \code{summary} methods)
#' @param ... any arguments that will be passed to the panel plotting functions
#' @rdname blocks-class
#' @aliases blocks-method boxplot
#' @export
setMethod("boxplot", signature = c(x = "blocks"),
          function(x, i, j, ...) {
            is.scalar <- function (y){length(y) == 1L && is.vector(y, mode = "numeric")}
            if(is.scalar(i) == TRUE){
              if(i == 0){
              i <- c(1:dim(x@array.block)[2])
              }
            }

            if(is.scalar(j) == TRUE){
              if(j == 0){
              j = c(1:dim(x@array.block)[3])
              }
              }

            if(is.scalar(j) == FALSE && is.scalar(i) == FALSE){
              z <- x@array.block[, i, j]
              nstat <- dim(z)[3]

              nb <- dim(z)[2]

              lb <- dim(z)[1]
              w <- as.vector(dimnames(z)[3][[1]])


            }
            if(is.scalar(j) == TRUE || is.scalar(i) == TRUE){
              z <- as.matrix(x@array.block[, i, j])
              if(is.scalar(j) == TRUE){
                nstat <- 1
                nb <- ncol(z)
                lb <- nrow(z)
                w <- dimnames(x@array.block)[3][[1]]}
              if(is.scalar(i) == TRUE){
                nstat <- ncol(z)
                nb <- 1
                lb <- nrow(z)
                w <- dimnames(x@array.block)[3][[1]]}
            }

            data.block <- matrix(NA, ncol = 2, nrow = lb*nb)
            for(k in 1:nstat){
              if(is.scalar(j) == FALSE && is.scalar(i) == FALSE){
                z.data <- as.matrix(z[,,k])
              }else{z.data <- z}
              for(l in 1:nb){

                data.block[(1+(l-1)*lb):(lb+(l-1)*lb), 1] <- z.data[,l]
                data.block[(1+(l-1)*lb):(lb+(l-1)*lb), 2] <- rep(l,lb)
              }
              cat("\n")
              cat("Box plot for ", w[k], "built \n")
              boxplot(data.block[,1] ~ data.block[,2],
                      data = data.block, xlab="Block", ylab = "Value", main = w[k], ...)
            }

          }

)
#' @param object object of class \code{blocks} for methods \code{show} and \code{summary}
#' @rdname blocks-class
#' @aliases blocks-method
#' @export
setMethod(f="show", signature="blocks", definition=function(object) {
  lb <- nrow(object@mat.block)
  nsp <- dim(object@array.block)[3]
  nb <- dim(object@array.block)[2]
  NA_data<- round((length(object@array.block[, nb, 1][is.na(object@array.block[, nb, 1])])/lb)*100, digits=3)

  cat("An object of class 'blocks', with", "\n")
  cat("number of spatial points = ", nsp, "\n")
  cat("number of blocks = ", nb, "\n")
  cat("length of each block = ", nrow(object@mat.block), "\n")
  cat("% of NA in the last block = ", NA_data, "\n")
  cat("\n")
  cat("Slot 'mat.block':")
  cat("\n")
  print(object@mat.block)
  cat("\n")
  cat("Slot 'array.block':")
  cat("\n")
  print(object@array.block)
  cat("\n")
  cat("Slot 'sel.staz':")
  cat("\n")
  print(object@sel.staz)
}
)
#' @rdname blocks-class
#' @aliases blocks-method
#' @export
setMethod(f="[", signature="blocks", definition=function(x, i, j) {

  x@array.block[, i, j]

}
)
#' @rdname blocks-class
#' @aliases blocks-method
#' @export
setMethod(f = "summary", signature = "blocks",
          definition = function(object, i, j) {

            is.scalar <- function (y){length(y) == 1L && is.vector(y, mode = "numeric")}
            if(is.scalar(i) == TRUE){
              if(i == 0){
                i <- c(1:dim(object@array.block)[2])
              }
            }

            if(is.scalar(j) == TRUE){
              if(j == 0){
                j = c(1:dim(object@array.block)[3])
              }
            }

            if(is.scalar(j) == FALSE && is.scalar(i) == FALSE){
              z <- object@array.block[, i, j]

              nstat <- dim(z)[3]

              nb <- dim(z)[2]

              lb <- dim(z)[1]
              w <- as.vector(dimnames(z)[3][[1]])
            }
            if(is.scalar(j) == TRUE || is.scalar(i) == TRUE){
              z <- as.matrix(object@array.block[, i, j])
              if(is.scalar(j) == TRUE){
                nstat <- 1
                nb <- ncol(z)
                lb <- nrow(z)
                w <- dimnames(object@array.block)[3][[1]]}
              if(is.scalar(i) == TRUE){
                nstat <- ncol(z)
                nb <- 1
                lb <- nrow(z)
                w <- dimnames(object@array.block)[3][[1]]}
            }

            for(k in 1:nstat){
            cat("\n")
            cat("\n")
            cat("Descriptive statistics for ", w[k], "\n")
            cat("\n")
            y <- matrix(data = 0, nrow = nb, ncol = 7)
            colnames(y) <- c(" Min. ", "  Q1  ", "Median", " Mean ", "  Q3  ", " Max. ", "St. Dev.")
            for(l in 1:nb){
              x <- object@array.block[,l,k]
              y[l, 1] <- round(min(x, na.rm = TRUE), digits=2)
              y[l, 2] <- round(quantile(x, probs = 0.25, na.rm = TRUE), digits=2)
              y[l, 3] <- round(median(x, na.rm = TRUE), digits=2)
              y[l, 4] <- round(mean(x, na.rm = TRUE), digits=2)
              y[l, 5] <- round(quantile(x, probs = 0.75, na.rm = TRUE), digits=2)
              y[l, 6] <- round(max(x, na.rm = TRUE), digits=2)
              y[l, 7] <- round(sd(x, na.rm = TRUE), digits=2)
            }

            print(y)
            }


          })
