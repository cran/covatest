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
#'
#' @rdname blocks-class
#' @exportClass blocks
setClass("blocks", slots = c(mat.block = "matrix",
                                       array.block = "array"))

#' @param lb integer, length of each block. The number of terms in each block
#' must be greater than 5 and smaller than the quarter part of the length of
#' each time series
#' @param ls integer, number of overlapped data between two consecutive blocks.
#' The number of overlapped terms between two consecutive blocks must in the
#' interval [0,lb/2]
#' @param matdata STFDF/STSDF or \code{data frame}; which contains the
#' coordinates of the spatial points, the identification code of the spatial
#' points, the indentification code of the temporal points and the values of
#' the variable, typically output from \code{dataprep}
#' @param stpairs object of class \code{couples}, containing the spatial
#' points and the corresponding temporal lags to be analyzed
#'
#' @details
#' The function requires the user to set some external arguments. In particular,
#' if the spatio-temporal data are given as a \code{data} \code{frame} it is
#' necessary to specify
#' \itemize{
#' \item the column in which the spatial ID is stored
#' \item the column in which the values of the variable are stored.
#' }
#' On the other hand, if the data are given as a STFDF/STSDF it is necessary to
#' specify
#' \itemize{
#' \item the number of variables in the STFDF/STSDF
#' \item the slot in which the values of the variable of interest are stored
#' (only if more than one variable is stored in the STFDF/STSDF).
#' }
#'
#' Moreover, a message informs the user of the number of blocks computed. The
#' user can choose to continue or not with this number of blocks.
#' @note
#' \itemize{
#' \item "Error in matdata[, clvr]: subscript out of bounds" appears if the
#' second external arguments required for the given \code{data} \code{frame}
#' does not exist in the argument \code{matdata}
#'
#' \item If "Error in matdata[, clvr]" appears, no data for some of the
#' spatial points, specified in \code{stpairs}, are available. The user has to
#' go back to \code{couples} and revise the vector of the selected spatial points
#'
#' \item A stop running message occurs if the length of the time series for each
#' spatial points is less than 29
#'
#' \item A message appears if the length of the time series for each
#' spatial point is greater than 29 and less than 89, since the length of the
#' time series is low and may not guarantee the reliability of the tests
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
#' Cappello, C., De Iaco, S., Posa, D., 2016, Testing the type of
#' non-separability and some classes of covariance models for space-time data.
#' (submitted)
#'
#' @examples
#' # In the example regarding the simmetry test (typetest = 0), the length of
#' # each block is equal to 40 (lb=40) and the number of overlapped data between
#' # two consecutive blocks is equal to 10 (ls=10). In this way, 24 blocks, for
#' # each time series of 730 data have been obtained [24 = (730 - 1) / (40 - 10)].
#' # This number of blocks is consistent with respect the number of comparisons
#' # (24) defined in couples. Moreover, it is necessary to specify some external
#' # arguments: the number of variables in the analyzed STFDF (rr_13) and to
#' # confirm whether or not one wants to proceed
#' # with the number of blocks obtained with the aforementioned combination of
#' # lb and ls
#' #
#' # To run the example, paste and copy the following lines
#' # (without the symbol '#') in the console
#' #
#' # coupl_sim <- couples(typetest = 0, typecode = character())
#' # blocks_sim <- blocks(lb = 40, ls = 10, matdata = rr_13, stpairs = coupl_sim)
#' # 1
#' # Y
#'
#' @seealso \code{\link{couples}}
#' @seealso \code{\link{dataprep}}
#' @rdname blocks-class
#' @export
blocks <- function(lb, ls, matdata, stpairs) {
  # blocks creates a matrix which cointains all blocks of the same length
  # that can be extracted from each time series
  selstaz <- stpairs@sel.staz

  # ==================================#
  # == Some checks on class of data ==#
  # ==================================#

  if (class(matdata) == "matrix" || class(matdata) == "data.frame") {
    iclsp.id <- as.integer(readline(prompt = "Enter the column in which the spatial id is stored: "))
    iclvr <- as.integer(readline(prompt = "Enter the column in which the values of the variable are stored: "))
  }
  flag <- 0
  if (is.vector(selstaz) == TRUE && length(selstaz) >= 2) {
    matrix.names.matblock <- vector(mode = typeof(selstaz),
                                    length = length(selstaz))
    for (i in 1:length(selstaz)) {


      ### data in GSLIB format###
      if (class(matdata) == "matrix" || class(matdata) == "data.frame") {
        selstaz.names <- matdata[, iclsp.id]
        selstaz.inter <- intersect(selstaz.names, selstaz)
        if (length(selstaz.inter) != length(selstaz)) {
          stop("No data for some of the selected spatial points. Please go back to the function 'couples' and revise the vector of the selected spatial points")
        }

        if (is.numeric(matdata[, iclvr]) == FALSE)
          stop("Check the column in which the values of the variable are stored. Data must be numeric")

        datistaz <- matdata[matdata[, iclsp.id] == selstaz[i],
                            iclvr]
      } else {

        ### data in gstat format###
        if (class(matdata) == "STFDF") {
          if (i == 1) {
            selstaz.names <- row.names(matdata@sp)
            selstaz.inter <- intersect(selstaz.names, selstaz)
            if (length(selstaz.inter) != length(selstaz)) {
              stop("No data for some of the selected spatial points.Please go back to the function 'couples' and revise the vector of the selected spatial points")

            }
            nvr <- readline(prompt = "Enter the number of variables in the STFDF: ")
            nvr <- as.integer(nvr)
            iclvr <- 1
            if (nvr > 1) {
              iclvr <- as.integer(readline(prompt = "Enter the slot in which the values of the variable of interest are stored: "))
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
                stop("No data for some of the selected spatial points.Please go back to the function 'couples' and revise the vector of the selected spatial points")

              nvr <- readline(prompt = "Enter the number of variables in the STSDF ")
              nvr <- as.integer(nvr)
              iclvr <- 1
              if (nvr > 1) {
                iclvr <- as.integer(readline(prompt = "Enter the slot in which the values of the variable of interest are stored: "))
              }
            }
            datistaz <- matrix(matdata[selstaz[i], ], ncol = (1 + nvr))[, iclvr]

          } else {
            stop("The class of data must be matrix (gslib format), data.frame, STFDF or STSDF")
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
            stop("Check the first two parameters (lb,ls). Data must be integer")
          }

          if (ls < 0 || ls > (lb / 2)) {
            stop("The number of overlapped terms between two consecutive blocks must be in [0,lb/2]")
          }

          if (lb <= 5) {
            stop("The number of terms in each block must be greater than 5")
          }

          if (lt <= 29) {
            stop("The length of the time series (equal to ", lt,
                 ") for each spatial point must be greater than 29")
          }
          if (lt <= 89 && lt > 29) {

            message("*****************************************************************************")
            message("* The length of the time series (equal to ",lt, ") for each spatial point   *")
            message("* is low and may not guarantee the reliability of the some tests.           *")
            message("* See the manual for more details                                           *")
            message("*****************************************************************************")

            ans_YN <- readline(prompt = "Would you like to continue? (Y/N)")
            if (ans_YN == "N" || ans_YN == "n") {
              stop("Stop running")
            }
          }

          if (lb > (lt / 2)) {
            stop("The number of terms in each block must be less than a quarter part of the length of each time series")
          }
          nflag <- 0
        }
        #=======================================#
        #== compute the number of blocks (nb) ==#
        #=======================================#

        ld <- lb - ls

        nb <- as.integer(((length(datistaz)) - 1) / ld)

        if (nb <= 5) {
          message("Warning: the number of blocks is not greater than 5")
        }
        if (nb <= 0) {
          stop("The number of blocks is not positive")
        }


        message("*****************************************************************************")
        message("* The number of blocks computed, on the basis of the arguments, is ", nb, " *")
        message("* The number of blocks has to be consistent with the number of contrasts    *")
        message("* See the manual for more details                                           *")
        message("*****************************************************************************")

        ans_YN <- readline(prompt = "Would you like to continue with this number of blocks? (Y/N)")
        if (ans_YN == "N" || ans_YN == "n") {
          stop("Stop running")
        }
        arrayblock <- array(data = NA, dim = c(lb, nb, length(selstaz)))

      }


      bmat1 <- matrix(0, nrow = lb, ncol = nb)

      n1 <- 1
      n2 <- lb

      lused <- lb + (nb - 2) * ld


      # if(((nb+1)*ld)-lt>(0.15*lb)&&flag==0){
      if ((lb - (ls + (lt - lused))) > (0.15 * lb) && flag == 0) {
        flag <- 1

        warning("The last block is not complete ( missing values are greater than 15% of the fixed length) ")

      }
      for (j in 1:nb) {
        bmat1[, j] <- c(datistaz[n1:n2])
        n1 <- n1 + ld
        n2 <- n2 + ld

      }

      if (i == 1) {
        matblock <- bmat1
      }
      if (i <= length(selstaz) - 1 & i > 1) {
        matblock <- cbind(matblock, bmat1)

      }

      arrayblock[, , i] <- bmat1
      matrix.names.matblock[i] <- paste("BlockStaz", selstaz[i], sep = "_")

    }


    matblock <- cbind(matblock, bmat1)
    arrayblock[, , length(selstaz)] <- bmat1

    array.block <- array(arrayblock, dim = c(lb, nb, length(selstaz)),
                         dimnames = list(NULL, NULL, matrix.names.matblock))


    new("blocks", mat.block = matblock, array.block = array.block)


  } else {
    stop("The number of spatial points selected in function 'couples' must be a vector with at least two components")
  }


}
#' @include sepindex.R couples.R
NULL
