#' Class "covastat"
#'
#' A class for the sample spatio-temporal covariances for the specified spatial
#' and temporal lags, given in \code{stpairs} (object of class \code{couples}),
#' for symmetry, separability and type of non separability tests.
#' Depending on the type of test, the empirical variance, the sample spatial
#' and temporal marginal covariances are also computed
#'
#' @slot G matrix; containing the spatio-temporal covariances for the specified
#' lags. For all tests, except for the symmetry test (\code{typetest = "sym"}), the
#' sample variance and the sample spatial and temporal marginal covariances are
#' also computed and stored in \code{G}
#' @slot cova.h matrix; containing the sample spatial marginal covariances
#' for the specified lags
#' @slot cova.u matrix; containing the sample temporal marginal covariances
#' for the specified lags
#' @slot f.G array; containing the computation of specific functions of the
#' elements of \code{G}, see references
#' @slot B matrix; containing the computation of the derivatives of each element
#' of \code{f.G} with respect to each element of \code{G}
#' @slot A contrast matrix
#' @slot typetest character; contains the code of the test to be performed
#'
#'
#' @rdname covastat-class
#' @exportClass covastat
setClass("covastat", slots = c(G = "matrix",
                               cova.h = "matrix",
                               cova.u = "matrix",
                               f.G = "array",
                               B = "array",
                               A = "matrix",
                               typetest = "character"))

#' @param matdata STFDF/STSDF or \code{data.frame}; which contains the
#' coordinates of the spatial points, the identification code of the spatial
#' points, the indentification code of the temporal points and the values of
#' the variable, typically output from \code{read.STdata}
#' @param pardata1 integer, it represents the column in which the spatial ID is
#' stored (if the spatio-temporal data set is given as data.frame) or the number of
#' variables in the STFDF/STSDF (if the data are given as a STFDF/STSDF)
#' @param pardata2 integer, it represents the column in which the values of the
#' variable are stored (if the spatio-temporal data set is given as data.frame) or
#' the slot in which the values of the variable of interest are stored
#' (if the data are given as a STFDF/STSDF). Note that for STFDF/STSDF the
#' argument is set, by default, equal to 1 if the number of variables is equal to 1
#' @param stpairs object of class \code{couples}, containing the spatial
#' points and the corresponding temporal lags to be analyzed
#' @param typetest character; set \code{typetest = "sym"} for symmetry test
#' (default choice), \code{typetest = "sep"} for separability test, \code{typetest = "tnSep"}
#' for type of non separability test
#'
#'
#' @details
#' A message appears on the user's console if the \code{G} vector
#' contains spatio-temporal negative covariances. The message returns the negative
#' value/values and it will help to identify the spatial and the temporal lags
#' involved
#'
#'
#' @note
#' \itemize{
#' \item A stop occurs if the number of spatial points fixed in \code{stpairs}
#' (object of class \code{couples}) is less than 2
#' \item If \code{typetest = "sym"} (symmetry test) \code{cova.h}, \code{cova.u},
#' \code{f.G} and \code{B} are not available
#' \item A stop occurs if more than 75\% of consecutive data are missing in the time
#' series, since a large number of missing values do not guarantee the reliability
#' of the tests
#' }
#'
#' @seealso \linkS4class{couples}
#' @seealso \code{\link{read.STdata}}
#'
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
#' covast.sym <- covastat(matdata = rr_13, pardata1 = 1, pardata2 = 1,
#' stpairs = couples.sym, typetest = "sym")
#'
#' ### method for covastat
#' #1. show
#' covast.sym
#'
#' @rdname covastat-class
#' @export
covastat <- function(matdata, pardata1, pardata2, stpairs, typetest = "sym") {

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {abs(x - round(x)) <
      tol}

  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}

  ### SOME CHECKS ON THE ARGUMENTS ###

  if (is.scalar(pardata1) == FALSE || is.scalar(pardata2) == FALSE) {
    message("Start error message. Some of the arguments are not numeric.")
    stop("End error message. Stop running.")
  }


  if(pardata1 != as.integer(pardata1) || pardata2 != as.integer(pardata2)){

    pardata1 <- as.integer(pardata1)
    pardata2 <- as.integer(pardata2)
    message("Warning message: the arguments expected to be integer are forced to be integer numbers.")
  }


  if (!inherits(stpairs, "couples")){
    message("Start error message. stpairs argument has to be of class couples.")
    stop("End error message. Stop running.")
  }

  if(stpairs@typetest != typetest){
    message("Warning message: the argument typetest is different from the one defined in stpairs.")
  }


  selstaz <- stpairs@sel.staz
  couples <- stpairs@couples.st
  nstaz <- length(selstaz)

  #========================================================================#
  #= type of tests: 0=symmetry (default choice), 1=separability,          =#
  #= 2= type of non separability, 3=type of variability                   =#
  #========================================================================#

  if (is.character(typetest) == FALSE) {
    message("Start error message. The argument for typetest is not admissible.")
    stop("End error message. Stop running.")
  }

  if (typetest != "sym" && typetest != "sep"  && typetest != "tnSep") {
    message("Start error message. The argument for typetest is not admissible.")
    stop("End error message. Stop running.")
  }

  if (typetest == "sym") {
    type.test <- 0
  }else{if (typetest == "sep"){
    type.test <- 1
  }else{
    type.test <- 2 # tnSep
  }
  }



  #========================================================================#
  #= Select data for the stations of interest                             =#
  #========================================================================#

  if (class(matdata) == "matrix" || class(matdata) == "data.frame") {
    iclsp.id <- as.integer(pardata1)
    iclvr <- as.integer(pardata2)
  }
  info.na <- NA
  info.nna29 <- NA
  info.nna89 <- NA
  if (is.vector(selstaz) && length(selstaz) >= 2) {
    for (i in 1:length(selstaz)) {

      #= data in GSLIB format =#
      if (class(matdata) == "matrix" || class(matdata) == "data.frame") {

        if (is.numeric(matdata[, iclvr]) == TRUE) {

          if (i == 1) {
            selstaz.names <- matdata[, iclsp.id]
            selstaz.inter <- intersect(selstaz.names, selstaz)
            if (length(selstaz.inter) != length(selstaz)) {
              message("Start error message. No data for some of the selected spatial points. Please go back to the function 'couples' and revise the vector of the selected spatial points.")
              stop("End error message. Stop running.")

            }
            }

          datistaz <- matdata[matdata[, iclsp.id] == selstaz[i], iclvr]
        } else {
          message("Start error message. Check the column in which the values of the variable are stored. Data must be numeric.")
          stop("End error message. Stop running.")
        }

        } else {

          #= data in gstat format =#
          if (class(matdata) == "STFDF") {
            if (i == 1) {
              selstaz.names <- row.names(matdata@sp)
              selstaz.inter <- intersect(selstaz.names, selstaz)
              if (length(selstaz.inter) != length(selstaz)) {
                message("Start error message. No data for some of the selected spatial points. Please go back to the function 'couples' and revise the vector of the selected spatial points.")
                stop("End error message. Stop running.")

              }
              nvr <- as.integer(pardata1)
              iclvr <- as.integer(pardata2)
              if (nvr == 1) {
                iclvr <- 1
              }
            }
            datistaz <- matrix(matdata[selstaz[i], ], ncol = (1 + nvr))[,
                                                                        iclvr]
          } else {
            #= data in gstat format =#
            if (class(matdata) == "STSDF") {
              matdata <- as(matdata, "STFDF")
              if (i == 1) {
                  selstaz.names <- row.names(matdata@sp)
                  selstaz.inter <- intersect(selstaz.names, selstaz)
                  if (length(selstaz.inter) != length(selstaz)) {
                    message("Start error message. No data for some of the selected spatial points. Please go back to the function 'couples' and revise the vector of the selected spatial points.")
                    stop("End error message. Stop running.")

                  }
                  nvr <- as.integer(pardata1)
                  iclvr <- as.integer(pardata2)
                  if (nvr == 1) {
                    iclvr <- 1
                  }

              }
              datistaz <- matrix(matdata[selstaz[i], ], ncol = (1 + nvr))[,
                                                                          iclvr]

            } else {
              message("Start error message. The class of data must be matrix (gslib format), data.frame, STFDF or STSDF.")
              stop("End error message. Stop running.")
            }
            }
          }
      if (i == 1) {
        lt <- length(datistaz)
        # ======================================#
        # ===== Some checks on lt          =====#
        # ======================================#
      if (lt <= 29) {
        message("Start error message. The length of the time series (equal to ", lt,") for each spatial point must be greater than 29.")
        stop("End error message. Stop running.")
      }
      if (lt <= 89 && lt > 29) {
        message("*****************************************************************************")
        message("* The length of the time series (equal to ",lt, ") for each spatial point   *")
        message("* is low and may not guarantee the reliability of the some tests.           *")
        message("* See the manual for more details.                                          *")
        message("*****************************************************************************")
        }}

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

      if(max.count.na > 0.75){
      if(is.na(info.na[1]) == TRUE){
        info.na <- selstaz[i]
      }else{info.na <- rbind(info.na, selstaz[i]) }
      }

      if (i == 1) {
        matdata.sel <- datistaz
      }
      if (i > 1) {
        matdata.sel <- cbind(matdata.sel, datistaz)
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
      message("Start error message. The following spatial points are non-admissible. Too many consecutive NAs (greater than 75%).")
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
      message("Warning message: the number of valid consecutive values of the following spatial points is low (<=89) and may not guarantee the reliability of the some tests.")
      for (i in 1:length(info.nna89)){
        message((info.nna89[i]))
      }
    }
    }else {
      message("Start error message. The number of spatial points selected in function 'couples' must be a vector with at least two components.")
    stop("End error message. Stop running.")
    }

  #==========================#
  #=  End selection of data =#
  #==========================#

  #====================================#
  #= Search for lags equal to zero    =#
  #====================================#

  couples.nrow <- nrow(couples)

  nct <- 0

  for (i in 1:nrow(couples)) {
    for (j in 3:ncol(couples)) {
      if (couples[i, j] != 0) {
        nct <- nct + 1
      }
    }
  }



  array.matdata.sel <- array(matdata.sel, dim = c(length(datistaz), 1, length(selstaz)))
  cova.nv <- matrix(data = "-", nrow = nrow(couples), ncol = ncol(couples))
  #========================================================================#
  #= End of search for lags equal to zero                                 =#
  #========================================================================#

  #===========================================================================#
  #= start if on type of test (O=symmetry, 1=separability,2=type of non_sep, =#
  #= 3=variability)                                                          =#
  #===========================================================================#

  #========================================================================#
  #= start test on type.test=0, 1,2,3                                      =#
  #========================================================================#
   vec.na <- matrix(NA, length(datistaz), 1)
   info.na.all <- matrix(NA, 1, 5)
    nflag.cova.nv <- 0
  if (type.test == 0 || type.test == 1 || type.test == 2 || type.test == 3) {

    #= Compute spatio-temporal covariance (with hs and ht different from zero)
    #= for each block =#

    if (type.test == 0 || type.test == 1 || type.test == 2 || type.test == 3) {
      cova <- matrix(0, nrow = nct, ncol = 1)

      couples.ncol <- as.integer((ncol(couples) - 2))

      couples.ncol.r <- 0

      for (i in 1:couples.nrow) {
        couples.nrow.r <- 0
        nf <- 0
        for (j in 1:couples.ncol) {


          if (couples[i, j + 2] != 0) {

            couples.nrow.r <- couples.nrow.r + 1
            cov.n <- as.integer(couples.nrow.r + couples.ncol.r)
            if (couples[i, j + 2] > 0) {
              cova[cov.n, ] <- cov(array.matdata.sel[-(nrow(matdata.sel) -
                                                         couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                             1]], array.matdata.sel[-(1:couples[i, j + 2]), , couples[i,
                                                                                                                                                                      2]], use = "pairwise.complete.obs")
              if(is.na(cova[cov.n, ]) == TRUE){
                  if(is.na(info.na.all[1,1]) == TRUE){
                    info.na.all[1,1] <- couples[i,1]
                    info.na.all[1,2] <- couples[i,2]
                    info.na.all[1,5] <- couples[i,j +2]
                    if(identical(vec.na,array.matdata.sel[-(nrow(matdata.sel) -
                                                            couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                                1]]) == TRUE){info.na.all[1,3] <- 1}
                    if(identical(vec.na,array.matdata.sel[-(1:couples[i, j + 2]), , couples[i,
                                                                                            2]]) == TRUE){info.na.all[1,4] <- 2}
                  }else{
                    info.na <- rbind(info.na, c(couples[i,1:2],NA,NA,NA))
                    info.na[1,5] <- couples[i,j +2]
                    if(identical(vec.na,array.matdata.sel[-(nrow(matdata.sel) -
                                                            couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                                1]]) == FALSE){info.na.all[1,3] <- 1}
                    if(identical(vec.na,array.matdata.sel[-(1:couples[i, j + 2]), , couples[i,
                                                                                            2]]) == FALSE){info.na.all[1,4] <- 2}


                  }

              }
              }

            if (couples[i, j + 2] < 0) {
              cova[cov.n, ] <- cov(array.matdata.sel[-(1:(-couples[i,
                                                                   j + 2])), , couples[i, 1]], array.matdata.sel[-(nrow(matdata.sel) +
                                                                                                                     couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                                                                                         2]], use = "pairwise.complete.obs")
              if(is.na(cova[cov.n, ]) == TRUE){
                if(is.na(info.na.all[1,1]) == TRUE){
                  info.na.all[1,1] <- couples[i,1]
                  info.na.all[1,2] <- couples[i,2]
                  info.na.all[1,5] <- couples[i,j +2]
                  if(identical(vec.na,array.matdata.sel[-(1:(-couples[i,
                                                                      j + 2])), , couples[i, 1]]) == TRUE){info.na.all[1,3] <- 1}
                  if(identical(vec.na,array.matdata.sel[-(nrow(matdata.sel) +
                                                          couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                              2]]) == TRUE){info.na.all[1,4] <- 2}
                }else{
                  info.na <- rbind(info.na, c(couples[i,1:2],NA,NA,NA))
                  info.na[1,5] <- couples[i,j +2]
                  if(identical(vec.na,array.matdata.sel[-(1:(-couples[i,
                                                                      j + 2])), , couples[i, 1]]) == FALSE){info.na.all[1,3] <- 1}
                  if(identical(vec.na,array.matdata.sel[-(nrow(matdata.sel) +
                                                          couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                              2]]) == FALSE){info.na.all[1,4] <- 2}


                }

              }
              }

            if (cova[cov.n,] < 0) {
              nflag.cova.nv <- nflag.cova.nv + 1
              cova.nv[i,1:2] <- selstaz[couples[i,1:2]]
              cova.nv[i,j+2] <- couples[i,j+2]
            }

          }

        }
        if (nf == 0) {
          couples.ncol.r <- couples.ncol.r + couples.nrow.r
        }
        nf <- 1
      }

      if(is.na(info.na.all[1,1]) == FALSE){
        message("Start error message. There are no enough data for computing the covariance: spatial couples, #point in the couple non-valid, #point in the couple non-valid, temporal lag non-valid.")
        for (i in 1:length(info.na.all)){
          print(info.na.all[i,])
        }
        message("Please exclude/change the non-valid spatial couples/points/temporal lag from the selection.")
        stop("End error message. Stop running.")
      }

    if (nflag.cova.nv != 0) {
      message("Warning message: ", nflag.cova.nv, " negative spatio-temporal covariance/es are detected.")
      message("In the following the spatial points and the temporal lags involved are visualized.")
        for (i in 1:couples.nrow) {
          if(cova.nv[i,1] != "-"  && cova.nv[i,2] != "-"){

                print(cova.nv[i,])

          }
        }

    }
    }


    #= This is to compute C(0,0) as a variance
    #= of all data cova00<-var(as.vector(array.matdata.sel), na.rm =TRUE)

    #= Compute the covariance C00 =#
    cova00_vec <- matrix(0, nrow = length(selstaz), ncol = 1)
    for (i in 1:length(selstaz)) {
      cova00_vec[i, 1] <- var(array.matdata.sel[, , i], na.rm = TRUE)

    }
    cova00 <- mean(cova00_vec, na.rm = TRUE)

    #========================================================================#
    #= Compute the spatial and temporal marginal covariance                 =#
    #========================================================================#
    nflag.cova.h.nv <- 0
    nflag.cova.u.nv <- 0

    if (type.test == 1 || type.test == 2 || type.test == 3) {
      #= Compute the spatial marginal covariance =#

      cova.h <- matrix(0, nrow(couples), 1)
      cova.h.nv <- matrix(NA, nrow(couples), 2)
      for (i in 1:nrow(couples)) {


        cova.h[i, ] <- cov(array.matdata.sel[, , couples[i, 1]], array.matdata.sel[,
                                                                                   , couples[i, 2]], use = "pairwise.complete.obs")

        if(cova.h[i, ] < 0){
          nflag.cova.h.nv <- nflag.cova.h.nv + 1
          cova.h.nv[i,] <- selstaz[couples[i,1:2]]
        }
      }

      if (nflag.cova.h.nv != 0) {
        message("Warning message: ", nflag.cova.h.nv, " negative spatial covariances are detected.")
        message("In the following the spatial points and the temporal lags involved are visualized.")
          for (i in 1:couples.nrow) {
            if(is.na(cova.h.nv[i, 1]) == FALSE && is.na(cova.h.nv[i, 2]) == FALSE){
              print(cova.h.nv[i, ])
            }
        }

      }

      #= Compute the temporal marginal covariance =#
      nstaz <- length(selstaz)
      cova.u.ncol <- as.integer(couples.ncol/2)
      cova.u <- matrix(0, cova.u.ncol, 1)
      cova.u.nv <- matrix(NA, cova.u.ncol, 1)
      cova.ui <- matrix(0, nstaz, 1)
      jj <- -1
      for (j in 1:(couples.ncol/2)) {
        jj <- jj + 2
        i <- 1
        while (couples[i, jj + 2] == 0 && i <= (nrow(couples) - 1)) {
          i <- i + 1
        }  #= temporal marginal covariance C(0,u) is computed as a mean of the covariance  C_i(u), i=1,..nstaz
        if (i <= nrow(couples) && couples[i, jj + 2] != 0) {
          for (z in 1:length(selstaz)) {

            cova.ui[z, ] <- cov(array.matdata.sel[-(nrow(matdata.sel) -
                                                      couples[i, jj + 2] + 1:nrow(matdata.sel)), , z], array.matdata.sel[-(1:couples[i,
                                                                                                                                     jj + 2]), , z], use = "pairwise.complete.obs")

          }
        }

        cova.u[j, ] <- mean(cova.ui, na.rm = TRUE)
        if(cova.u[j, ] < 0){
          nflag.cova.u.nv <- nflag.cova.u.nv + 1
          cova.u.nv[j,] <- couples[i,jj+2]
        }
      }

    }

    if (nflag.cova.u.nv != 0) {
      message("Warning message: ", nflag.cova.u.nv, " negative temporal covariances are detected.")
      message("In the following the spatial points and the temporal lags involved are visualized.")
      for (i in 1:cova.u.ncol) {
          if(is.na(cova.u.nv[i, ]) == FALSE){
            print(cova.u.nv[i, ])
          }
      }

    }

    #========================================================================#
    #= End of computing the spatial and temporal marginal covariance        =#
    #========================================================================#

    #=====================================#
    #= Start computing f(G) and matrix B =#
    #=====================================#

    #========================================================================#
    #= Start computing f(G) and matrix B for type.test = 1, 2                =#
    #========================================================================#
    if (type.test == 1 || type.test == 2) {

      #= f(G) =#

      f.cova <- matrix(0, nrow = nct + (couples.ncol/2), ncol = 1)
      couples.ncol.r <- 0
      for (i in 1:couples.nrow) {
        couples.nrow.r <- 0
        nf <- 0
        for (j in 1:couples.ncol) {
          if (couples[i, j + 2] != 0) {
            couples.nrow.r <- couples.nrow.r + 1
            cov.n <- as.integer(couples.nrow.r + couples.ncol.r)
            f.cova[cov.n, ] <- cova[cov.n, ]/cova.h[i, ]
          }
        }
        if (nf == 0) {
          couples.ncol.r <- couples.ncol.r + couples.nrow.r
        }
        nf <- 1
      }
      for (i in 1:(couples.ncol/2)) {
        f.cova[nct + i, ] <- cova.u[i, ]/cova00
      }

      #= Compute the matrix B =#

      B <- matrix(0, nrow = (1 + nct + nrow(couples) + (couples.ncol/2)),
                  ncol = nct + (couples.ncol/2))

      for (i in (nct + 1):(nct + couples.ncol/2)) {

        B[1, i] <- -f.cova[i, 1]/cova00

      }

      for (i in 1:nct) {

        B[i + 1, i] <- f.cova[i, 1]/cova[i, 1]
     }

      jj <- 0
      couples.ncol.r <- 0
      for (i in 1:couples.nrow) {
        couples.nrow.r <- 0
        nf <- 0
        for (j in 1:couples.ncol) {
          if (couples[i, j + 2] != 0) {
            couples.nrow.r <- couples.nrow.r + 1
            cov.n <- as.integer(couples.nrow.r + couples.ncol.r)
            jj <- jj + 1
            B[i + 1 + nct, jj] <- -f.cova[cov.n, 1]/cova.h[i, ]
          }
        }
        if (nf == 0) {
          couples.ncol.r <- couples.ncol.r + couples.nrow.r
        }
        nf <- 1

      }
      for (i in 1:(couples.ncol/2)) {
        B[i + 1 + nct + nrow(couples), nct + i] <- 1/cova00

      }

    }
    #=======================================================#
    #= End computing f(G) and matrix B for type.test =1, 2  =#
    #=======================================================#


    #========================================================================#
    #= Start computing matrix B for type.test =3 (f(G) is not required)      =#
    #========================================================================#
    if (type.test == 3) {

      #= Compute the matrix B =#

      B <- matrix(0, nrow = nct + nrow(couples) + (couples.ncol/2), ncol = 2 *
                    nct)

      ii <- -1

      for (i in 1:nct) {

        ii <- ii + 2

        B[i, ii] <- 1
        B[i, ii + 1] <- 1

      }

      jj <- -1
      couples.ncol.r <- 0
      for (i in 1:couples.nrow) {
        couples.nrow.r <- 0
        for (j in 1:couples.ncol) {
          if (couples[i, j + 2] != 0) {
            jj <- jj + 2
            B[i + nct, jj] <- -1
          }
        }
      }
      couples.ncol.r <- 0
      kk <- -1
      for (j in 1:((couples.ncol)/2)) {
        kk <- kk + 2
        lagt <- sort(couples[, kk + 2], decreasing = TRUE)
        jj <- -1
        for (k in 1:couples.nrow) {
          for (i in 1:couples.ncol) {
            if (couples[k, i + 2] != 0) {
              jj <- jj + 2
            }
            if (couples[k, i + 2] == lagt[1]) {
              B[j + couples.nrow + nct, jj + 1] <- -1
            }
          }

        }
      }

    }
    #==================================================================#
    #= End computing matrix B for type.test =3 (f(G) is not required)  =#
    #==================================================================#

    #========================================================================#
    #= End computing f(G) and matrix B                                      =#
    #========================================================================#

    #======================#
    #= Start computing G  =#
    #======================#

    if (type.test == 1 || type.test == 2) {
      cova <- rbind(cova00, cova, cova.h, cova.u)
      row.names(cova) <- NULL
    }
    if (type.test == 3) {
      cova <- rbind(cova, cova.h, cova.u)
    }

    #====================#
    #= End computing G  =#
    #====================#
  }
  #========================================================================#
  #= end test on type.test=0, 1, 2, 3                                      =#
  #========================================================================#

  #========================================================================#
  #= Start building the contrast matrix                                   =#
  #========================================================================#

  #= start test on type.test=0 =#
  if (type.test == 0) {
    nct <- 0
    for (i in 1:nrow(couples)) {
      for (j in 3:ncol(couples)) {
        if (couples[i, j] != 0) {
          nct <- nct + 1
        }
      }
    }
    A.0.nrow <- as.integer(nct/2)
    A.0 <- matrix(0, nrow = A.0.nrow, ncol = nct)
    n2 <- 1
    for (i in 1:A.0.nrow) {
      A.0[i, n2] <- (1)
      A.0[i, n2 + 1] <- (-1)
      n2 <- n2 + 2
    }
  }
  #= end test on type.test=0 =#

  #= start test on type.test=1 and type.test=2 =#
  if (type.test == 1 || type.test == 2) {
    nct <- 0

    for (i in 1:nrow(couples)) {
      for (j in 3:ncol(couples)) {
        if (couples[i, j] != 0) {
          nct <- nct + 1
        }
      }
    }

    couples.ncol <- as.integer((ncol(couples) - 2))

    A.1 <- matrix(0, nrow = nct, ncol = nct + (couples.ncol/2))

    #= in order to include 1 =#

    for (i in 1:nct) {
      A.1[i, i] <- 1
    }

    #= in order to include -1 =#
    jj <- 0
    for (i in 1:nrow(couples)) {
      for (j in 1:couples.ncol) {
        if (couples[i, j + 2] != 0) {
          jj <- jj + 1
          kk <- as.integer((j - 1)/2) + 1
          A.1[jj, kk + nct] <- -1
        }
      }
    }
  }

  #= end test on type.test=1 and type.test=2 =#

  #============================================#
  #= Start test on type.test=2(with partition) =#
  #= (to be completed)                        =#
  #============================================#
  # pdomain <- 0
  # if(type.test==2 && pdomain==1){
  # message('For this type of test, building the contrast matrix is not required
  #          to perform the specific test.') }
  #  End test on type.test=2(with partition)
  # (to be completed)

  #= start test on type.test=3 =#

  if (type.test == 3) {
    nct <- 0

    for (i in 1:nrow(couples)) {
      for (j in 3:ncol(couples)) {
        if (couples[i, j] != 0) {
          nct <- nct + 1
        }
      }
    }

    couples.ncol <- as.integer((ncol(couples) - 2))
    A.3 <- matrix(0, nrow = 2 * nct, ncol = nct + nrow(couples) + (couples.ncol/2))

    #= in order to include 1 =#

    jj <- -1
    for (i in 1:nct) {
      jj <- jj + 2
      A.3[jj, i] <- 1
      A.3[jj + 1, i] <- 1
    }

    #= in order to include -1=#
    jj <- -1
    for (i in 1:nrow(couples)) {
      for (j in 1:couples.ncol) {
        if (couples[i, j + 2] != 0) {
          jj <- jj + 2
          kk <- as.integer((j - 1)/2) + 1
          A.3[jj, i + nct] <- -1
          A.3[jj + 1, kk + nct + nrow(couples)] <- -1
        }
      }
    }
    A.3bis.nrow <- as.integer(nct)
    A.3bis <- matrix(0, nrow = A.3bis.nrow, ncol = nct * 2)
    n2 <- 1
    for (i in 1:A.3bis.nrow) {
      A.3bis[i, n2] <- (1)
      A.3bis[i, n2 + 1] <- (-1)
      n2 <- n2 + 2
    }
  }
  #= end test on type.test=3 =#



  #==========================================#
  #= Start defining the class covastat      =#
  #==========================================#

  if (type.test == 0) {

    cova.h <- matrix(NA, 1, 1)
    cova.u <- matrix(NA, 1, 1)
    f.G <- matrix(NA, 1, 1)
    B <- matrix(NA, 1, 1)
    A <- A.0
  }
  if (type.test == 1 || type.test == 2) {
    f.G <- f.cova
    A <- A.1
  }
  if (type.test == 3) {
    f.G <- matrix(NA, 1, 1)
    A <- A.3
  }


  new("covastat", G = cova, cova.h = cova.h, cova.u = cova.u, f.G = f.G,
      B = B, A = A, typetest = typetest)


  #= End defining the class covastat =#
}
#' @include sepindex.R couples.R blocks.R covablocks.R
NULL
#' @param object object of class \code{covastat} for method \code{show}
#' @rdname covastat-class
#' @aliases covastat-class
#' @aliases covastat-method
#' @aliases show
#' @export
setMethod(f="show", signature="covastat", definition=function(object) {
  cat("An object of class covastat", "\n")
  cat("\n")
  cat("Slot 'G':")
  cat("\n")
  print(object@G)
  cat("\n")
  cat("Slot 'cova.h':")
  cat("\n")
  if(object@typetest == "sym"){
    print("This slot is not available for the required typetest")
  }else{
    print(object@cova.h)
  }
  cat("\n")
  cat("Slot 'cov.u':")
  cat("\n")
  if(object@typetest == "sym"){
    print("This slot is not available for the required typetest")
  }else{
    print(object@cova.u)
  }
  cat("\n")
  cat("Slot 'f.G':")
  cat("\n")
  if(object@typetest == "sym"){
    print("This slot is not available for the required typetest")
  }else{
    print(object@f.G)
  }
  cat("\n")
  cat("Slot 'B':")
  cat("\n")
  if(object@typetest == "sym"){
    print("This slot is not available for the required typetest")
  }else{
    print(object@B)
  }
  cat("\n")
  cat("Slot 'A':")
  cat("\n")
  print(object@A)
  cat("\n")
  cat("Slot 'typetest':")
  cat("\n")
  print(object@typetest)
}
)
