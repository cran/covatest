#' Class "covastat"
#'
#' A class for the sample spatio-temporal covariances for the specified spatial
#' and temporal lags, given in \code{stpairs} (object of class \code{couple}).
#' Depending on the type of test, the empirical variance, the sample spatial
#' and temporal marginal covariances are also computed.
#'
#' @slot G matrix; containing the spatio-temporal covariances for the specified
#' lags. For all tests, except for the symmetry test (\code{typetest=0}), the
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
#' @slot beta.data vector; containing the different values of the parameter beta,
#' available only for the test on the Gneiting class of model (\code{typetest=5})
#' @slot typetest numeric; contains the code of the test to be performed
#'
#' @note {
#' A stop occurs if the number of spatial points fixed in \code{stpairs}
#' (object of class \code{couples}) is less than 2.
#' \itemize{
#' \item If \code{typetest} is equal to 0 (symmetry test) \code{cova.h}, \code{cova.u},
#' \code{f.G} and \code{B} are not available
#'
#' \item If \code{typetest} is equal to 4 (test on the integrated product class
#' of models) \code{cova.h} and \code{cova.u} are not available
#'
#' \item If \code{typetest} is equal to 5 (test on the Gneiting class of models),
#' \code{cova.h} is not available
#' }
#' }
#'
#' @rdname covastat-class
#' @exportClass covastat
setClass("covastat", slots = c(G = "matrix",
                               cova.h = "matrix",
                               cova.u = "matrix",
                               f.G = "array",
                               B = "array",
                               A = "matrix",
                               typetest = "numeric",
                               beta.data = "ANY"))

#' @param matdata STFDF/STSDF or \code{data.frame}; which contains the
#' coordinates of the spatial points, the identification code of the spatial
#' points, the indentification code of the temporal points and the values of
#' the variable, typically output from \code{dataprep}
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
#' @param typetest integer; set \code{typetest=0} for symmetry test (default
#' choice), \code{typetest=1} for separability test, \code{typetest=2} for type
#' of non separability test, \code{typetest=3} for the test on the product-sum
#' class of models, \code{typetest=4} for the test on the integrated product
#' class of models, \code{typetest=5} for the test on the Gneiting class of
#' models
#' @param beta.data vector; this argument is required only for \code{typetest=5},
#' otherwise it has to be set equal to NULL (default choice). It contains the
#' different values of the parameter beta, which can assume values in the range 0-1
#'
#' @details
#' A message appears on the user's console if the \code{G} vector
#' contains spatio-temporal negative covariances. The message returns the negative
#' value/values and it will help to identify the spatial and the temporal lags
#' involved.
#'
#' @seealso \linkS4class{couples}
#' @seealso \code{\link{dataprep}}
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
#' Cappello, C., De Iaco, S., Posa, D., 2017, Testing the type of
#' non-separability and some classes of covariance models for space-time data.
#' Stochastic Environmental Research and Risk Assessment,
#' doi 10.1007/s00477-017-1472-2
#'
#' @examples
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
#' t.couples.in = t.couples.in.sym, typetest = 0, typecode = character())
#'
#' covast.sym <- covastat(matdata = rr_13, pardata1 = 1, pardata2 = 1,
#' stpairs = couples.sym, typetest = 0, beta.data = NULL)
#'
#' ###method for covastat
#' #1. show
#' covast.sym
#'
#' @rdname covastat-class
#' @export
covastat <- function(matdata, pardata1, pardata2, stpairs, typetest = 0, beta.data = NULL) {

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {abs(x - round(x)) <
      tol}

  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}

  ### SOME CHECKS ON THE ARGUMENTS ###

  if (is.scalar(pardata1) == FALSE || is.scalar(pardata2) == FALSE) {
    stop("Some of the arguments are not numeric. Stop running")
  }


  if(pardata1 != as.integer(pardata1) || pardata2 != as.integer(pardata2)){

    pardata1 <- as.integer(pardata1)
    pardata2 <- as.integer(pardata2)
    warning("The arguments expected to be integer are forced to be integer numbers.")
  }


  if (!inherits(stpairs, "couples")){
    stop("stpairs argument has to be of class couples")
  }

  if(typetest == 5 && !is.numeric(beta.data) == TRUE){
    stop("beta.data argument has to be numeric")
  }

  if(typetest == 5 && is.null(beta.data) == TRUE){
    stop("beta.data argument has to be set")
  }

  if(typetest == 5 && max(beta.data) > 1 ){
    stop("beta.data argument has to be in the interval 0-1")
  }

  if(typetest == 5 && min(beta.data) < 0 ){
    stop("beta.data argument has to be in the interval 0-1")
  }

  if(stpairs@typetest != typetest){
    warning("Reminder: the argument typetest is different from the one defined in stpairs")
  }

  if(stpairs@typetest == 0 && typetest>=3){
    stop("The argument typetest is not consistent with respect to the one defined in stpairs. Please change typetest or define a new stpairs")
  }

   if(any(stpairs@tl.couples < 0) && typetest>=3){
     stop("The argument typetest is not consistent with respect to the one defined in stpairs. Please change typetest or define a new stpairs")
   }



  selstaz <- stpairs@sel.staz
  couples <- stpairs@couples.st
  nstaz <- length(selstaz)

  #========================================================================#
  #= type of tests: 0=symmetry (default choice), 1=separability,          =#
  #= 2= type of non separability, 3=type of variability,                  =#
  #= 4-7= type of model (4= product-sum, 5= integrated product            =#
  #= 6= Gneiting, 7= Cressie-Huang)                                       =#
  #========================================================================#


  if (is.scalar(typetest) == FALSE || typetest < 0 || typetest > 5) {
    stop("The argument for typetest is not admissible.")
  }

  if (typetest >= 3) {
    typetest <- typetest + 1
  }

  #========================================================================#
  #= Select data for the stations of interest                             =#
  #========================================================================#

  if (class(matdata) == "matrix" || class(matdata) == "data.frame") {
    iclsp.id <- as.integer(pardata1)
    iclvr <- as.integer(pardata2)
  }

  if (is.vector(selstaz) && length(selstaz) >= 2) {
    for (i in 1:length(selstaz)) {

      #= data in GSLIB format =#
      if (class(matdata) == "matrix" || class(matdata) == "data.frame") {

        if (is.numeric(matdata[, iclvr]) == TRUE) {

          if (i == 1) {
            selstaz.names <- matdata[, iclsp.id]
            selstaz.inter <- intersect(selstaz.names, selstaz)
            if (length(selstaz.inter) != length(selstaz)) {
              stop("No data for some of the selected spatial points.Please go back to the function 'couples' and revise the vector of the selected spatial points")

            }
            }

          datistaz <- matdata[matdata[, iclsp.id] == selstaz[i], iclvr]
        } else {
          stop("Check the column in which the values of the
               variable are stored. Data must be numeric")
        }

        } else {

          #= data in gstat format =#
          if (class(matdata) == "STFDF") {
            if (i == 1) {
              selstaz.names <- row.names(matdata@sp)
              selstaz.inter <- intersect(selstaz.names, selstaz)
              if (length(selstaz.inter) != length(selstaz)) {
                stop("No data for some of the selected spatial points.Please go back to the function 'couples' and revise the vector of the selected spatial points")

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
                    stop("No data for some of the selected spatial points.Please go back to the function 'couples' and revise the vector of the selected spatial points")

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
              stop("The class of data must be matrix (gslib format),
                   data.frame, STSDF or STSDF")
            }
            }
          }


      if (i == 1) {
        matdata.sel <- datistaz
      }
      if (i > 1) {
        matdata.sel <- cbind(matdata.sel, datistaz)
      }


        }

  } else {
    stop("The number of spatial points selected in function 'couples' must be a vector with at least two components")
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
  #= 3=variability, 4-7=typemodel)                                           =#
  #===========================================================================#

  #========================================================================#
  #= start test on typetest=0, 1,2,3 and 4-7                              =#
  #========================================================================#

    nflag.cova.nv <- 0
  if (typetest == 0 || typetest == 1 || typetest == 2 || typetest == 3 ||
      typetest >= 4) {

    #= Compute spatio-temporal covariance (with hs and ht different from zero)
    #= for each block =#

    if (typetest == 0 || typetest == 1 || typetest == 2 || typetest == 3) {
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
            }

            if (couples[i, j + 2] < 0) {
              cova[cov.n, ] <- cov(array.matdata.sel[-(1:(-couples[i,
                                                                   j + 2])), , couples[i, 1]], array.matdata.sel[-(nrow(matdata.sel) +
                                                                                                                     couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                                                                                         2]], use = "pairwise.complete.obs")
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





    if (nflag.cova.nv != 0) {
      message(nflag.cova.nv, " negative spatio-temporal covariance/es detected.")
      message("In the following the spatial points and the temporal lags involved are visualized.")
        for (i in 1:couples.nrow) {
          if(cova.nv[i,1] != "-"  && cova.nv[i,2] != "-"){

                print(cova.nv[i,])

          }
        }

    }
    }

        if (typetest >= 4) {
      cova <- matrix(0, nrow = nct, ncol = 1)

      couples.ncol <- as.integer((ncol(couples) - 2))
      covacouple <- matrix(0, nrow = nrow(couples), ncol = couples.ncol)

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
            }


            if (couples[i, j + 2] < 0) {
              cova[cov.n, ] <- cov(array.matdata.sel[-(1:(-couples[i,
                                                                   j + 2])), , couples[i, 1]], array.matdata.sel[-(nrow(matdata.sel) +
                                                                                                                     couples[i, j + 2] + 1:nrow(matdata.sel)), , couples[i,
                                                                                                                                                                         2]], use = "pairwise.complete.obs")
            }
            covacouple[i, j] <- cova[cov.n, ]

            if (cova[cov.n, ] < 0) {
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

      if (nflag.cova.nv != 0) {
        message(nflag.cova.nv, " negative spatio-temporal covariance/es detected.")
        message("In the following the spatial points and the temporal lags involved are visualized.")
          for (i in 1:couples.nrow) {
            if(cova.nv[i,1] != "-" && cova.nv[i,2] != "-"){
              print(cova.nv[i,])
          }
        }

      }

      #= Check on columns and rows with non-zero values
      count_nozero_col <- matrix(0, nrow = (nrow(couples)/3), ncol = ((ncol(couples) -
                                                                         2)/2))
      count_nozero_row <- matrix(0, nrow = (nrow(couples)/3), ncol = ((ncol(couples) -
                                                                         2)/2))
      ii <- -2
      jj <- 0
      iflag <- matrix(0, nrow = (nrow(couples)/3), ncol = (((ncol(couples) -
                                                               2)/3)/2))
      for (i in 1:(nrow(couples)/3)) {
        kk <- -5
        ii <- ii + 3
        ki <- -2
        for (k in 1:(((ncol(couples) - 2)/3)/2)) {
          kk <- kk + 6
          ki <- ki + 3
          col_tlag <- c((kk + 2), (kk + 4), (kk + 6))
          count_nozero_col[i, ki:(ki + 2)] <- colSums(couples[ii:(ii +
                                                                    2), col_tlag] != 0)
          count_nozero_row[i, ki:(ki + 2)] <- rowSums(couples[ii:(ii +
                                                                    2), col_tlag] != 0)

          if (match(3, count_nozero_col[i, ki:(ki + 2)], nomatch = 0) !=
              0 && match(3, count_nozero_row[i, ki:(ki + 2)], nomatch = 0) !=
              0) {
            iflag[i, k] <- 1
            jj <- jj + 1
          }

        }
      }
      if (jj == 0) {
        stop("No temporal lags have been specified.")
      }


      #= End check on columns and rows with non-zero values

      #= Start permutation of each column of mat.cova:store sequentally the
      #= covariance for each set of comparisons (at least 3 spatial lags and 3
      #= temporal lags)


      cova2 <- cova
      cova <- matrix(0, nrow = nrow(cova), ncol = ncol(cova))

      z2 <- 1
      z <- 1
      #= read the elements of count_nozero_row

      for (ii in 1:(nrow(couples)/3)) {

        for (jj in 1:(((ncol(couples) - 2)/2)/3)) {


          if (iflag[ii, jj] == 1) {
            z3 <- z
            for (j in 1:3) {
              count.skip <- 0

              if (jj == 1 && j > 1) {
                for (zz in 1:(((ncol(couples) - 2)/2)/3)) {
                  if (jj != zz) {
                    count.skip <- count.skip + count_nozero_row[ii,
                                                                j - 1 + (zz - 1) * 3]

                  }
                }
              }
              if (jj != 1) {
                for (zz in 1:(((ncol(couples) - 2)/2)/3)) {
                  if (zz < jj) {
                    count.skip <- count.skip + count_nozero_row[ii,
                                                                j + (zz - 1) * 3]

                  }
                  if (zz > jj && j > 1) {
                    count.skip <- count.skip + count_nozero_row[ii,
                                                                j - 1 + (zz - 1) * 3]

                  }
                }
              }
              cova[z2:(z2 + count_nozero_row[ii, j + (jj - 1) * 3] -
                         1), ] <- cova2[(count.skip + z3):(count.skip + z3 +
                                                             count_nozero_row[ii, j + (jj - 1) * 3] - 1), ]
              z2 <- z2 + count_nozero_row[ii, j + (jj - 1) * 3]
              z3 <- z3 + count_nozero_row[ii, j + (jj - 1) * 3] + count.skip


            }



          }

        }
        z <- z2
      }


      #= End permutation of each column of mat.cova:store sequentally the
      #= covariance for each set of comparisons (at least 3 spatial lags and 3
      #= temporal lags)



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

    if (typetest == 1 || typetest == 2 || typetest == 3) {
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
        message(nflag.cova.h.nv, " negative spatial covariances have been detected.")
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
      message(nflag.cova.u.nv, " negative temporal covariances have been detected.")
      message("In the following the spatial points and the temporal lags involved are visualized.")
      for (i in 1:cova.u.ncol) {
          if(is.na(cova.u.nv[i, ]) == FALSE){
            print(cova.u.nv[i, ])
          }
      }

    }

    lstaz_zero <- 0

    if (typetest >= 4) {

      #============================================================================#
      #= Detect spatial points non used for comparisons (corresponding to rows of =#
      #= zeros in count_nozero_col)                                               =#
      #============================================================================#
      nspaz <- 0
      for (i in 1:(nrow(couples)/3)) {

        if (sum(count_nozero_col[i, ] != 0) == 0) {


          lstaz_zero <- lstaz_zero + length(unique(couples[(1 + (i -
                                                                   1) * 3):(3 + (i - 1) * 3), 1:2]))

        }
      }
      staz_zero <- as.vector(matrix(0, nrow = lstaz_zero, ncol = 1))
      for (i in 1:(nrow(couples)/3)) {


        if (sum(count_nozero_col[i, ] != 0) != 0) {
          nspaz <- nspaz + 1
        }
        lstaz <- 0
        if (sum(count_nozero_col[i, ] != 0) == 0) {

          lstaz_zero <- lstaz + length(unique(couples[(1 + (i - 1) *
                                                         3):(3 + (i - 1) * 3), 1:2]))
          staz_zero[(lstaz + 1):lstaz_zero] <- as.vector(unique(couples[(1 +
                                                                           (i - 1) * 3):(3 + (i - 1) * 3), 1:2]))
          lstaz <- lstaz_zero


        }
      }
      #= End detecting spatial points non used for comparisons (corresponding to
      #= rows of zeros in count_nozero_col)

      #= Compute the spatial marginal covariance =#
      cova.hsel <- matrix(0, nspaz * 3, 1)
      cova.h <- matrix(0, nrow(couples), 1)
      cova.h.nv <- matrix(NA, nspaz * 3, 2)
      j <- 0
      for (i in 1:nrow(couples)) {
        cova.h[i, ] <- 0

        if (sum(count_nozero_col[as.integer(1 + ((i - 1)/3)), ] != 0) !=
            0) {
          j <- j + 1
          cova.hsel[j, ] <- cov(array.matdata.sel[, , couples[i, 1]],
                                array.matdata.sel[, , couples[i, 2]], use = "pairwise.complete.obs")
          cova.h[i, ] <- cova.hsel[j, ]
          if(cova.hsel[j, ] < 0){
            nflag.cova.h.nv <- nflag.cova.h.nv + 1
            cova.h.nv[j,] <- selstaz[couples[i,1:2]]
          }
        }

      }

      if (nflag.cova.h.nv != 0) {
        message(nflag.cova.h.nv, " negative spatial covariances have been detected.")
        message("In the following the spatial points and the temporal lags involved are visualized.")
        for (i in 1:(nspaz * 3)) {
            if(is.na(cova.h.nv[i, 1]) == FALSE && is.na(cova.h.nv[i, 2]) == FALSE){
              print(cova.h.nv[i, ])
            }
        }

      }

      #= Compute the temporal marginal covariance =#
      nstaz <- length(selstaz)

      ntemp <- 0

      for (i in 1:(couples.ncol/2)) {

        if (sum(count_nozero_col[, i] != 0) != 0) {
          ntemp <- ntemp + 1
        }
      }

      nstaz <- length(selstaz)
      cova.u.ncol <- as.integer(couples.ncol/2)
      cova.u <- matrix(0, cova.u.ncol, 1)
      cova.ui <- matrix(0, nstaz, 1)

      cova.u.ncolsel <- ntemp
      cova.usel <- matrix(0, cova.u.ncolsel, 1)
      cova.usel.nv <- matrix(NA, cova.u.ncolsel, 1)
      cova.uisel <- matrix(0, length(setdiff(selstaz, selstaz[staz_zero])),
                           1)

      jjj <- 0
      jj <- -1
      zz <- 0
      for (j in 1:(couples.ncol/2)) {

        jj <- jj + 2
        i <- 1
        while (couples[i, jj + 2] == 0 && i <= (nrow(couples) - 1)) {
          i <- i + 1
        }  #= temporal marginal covariance C(0,u) is computed as a mean of the covariance  C_i(u), i=1,..nstaz
        if (i <= nrow(couples) && couples[i, jj + 2] != 0) {
          zz <- 0
          for (z in 1:length(selstaz)) {

            if (length(staz_zero) == 0) {
              cova.uisel[z, ] <- cov(array.matdata.sel[-(nrow(matdata.sel) -
                                                           couples[i, jj + 2] + 1:nrow(matdata.sel)), , z], array.matdata.sel[-(1:couples[i,
                                                                                                                                          jj + 2]), , z], use = "pairwise.complete.obs")
            } else {

              if (length(intersect(selstaz[z], selstaz[staz_zero])) ==
                  0) {
                zz <- zz + 1

                cova.uisel[zz, ] <- cov(array.matdata.sel[-(nrow(matdata.sel) -
                                                              couples[i, jj + 2] + 1:nrow(matdata.sel)), , z], array.matdata.sel[-(1:couples[i,
                                                                                                                                             jj + 2]), , z], use = "pairwise.complete.obs")
              }
            }
          }
        }

        cova.u[j, ] <- 0
        if (sum(count_nozero_col[, j] != 0) != 0) {
          jjj <- jjj + 1
          cova.usel[jjj, ] <- mean(cova.uisel, na.rm = TRUE)
          cova.u[j, ] <- cova.usel[jjj, ]
          if(cova.u[j, ] < 0){
            nflag.cova.u.nv <- nflag.cova.u.nv + 1
            cova.u.nv[j,] <- stpairs@tl.couples[j]
          }
        }
      }

      if (nflag.cova.u.nv != 0) {
        message(nflag.cova.u.nv, " negative temporal covariances have been detected.")
        message("In the following the spatial points and the temporal lags involved are visualized.")
          for (i in 1:cova.u.ncol) {
            if(is.na(cova.u.nv[i, ]) == FALSE){
              print(cova.u.nv[i, ])
            }
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
    #= Start computing f(G) and matrix B for typetest = 1, 2                =#
    #========================================================================#
    if (typetest == 1 || typetest == 2) {

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
    #= End computing f(G) and matrix B for typetest =1, 2  =#
    #=======================================================#


    #========================================================================#
    #= Start computing matrix B for typetest =3 (f(G) is not required)      =#
    #========================================================================#
    if (typetest == 3) {

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
    #= End computing matrix B for typetest =3 (f(G) is not required)  =#
    #==================================================================#

    #========================================================================#
    #= Start computing f(G) and matrix B for typetest >= 4                  =#
    #========================================================================#

    if (typetest >= 4) {


      #========================================================================#
      #= Start computing f(G) and matrix B for typetest = 4                   =#
      #========================================================================#


      if (typetest == 4) {


        #= f(G) for typetest= 4 =#
        f.cova.sub <- matrix(0, nrow = 2, ncol = 1)
        jj <- 0
        ii <- -2
        #= iflag<-matrix(0, nrow=(nrow(sp.couples)/3), ncol=((n.temp/3)/2))
        flag1 <- 0
        flag2 <- 0
        for (i in 1:(nrow(couples)/3)) {
          kk <- -5
          ii <- ii + 3
          ki <- -2
          for (k in 1:(((ncol(couples) - 2)/3)/2)) {
            kk <- kk + 6
            ki <- ki + 3
            col_tlag <- c((kk + 2), (kk + 4), (kk + 6))
            #= count_nozero_col<- colSums(couples[ii:(ii+2), col_tlag] != 0)
            #= count_nozero_row<- rowSums(couples[ii:(ii+2), col_tlag] != 0)

            #= Start cicle on iflag which identifies the subblock with at least one row
            #= and one column full filled

            if (iflag[i, k] == 1) {
              #= jj<-jj+1

              zz <- ki - 1
              for (z in 1:3) {
                zz <- zz + 1
                if (count_nozero_col[i, zz] == 3) {
                  f.cova.sub[1, ] <- (covacouple[(ii + 1), col_tlag[z] -
                                                   2] - covacouple[(ii), col_tlag[z] - 2])/(cova.h[2 +
                                                                                                     (i - 1) * 3, ] - cova.h[1 + (i - 1) * 3, ])
                  f.cova.sub[2, ] <- (covacouple[(ii + 2), col_tlag[z] -
                                                   2] - covacouple[(ii + 1), col_tlag[z] - 2])/(cova.h[3 +
                                                                                                         (i - 1) * 3, ] - cova.h[2 + (i - 1) * 3, ])
                  if (flag1 == 1) {
                    f.cova1 <- rbind(f.cova1, f.cova.sub)
                  }
                  if (flag1 == 0) {
                    f.cova1 <- f.cova.sub
                    flag1 <- 1
                  }
                }


              }
              zz <- ki - 1
              for (z in 1:3) {
                zz <- zz + 1
                if (count_nozero_row[i, zz] == 3) {
                  f.cova.sub[1, ] <- (covacouple[(ii + z - 1), (kk +
                                                                  4 - 2)] - covacouple[(ii + z - 1), (kk + 2 - 2)])/(cova.u[2 +
                                                                                                                              (k - 1) * 3, ] - cova.u[1 + (k - 1) * 3, ])

                  f.cova.sub[2, ] <- (covacouple[(ii + z - 1), (kk +
                                                                  6 - 2)] - covacouple[(ii + z - 1), (kk + 4 - 2)])/(cova.u[3 +
                                                                                                                              (k - 1) * 3, ] - cova.u[2 + (k - 1) * 3, ])
                  if (flag2 == 1) {
                    f.cova2 <- rbind(f.cova2, f.cova.sub)
                  }
                  if (flag2 == 0) {
                    f.cova2 <- f.cova.sub
                    flag2 <- 1
                  }
                }

              }

            }

            #= end cicle on iflag which identifies the subblock with at least one row and
            #= one column full filled

          }
        }

        f.cova <- as.matrix(rbind(f.cova1, f.cova2), ncol = 1)


        #= Compute the matrix B for typetest = 4 =#

        count_colfull <- sum(count_nozero_col[, ] == 3)
        count_rowfull <- sum(count_nozero_row[, ] == 3)
        B <- matrix(0, nrow = (nct + nspaz * 3 + ntemp), ncol = ((count_colfull +
                                                                    count_rowfull) * 2))

        zz <- -1
        nct_block <- 0
        jj <- 0
        zz2 <- 0
        ii2 <- 0
        iii <- 0
        for (i in 1:(nrow(couples)/3)) {
          kk <- -5
          kkk <- 0

          zflag <- 0
          ki <- -2
          for (k in 1:(((ncol(couples) - 2)/3)/2)) {
            kk <- kk + 6
            ki <- ki + 3
            jflag <- 0

            #= Start cicle on iflag which identifies the subblock with at least one row
            #= and one column full filled

            if (iflag[i, k] == 1) {


              #= Start cicle on each subblock in order to read the non zero elements in
              #= couples

              jj <- nct_block


              ii <- -1 + count_colfull * 2 + ii2
              for (j in 1:3) {
                #= j in 1:3row because the test is applied on groups of 3 spatial lags
                zz <- -1 + zz2

                for (z in 1:3) {
                  #= z in 1:3column because the test is applied on groups of 3 temporal lags


                  if (couples[(j + (i - 1) * 3), (kk + (z - 1) * 2 +
                                                  2)] != 0) {
                    nct_block <- nct_block + 1
                    jj <- jj + 1
                    #= Check on full columns
                    if (count_nozero_col[i, z + ki - 1] == 3) {
                      zz <- zz + 2
                      if (j == 1) {

                        B[jj, zz] <- -1/(cova.h[2 + (i - 1) * 3, ] -
                                           cova.h[1 + (i - 1) * 3, ])

                        if (j == 1 && zflag == 0) {
                          iii <- iii + 1
                          zflag <- 1
                        }

                        B[nct + 1 + (i - 1) * 3, zz] <- (covacouple[(j +
                                                                       1 + (i - 1) * 3), (kk + (z - 1) * 2)] - covacouple[(j +
                                                                                                                             (i - 1) * 3), (kk + (z - 1) * 2)])/((cova.h[2 +
                                                                                                                                                                           (i - 1) * 3, ] - cova.h[1 + (i - 1) * 3, ])^2)

                        B[nct + 2 + (i - 1) * 3, zz] <- -(covacouple[(j +
                                                                        1 + (i - 1) * 3), (kk + (z - 1) * 2)] - covacouple[(j +
                                                                                                                              (i - 1) * 3), (kk + (z - 1) * 2)])/((cova.h[2 +
                                                                                                                                                                            (i - 1) * 3, ] - cova.h[1 + (i - 1) * 3, ])^2)


                        B[nct + 2 + (i - 1) * 3, zz + 1] <- (covacouple[(j +
                                                                           2 + (i - 1) * 3), (kk + (z - 1) * 2)] - covacouple[(j +
                                                                                                                                 1 + (i - 1) * 3), (kk + (z - 1) * 2)])/((cova.h[3 +
                                                                                                                                                                                   (i - 1) * 3, ] - cova.h[2 + (i - 1) * 3, ])^2)

                        B[nct + 3 + (i - 1) * 3, zz + 1] <- -(covacouple[(j +
                                                                            2 + (i - 1) * 3), (kk + (z - 1) * 2)] - covacouple[(j +
                                                                                                                                  1 + (i - 1) * 3), (kk + (z - 1) * 2)])/((cova.h[3 +
                                                                                                                                                                                    (i - 1) * 3, ] - cova.h[2 + (i - 1) * 3, ])^2)
                      }

                      if (j == 2) {

                        B[jj, zz] <- 1/(cova.h[2 + (i - 1) * 3, ] -
                                          cova.h[1 + (i - 1) * 3, ])
                        B[jj, zz + 1] <- -1/(cova.h[3 + (i - 1) * 3,
                                                    ] - cova.h[2 + (i - 1) * 3, ])
                      }

                      if (j == 3) {

                        B[jj, zz + 1] <- 1/(cova.h[3 + (i - 1) * 3,
                                                   ] - cova.h[2 + (i - 1) * 3, ])
                      }
                    }


                    #= Check on full rows

                    if (count_nozero_row[i, j + ki - 1] == 3) {

                      if (z == 1) {

                        ii <- ii + 2
                        B[jj, ii] <- -1/(cova.u[2 + (k - 1) * 3, ] -
                                           cova.u[1 + (k - 1) * 3, ])
                        if (z == 1 && jflag == 0) {
                          kkk <- kkk + 1
                          jflag <- 1
                        }
                        B[nct + nspaz * 3 + 1 + (k - 1) * 3, ii] <- (covacouple[(j +
                                                                                   (i - 1) * 3), (kk + 2 + (z - 1) * 2)] - covacouple[(j +
                                                                                                                                         (i - 1) * 3), (kk + (z - 1) * 2)])/((cova.u[2 +
                                                                                                                                                                                       (k - 1) * 3, ] - cova.u[1 + (k - 1) * 3, ])^2)

                        B[nct + nspaz * 3 + 2 + (k - 1) * 3, ii] <- -(covacouple[(j +
                                                                                    (i - 1) * 3), (kk + 2 + (z - 1) * 2)] - covacouple[(j +
                                                                                                                                          (i - 1) * 3), (kk + (z - 1) * 2)])/((cova.u[2 +
                                                                                                                                                                                        (k - 1) * 3, ] - cova.u[1 + (k - 1) * 3, ])^2)


                        B[nct + nspaz * 3 + 2 + (k - 1) * 3, ii + 1] <- (covacouple[(j +
                                                                                       (i - 1) * 3), (kk + 4 + (z - 1) * 2)] - covacouple[(j +
                                                                                                                                             (i - 1) * 3), (kk + 2 + (z - 1) * 2)])/((cova.u[3 +
                                                                                                                                                                                               (k - 1) * 3, ] - cova.u[2 + (k - 1) * 3, ])^2)

                        B[nct + nspaz * 3 + 3 + (k - 1) * 3, ii + 1] <- -(covacouple[(j +
                                                                                        (i - 1) * 3), (kk + 4 + (z - 1) * 2)] - covacouple[(j +
                                                                                                                                              (i - 1) * 3), (kk + 2 + (z - 1) * 2)])/((cova.u[3 +
                                                                                                                                                                                                (k - 1) * 3, ] - cova.u[2 + (k - 1) * 3, ])^2)
                      }

                      if (z == 2) {

                        B[jj, ii] <- 1/(cova.u[2 + (k - 1) * 3, ] -
                                          cova.u[1 + (k - 1) * 3, ])
                        B[jj, ii + 1] <- -1/(cova.u[3 + (k - 1) * 3,
                                                    ] - cova.u[2 + (k - 1) * 3, ])
                      }

                      if (z == 3) {

                        B[jj, ii + 1] <- 1/(cova.u[3 + (k - 1) * 3,
                                                   ] - cova.u[2 + (k - 1) * 3, ])
                      }
                    }



                  }
                }



              }

              ii2 <- ii - count_colfull * 2 + 1

            }


            zz2 <- zz + 1
            #= end cicle on iflag which identifies the subblock with at least one row and
            #= one column full filled

          }
        }


      }
      #========================================================================#
      #= End computing f(G) and matrix B for typetest=4                       =#
      #========================================================================#

      #========================================================================#
      #= Start computing f(G) and matrix B for typetest = 5                   =#
      #========================================================================#

      if (typetest == 5) {

        #= f(G) for typetest= 5 =#


        f.cova.sub <- matrix(0, nrow = 2, ncol = 1)
        jj <- 0
        ii <- -2
        #= iflag<-matrix(0, nrow=(nrow(sp.couples)/3), ncol=((n.temp/3)/2))
        flag1 <- 0
        flag2 <- 0
        for (i in 1:(nrow(couples)/3)) {
          kk <- -5
          ii <- ii + 3
          ki <- -2
          for (k in 1:(((ncol(couples) - 2)/3)/2)) {
            kk <- kk + 6
            ki <- ki + 3
            col_tlag <- c((kk + 2), (kk + 4), (kk + 6))
            #= count_nozero_col<- colSums(couples[ii:(ii+2), col_tlag] != 0)
            #= count_nozero_row<- rowSums(couples[ii:(ii+2), col_tlag] != 0)

            if (iflag[i, k] == 1) {
              #= jj<-jj+1

              zz <- ki - 1
              for (z in 1:3) {
                zz <- zz + 1
                if (count_nozero_col[i, zz] == 3) {
                  f.cova.sub[1, ] <- 1/(covacouple[(ii + 1), col_tlag[z] -
                                                     2]) - 1/(covacouple[(ii), col_tlag[z] - 2])

                  f.cova.sub[2, ] <- 1/(covacouple[(ii + 2), col_tlag[z] -
                                                     2]) - 1/(covacouple[(ii + 1), col_tlag[z] - 2])

                  if (flag1 == 1) {
                    f.cova1 <- rbind(f.cova1, f.cova.sub)
                  }
                  if (flag1 == 0) {
                    f.cova1 <- f.cova.sub
                    flag1 <- 1
                  }
                }


              }
              zz <- ki - 1
              for (z in 1:3) {
                zz <- zz + 1
                if (count_nozero_row[i, zz] == 3) {

                  f.cova.sub[1, ] <- 1/(covacouple[(ii + z - 1), (kk +
                                                                    4 - 2)]) - 1/(covacouple[(ii + z - 1), (kk + 2 -
                                                                                                              2)])

                  f.cova.sub[2, ] <- 1/(covacouple[(ii + z - 1), (kk +
                                                                    6 - 2)]) - 1/(covacouple[(ii + z - 1), (kk + 4 -
                                                                                                              2)])

                  if (flag2 == 1) {
                    f.cova2 <- rbind(f.cova2, f.cova.sub)
                  }
                  if (flag2 == 0) {
                    f.cova2 <- f.cova.sub
                    flag2 <- 1
                  }

                }

              }




            }
            #= end cicle on iflag which identifies the subblock with at least one row and
            #= one column full filled

          }
        }

        f.cova <- as.matrix(rbind(f.cova1, f.cova2), ncol = 1)


        #= Compute the matrix B for typetest = 5 =#


        count_colfull <- sum(count_nozero_col[, ] == 3)
        count_rowfull <- sum(count_nozero_row[, ] == 3)
        B <- matrix(0, nrow = (nct), ncol = ((count_colfull + count_rowfull) *
                                               2))

        zz <- -1
        nct_block <- 0
        jj <- 0
        zz2 <- 0
        ii2 <- 0

        for (i in 1:(nrow(couples)/3)) {
          kk <- -5
          ki <- -2
          for (k in 1:(((ncol(couples) - 2)/3)/2)) {
            kk <- kk + 6
            ki <- ki + 3


            #= Start cicle on iflag which identifies the subblock with at least one row
            #= and one column full filled

            if (iflag[i, k] == 1) {

              #= Start cicle on each subblock in order to read the non zero elements in
              #= couples

              jj <- nct_block


              ii <- -1 + count_colfull * 2 + ii2
              for (j in 1:3) {
                #= j in 1:3row because the test is applied on groups of 3 spatial lags
                zz <- -1 + zz2

                for (z in 1:3) {
                  #= z in 1:3column because the test is applied on groups of 3 temporal lags


                  if (couples[(j + (i - 1) * 3), (kk + (z - 1) * 2 +
                                                  2)] != 0) {
                    nct_block <- nct_block + 1
                    jj <- jj + 1
                    #= Check on full columns


                    if (count_nozero_col[i, z + ki - 1] == 3) {
                      zz <- zz + 2
                      if (j == 1) {

                        B[jj, zz] <- 1/((covacouple[(j + (i - 1) * 3),
                                                    (kk + (z - 1) * 2)])^2)


                      }

                      if (j == 2) {

                        B[jj, zz] <- -1/((covacouple[(j + (i - 1) *
                                                        3), (kk + (z - 1) * 2)])^2)
                        B[jj, zz + 1] <- 1/((covacouple[(j + (i - 1) *
                                                           3), (kk + (z - 1) * 2)])^2)

                      }

                      if (j == 3) {

                        B[jj, zz + 1] <- -1/((covacouple[(j + (i - 1) *
                                                            3), (kk + (z - 1) * 2)])^2)

                      }
                    }

                    #= Check on full rows

                    if (count_nozero_row[i, j + ki - 1] == 3) {

                      if (z == 1) {

                        ii <- ii + 2
                        B[jj, ii] <- 1/((covacouple[(j + (i - 1) * 3),
                                                    (kk + (z - 1) * 2)])^2)

                      }

                      if (z == 2) {

                        B[jj, ii] <- -1/((covacouple[(j + (i - 1) *
                                                        3), (kk + (z - 1) * 2)])^2)
                        B[jj, ii + 1] <- 1/((covacouple[(j + (i - 1) *
                                                           3), (kk + (z - 1) * 2)])^2)
                      }

                      if (z == 3) {

                        B[jj, ii + 1] <- -1/((covacouple[(j + (i - 1) *
                                                            3), (kk + (z - 1) * 2)])^2)
                      }
                    }


                  }
                }


              }

              ii2 <- ii - count_colfull * 2 + 1

            }


            zz2 <- zz + 1
            #= end cicle on iflag which identifies the subblock with at least one row and
            #= one column full filled

          }
        }


      }
      #========================================================================#
      #= End computing f(G) and matrix B for typetest=5                       =#
      #========================================================================#

      #========================================================================#
      # Start computing f(G) and matrix B for typetest = 6                    =#
      #========================================================================#

      if (typetest == 6) {

        nbeta <- length(beta.data)
        count_colfull <- sum(count_nozero_col[, ] == 3)
        count_rowfull <- sum(count_nozero_row[, ] == 3)
        arrayB <- array(0, dim = c(nct + ntemp, ((count_colfull + count_rowfull) *
                                                   2), nbeta))
        arrayF <- array(0, dim = c(((count_colfull + count_rowfull) *
                                      2), 1, nbeta))
        for (ciclebeta in 1:nbeta) {

          beta <- round(beta.data[ciclebeta], digits = 1)

          #= f(G) for typetest= 6 =#




          f.cova.sub <- matrix(0, nrow = 2, ncol = 1)
          jj <- 0
          ii <- -2
          flag1 <- 0
          flag2 <- 0
          for (i in 1:(nrow(couples)/3)) {
            kk <- -5
            ii <- ii + 3
            ki <- -2
            for (k in 1:(((ncol(couples) - 2)/3)/2)) {
              kk <- kk + 6
              ki <- ki + 3
              col_tlag <- c((kk + 2), (kk + 4), (kk + 6))

              if (iflag[i, k] == 1) {
                #= jj<-jj+1

                zz <- ki - 1
                for (z in 1:3) {
                  zz <- zz + 1
                  if (count_nozero_col[i, zz] == 3) {
                    f.cova.sub[1, ] <- ((covacouple[(ii + 1), col_tlag[z] -
                                                      2]))/(covacouple[(ii), col_tlag[z] - 2])

                    f.cova.sub[2, ] <- ((covacouple[(ii + 2), col_tlag[z] -
                                                      2]))/(covacouple[(ii + 1), col_tlag[z] - 2])


                    if (flag1 == 1) {
                      f.cova1 <- rbind(f.cova1, f.cova.sub)
                    }
                    if (flag1 == 0) {
                      f.cova1 <- f.cova.sub
                      flag1 <- 1
                    }
                  }


                }
                zz <- ki - 1
                for (z in 1:3) {
                  zz <- zz + 1
                  if (count_nozero_row[i, zz] == 3) {


                    f.cova.sub[1, ] <- (log(cova.u[2 + (k - 1) * 3,
                                                   ]/(covacouple[(ii + z - 1), (kk + 4 - 2)]), base = exp(1)))^(-2/(beta)) -
                      (log(cova.u[1 + (k - 1) * 3, ]/(covacouple[(ii +
                                                                    z - 1), (kk + 2 - 2)]), base = exp(1)))^(-2/(beta))

                    f.cova.sub[2, ] <- (log(cova.u[3 + (k - 1) * 3,
                                                   ]/(covacouple[(ii + z - 1), (kk + 6 - 2)]), base = exp(1)))^(-2/(beta)) -
                      (log(cova.u[2 + (k - 1) * 3, ]/(covacouple[(ii +
                                                                    z - 1), (kk + 4 - 2)]), base = exp(1)))^(-2/(beta))



                    if (flag2 == 1) {
                      f.cova2 <- rbind(f.cova2, f.cova.sub)
                    }
                    if (flag2 == 0) {
                      f.cova2 <- f.cova.sub
                      flag2 <- 1
                    }
                  }

                }




              }
              #= end cicle on iflag which identifies the subblock with at least one row and
              #= one column full filled

            }
          }
          f.cova <- as.matrix(rbind(f.cova1, f.cova2), ncol = 1)

          #= Compute the matrix B for typetest = 6 =#

          count_colfull <- sum(count_nozero_col[, ] == 3)
          count_rowfull <- sum(count_nozero_row[, ] == 3)

          B <- matrix(0, nrow = (nct + ntemp), ncol = ((count_colfull +
                                                          count_rowfull) * 2))

          zz <- -1
          nct_block <- 0
          jj <- 0
          zz2 <- 0
          ii2 <- 0
          iii <- 0
          for (i in 1:(nrow(couples)/3)) {
            kk <- -5
            kkk <- 0

            zflag <- 0
            ki <- -2
            for (k in 1:(((ncol(couples) - 2)/3)/2)) {
              kk <- kk + 6
              ki <- ki + 3
              jflag <- 0

              #= Start cicle on iflag which identifies the subblock with at least one row
              #= and one column full filled

              if (iflag[i, k] == 1) {


                #= Start cicle on each subblock in order to read the non zero elements in
                #= couples

                jj <- nct_block


                ii <- -1 + count_colfull * 2 + ii2
                for (j in 1:3) {
                  #= j in 1:3row because the test is applied on groups of 3 spatial lags
                  zz <- -1 + zz2

                  for (z in 1:3) {
                    #= z in 1:3column because the test is applied on groups of 3 temporal lags


                    if (couples[(j + (i - 1) * 3), (kk + (z - 1) * 2 +
                                                    2)] != 0) {
                      nct_block <- nct_block + 1
                      jj <- jj + 1

                      #= Check on full columns
                      if (count_nozero_col[i, z + ki - 1] == 3) {
                        zz <- zz + 2
                        if (j == 1) {

                          B[jj, zz] <- -(covacouple[(j + 1 + (i - 1) *
                                                       3), (kk + (z - 1) * 2)])/((covacouple[(j +
                                                                                                (i - 1) * 3), (kk + (z - 1) * 2)])^2)


                        }

                        if (j == 2) {

                          B[jj, zz] <- 1/((covacouple[(j - 1 + (i -
                                                                  1) * 3), (kk + (z - 1) * 2)]))
                          B[jj, zz + 1] <- -(covacouple[(j + 1 + (i -
                                                                    1) * 3), (kk + (z - 1) * 2)])/((covacouple[(j +
                                                                                                                  (i - 1) * 3), (kk + (z - 1) * 2)])^2)
                        }

                        if (j == 3) {

                          B[jj, zz + 1] <- 1/(covacouple[(j - 1 + (i -
                                                                     1) * 3), (kk + (z - 1) * 2)])
                        }
                      }

                      #= Check on full rows

                      if (count_nozero_row[i, j + ki - 1] == 3) {

                        if (z == 1) {

                          ii <- ii + 2
                          B[jj, ii] <- -(2/beta) * ((log(cova.u[1 +
                                                                  (k - 1) * 3, ]/covacouple[(j + (i - 1) *
                                                                                               3), (kk + (z - 1) * 2)], base = exp(1)))^(-(2/beta) -
                                                                                                                                           1)) * (1/covacouple[(j + (i - 1) * 3), (kk +
                                                                                                                                                                                     (z - 1) * 2)])

                          #= if(z==1&&jflag==0){ kkk<-kkk+1 jflag<-1 }
                          B[nct + 1 + (k - 1) * 3, ii] <- (2/beta) *
                            ((log(cova.u[1 + (k - 1) * 3, ]/covacouple[(j +
                                                                          (i - 1) * 3), (kk + (z - 1) * 2)], base = exp(1)))^(-(2/beta) -
                                                                                                                                1)) * (1/cova.u[1 + (k - 1) * 3, ])

                          B[nct + 2 + (k - 1) * 3, ii] <- -(2/beta) *
                            ((log(cova.u[2 + (k - 1) * 3, ]/covacouple[(j +
                                                                          (i - 1) * 3), (kk + 2 + (z - 1) * 2)],
                                  base = exp(1)))^(-(2/beta) - 1)) * (1/cova.u[2 +
                                                                                 (k - 1) * 3, ])


                          B[nct + 2 + (k - 1) * 3, ii + 1] <- (2/beta) *
                            ((log(cova.u[2 + (k - 1) * 3, ]/covacouple[(j +
                                                                          (i - 1) * 3), (kk + 2 + (z - 1) * 2)],
                                  base = exp(1)))^(-(2/beta) - 1)) * (1/cova.u[2 +
                                                                                 (k - 1) * 3, ])

                          B[nct + 3 + (k - 1) * 3, ii + 1] <- -(2/beta) *
                            ((log(cova.u[3 + (k - 1) * 3, ]/covacouple[(j +
                                                                          (i - 1) * 3), (kk + 4 + (z - 1) * 2)],
                                  base = exp(1)))^(-(2/beta) - 1)) * (1/cova.u[3 +
                                                                                 (k - 1) * 3, ])
                        }

                        if (z == 2) {

                          B[jj, ii] <- (2/beta) * ((log(cova.u[2 + (k -
                                                                      1) * 3, ]/covacouple[(j + (i - 1) * 3),
                                                                                           (kk + (z - 1) * 2)], base = exp(1)))^(-(2/beta) -
                                                                                                                                   1)) * (1/covacouple[(j + (i - 1) * 3), (kk +
                                                                                                                                                                             (z - 1) * 2)])
                          B[jj, ii + 1] <- -(2/beta) * ((log(cova.u[2 +
                                                                      (k - 1) * 3, ]/covacouple[(j + (i - 1) *
                                                                                                   3), (kk + (z - 1) * 2)], base = exp(1)))^(-(2/beta) -
                                                                                                                                               1)) * (1/covacouple[(j + (i - 1) * 3), (kk +
                                                                                                                                                                                         (z - 1) * 2)])
                        }

                        if (z == 3) {

                          B[jj, ii + 1] <- (2/beta) * ((log(cova.u[3 +
                                                                     (k - 1) * 3, ]/covacouple[(j + (i - 1) *
                                                                                                  3), (kk + (z - 1) * 2)], base = exp(1)))^(-(2/beta) -
                                                                                                                                              1)) * (1/covacouple[(j + (i - 1) * 3), (kk +
                                                                                                                                                                                        (z - 1) * 2)])
                        }
                      }



                    }
                  }


                }

                ii2 <- ii - count_colfull * 2 + 1

              }


              zz2 <- zz + 1
              #= end cicle on iflag which identifies the subblock with at least one row and
              #= one column full filled

            }
          }




          arrayB[, , ciclebeta] <- B[, ]
          arrayF[, , ciclebeta] <- f.cova[, ]

        }
      }
      #========================================================================#
      #= End computing f(G) and matrix B for typetest=6                       =#
      #========================================================================#

      #========================================================================#
      #= Start computing f(G) and matrix B for typetest = 7                   =#
      #========================================================================#

      if (typetest == 7) {

        #= f(G) for typetest= 7 =#

        f.cova.sub <- matrix(0, nrow = 2, ncol = 1)
        jj <- 0
        ii <- -2
        #= iflag<-matrix(0, nrow=(nrow(sp.couples)/3), ncol=((n.temp/3)/2))
        flag1 <- 0
        flag2 <- 0
        for (i in 1:(nrow(couples)/3)) {
          kk <- -5
          ii <- ii + 3
          ki <- -2
          for (k in 1:(((ncol(couples) - 2)/3)/2)) {
            kk <- kk + 6
            ki <- ki + 3
            col_tlag <- c((kk + 2), (kk + 4), (kk + 6))



            if (iflag[i, k] == 1) {


              zz <- ki - 1
              for (z in 1:3) {
                zz <- zz + 1
                if (count_nozero_col[i, zz] == 3) {

                  f.cova.sub[1, ] <- ((covacouple[(ii + 1), col_tlag[z] -
                                                    2]))/((covacouple[(ii), col_tlag[z] - 2]))


                  if (z == 3) {


                  }
                  f.cova.sub[2, ] <- ((covacouple[(ii + 2), col_tlag[z] -
                                                    2]))/((covacouple[(ii + 1), col_tlag[z] - 2]))



                  if (flag1 == 1) {
                    f.cova1 <- rbind(f.cova1, f.cova.sub)
                  }
                  if (flag1 == 0) {
                    f.cova1 <- f.cova.sub
                    flag1 <- 1
                  }
                }


              }
              zz <- ki - 1
              for (z in 1:3) {
                zz <- zz + 1
                if (count_nozero_row[i, zz] == 3) {


                  f.cova.sub[1, ] <- (log(cova.u[2 + (k - 1) * 3, ]/(covacouple[(ii +
                                                                                   z - 1), (kk + 4 - 2)]), base = exp(1)))^(-1) - (log(cova.u[1 +
                                                                                                                                                (k - 1) * 3, ]/(covacouple[(ii + z - 1), (kk + 2 -
                                                                                                                                                                                            2)]), base = exp(1)))^(-1)





                  f.cova.sub[2, ] <- (log(cova.u[3 + (k - 1) * 3, ]/(covacouple[(ii +
                                                                                   z - 1), (kk + 6 - 2)]), base = exp(1)))^(-1) - (log(cova.u[2 +
                                                                                                                                                (k - 1) * 3, ]/(covacouple[(ii + z - 1), (kk + 4 -
                                                                                                                                                                                            2)]), base = exp(1)))^(-1)


                  if (flag2 == 1) {
                    f.cova2 <- rbind(f.cova2, f.cova.sub)
                  }
                  if (flag2 == 0) {
                    f.cova2 <- f.cova.sub
                    flag2 <- 1
                  }
                }

              }





            }
            #= end cicle on iflag which identifies the subblock with at least one row and
            #= one column full filled

          }
        }

        f.cova <- as.matrix(rbind(f.cova1, f.cova2), ncol = 1)



        #= Compute the matrix B for typetest = 7 =#

        count_colfull <- sum(count_nozero_col[, ] == 3)
        count_rowfull <- sum(count_nozero_row[, ] == 3)

        B <- matrix(0, nrow = (nct + ntemp), ncol = ((count_colfull +
                                                        count_rowfull) * 2))

        zz <- -1
        nct_block <- 0
        jj <- 0
        zz2 <- 0
        ii2 <- 0

        for (i in 1:(nrow(couples)/3)) {
          kk <- -5
          kkk <- 0

          iii <- 0
          ki <- -2
          for (k in 1:(((ncol(couples) - 2)/3)/2)) {
            kk <- kk + 6
            ki <- ki + 3
            zflag <- 0
            jflag <- 0

            #= Start cicle on iflag which identifies the subblock with at least one row
            #= and one column full filled

            if (iflag[i, k] == 1) {


              #= Start cicle on each subblock in order to read the non zero elements in
              #= couples

              jj <- nct_block


              ii <- -1 + count_colfull * 2 + ii2
              for (j in 1:3) {
                #= j in 1:3row because the test is applied on groups of 3 spatial lags
                zz <- -1 + zz2

                for (z in 1:3) {
                  #= z in 1:3column because the test is applied on groups of 3 temporal lags


                  if (couples[(j + (i - 1) * 3), (kk + (z - 1) * 2 +
                                                  2)] != 0) {
                    nct_block <- nct_block + 1
                    jj <- jj + 1
                    #= Check on full columns
                    if (count_nozero_col[i, z + ki - 1] == 3) {

                      zz <- zz + 2
                      if (j == 1) {

                        B[jj, zz] <- -(covacouple[(j + 1 + (i - 1) *
                                                     3), (kk + (z - 1) * 2)])/((covacouple[(j +
                                                                                              (i - 1) * 3), (kk + (z - 1) * 2)])^2)


                      }

                      if (j == 2) {

                        B[jj, zz] <- 1/((covacouple[(j - 1 + (i - 1) *
                                                       3), (kk + (z - 1) * 2)]))
                        B[jj, zz + 1] <- -(covacouple[(j + 1 + (i -
                                                                  1) * 3), (kk + (z - 1) * 2)])/((covacouple[(j +
                                                                                                                (i - 1) * 3), (kk + (z - 1) * 2)])^2)
                      }

                      if (j == 3) {

                        B[jj, zz + 1] <- 1/(covacouple[(j - 1 + (i -
                                                                   1) * 3), (kk + (z - 1) * 2)])
                      }
                    }


                    #= Check on full rows

                    if (count_nozero_row[i, j + ki - 1] == 3) {

                      if (z == 1) {

                        ii <- ii + 2
                        B[jj, ii] <- -((log(cova.u[1 + (k - 1) * 3,
                                                   ]/covacouple[(j + (i - 1) * 3), (kk + (z -
                                                                                            1) * 2)], base = exp(1)))^(-2)) * (1/covacouple[(j +
                                                                                                                                               (i - 1) * 3), (kk + (z - 1) * 2)])

                        B[nct + 1 + (k - 1) * 3, ii] <- ((log(cova.u[1 +
                                                                       (k - 1) * 3, ]/covacouple[(j + (i - 1) * 3),
                                                                                                 (kk + (z - 1) * 2)], base = exp(1)))^(-2)) *
                          (1/cova.u[1 + (k - 1) * 3, ])

                        B[nct + 2 + (k - 1) * 3, ii] <- -((log(cova.u[2 +
                                                                        (k - 1) * 3, ]/covacouple[(j + (i - 1) * 3),
                                                                                                  (kk + 2 + (z - 1) * 2)], base = exp(1)))^(-2)) *
                          (1/cova.u[2 + (k - 1) * 3, ])


                        B[nct + 2 + (k - 1) * 3, ii + 1] <- ((log(cova.u[2 +
                                                                           (k - 1) * 3, ]/covacouple[(j + (i - 1) * 3),
                                                                                                     (kk + 2 + (z - 1) * 2)], base = exp(1)))^(-2)) *
                          (1/cova.u[2 + (k - 1) * 3, ])

                        B[nct + 3 + (k - 1) * 3, ii + 1] <- -((log(cova.u[3 +
                                                                            (k - 1) * 3, ]/covacouple[(j + (i - 1) * 3),
                                                                                                      (kk + 4 + (z - 1) * 2)], base = exp(1)))^(-2)) *
                          (1/cova.u[3 + (k - 1) * 3, ])
                      }

                      if (z == 2) {

                        B[jj, ii] <- ((log(cova.u[2 + (k - 1) * 3, ]/covacouple[(j +
                                                                                   (i - 1) * 3), (kk + (z - 1) * 2)], base = exp(1)))^(-2)) *
                          (1/covacouple[(j + (i - 1) * 3), (kk + (z -
                                                                    1) * 2)])
                        B[jj, ii + 1] <- -((log(cova.u[2 + (k - 1) *
                                                         3, ]/covacouple[(j + (i - 1) * 3), (kk + (z -
                                                                                                     1) * 2)], base = exp(1)))^(-2)) * (1/covacouple[(j +
                                                                                                                                                        (i - 1) * 3), (kk + (z - 1) * 2)])
                      }

                      if (z == 3) {

                        B[jj, ii + 1] <- ((log(cova.u[3 + (k - 1) *
                                                        3, ]/covacouple[(j + (i - 1) * 3), (kk + (z -
                                                                                                    1) * 2)], base = exp(1)))^(-2)) * (1/covacouple[(j +
                                                                                                                                                       (i - 1) * 3), (kk + (z - 1) * 2)])
                      }
                    }



                  }
                }


              }

              ii2 <- ii - count_colfull * 2 + 1

            }


            zz2 <- zz + 1
            #= end cicle on iflag which identifies the subblock with at least one row and
            #= one column full filled

          }
        }


      }
      #========================================================================#
      #= End computing f(G) and matrix B for typetest=7                       =#
      #========================================================================#

    }
    #========================================================================#
    #= End computing f(G) and matrix B                                      =#
    #========================================================================#


    #======================#
    #= Start computing G  =#
    #======================#

    if (typetest == 1 || typetest == 2) {
      cova <- rbind(cova00, cova, cova.h, cova.u)
      row.names(cova) <- NULL
    }
    if (typetest == 3) {
      cova <- rbind(cova, cova.h, cova.u)
    }
    if (typetest == 4) {
      cova <- rbind(cova, cova.hsel, cova.usel)
    }
    if (typetest == 5) {
      cova <- rbind(cova)
    }
    if (typetest == 6 || typetest == 7) {
      cova <- rbind(cova, cova.usel)
    }

    #====================#
    #= End computing G  =#
    #====================#
  }
  #========================================================================#
  #= end test on typetest=0, 1, 2, 3 and 4-7                              =#
  #========================================================================#



  #========================================================================#
  #= Start building the contrast matrix                                   =#
  #========================================================================#

  #= start test on typetest=0 =#


  if (typetest == 0) {
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
  #= end test on typetest=0 =#


  #= start test on typetest=1 and typetest=2 =#
  if (typetest == 1 || typetest == 2) {
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

  #= end test on typetest=1 and typetest=2 =#

  #============================================#
  #= Start test on typetest=2(with partition) =#
  #= (to be completed)                        =#
  #============================================#
  # pdomain <- 0
  # if(typetest==2 && pdomain==1){
  # message('For this type of test, building the contrast matrix is not required
  #          to perform the specific test.') }
  #  End test on typetest=2(with partition)
  # (to be completed)


  #= start test on typetest=3 =#

  if (typetest == 3) {
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
  #= end test on typetest=3 =#

  #==============================#
  #= Start test on typetest>=4  =#
  #==============================#

  if (typetest >= 4) {

    #= Check on columns and rows with non-zero values
    count_nozero_col <- matrix(0, nrow = (nrow(couples)/3), ncol = ((ncol(couples) -
                                                                       2)/2))
    count_nozero_row <- matrix(0, nrow = (nrow(couples)/3), ncol = ((ncol(couples) -
                                                                       2)/2))
    ii <- -2


    for (i in 1:(nrow(couples)/3)) {
      kk <- -5
      ii <- ii + 3
      ki <- -2
      for (k in 1:(((ncol(couples) - 2)/3)/2)) {
        kk <- kk + 6
        ki <- ki + 3
        col_tlag <- c((kk + 2), (kk + 4), (kk + 6))
        count_nozero_col[i, ki:(ki + 2)] <- colSums(couples[ii:(ii +
                                                                  2), col_tlag] != 0)
        count_nozero_row[i, ki:(ki + 2)] <- rowSums(couples[ii:(ii +
                                                                  2), col_tlag] != 0)

      }
    }


    count_colfull <- sum(count_nozero_col[, ] == 3)
    count_rowfull <- sum(count_nozero_row[, ] == 3)

    #= End check on columns and rows with non-zero values

    #= if(typetest>=4){
    nct <- (count_colfull + count_rowfull) * 2



    A.4.nrow <- as.integer(nct/2)

    A.4 <- matrix(0, nrow = A.4.nrow, ncol = nct)

    n2 <- 1
    for (i in 1:A.4.nrow) {
      A.4[i, n2] <- (1)
      A.4[i, n2 + 1] <- (-1)
      n2 <- n2 + 2

    }

  }
  #= end test on typetest>=4 =#

  #==========================================#
  #= Start defining the class covastat      =#
  #==========================================#

  if (typetest == 0) {

    cova.h <- matrix(NA, 1, 1)
    cova.u <- matrix(NA, 1, 1)
    f.G <- matrix(NA, 1, 1)
    B <- matrix(NA, 1, 1)
    A <- A.0

  }

  if (typetest == 1 || typetest == 2) {
    f.G <- f.cova
    A <- A.1
  }



  if (typetest == 3) {


    f.G <- matrix(NA, 1, 1)
    A <- A.3
  }

  if (typetest == 4) {

    f.G <- f.cova
    cova.h <- cova.hsel
    cova.u <- cova.usel
    A <- A.4
  }
  if (typetest == 5) {

    f.G <- f.cova
    cova.h <- matrix(NA, 1, 1)
    cova.u <- matrix(NA, 1, 1)
    A <- A.4
  }

  if (typetest == 6) {

    cova.h <- matrix(NA, 1, 1)
    cova.u <- cova.usel
    f.G <- arrayF
    B <- arrayB
    A <- A.4
  }

  if (typetest == 7) {

    f.G <- f.cova
    cova.h <- matrix(NA, 1, 1)
    cova.u <- cova.usel
    A <- A.4
  }

  if (typetest >= 4) {
      typetest <- typetest -1
    }

  new("covastat", G = cova, cova.h = cova.h, cova.u = cova.u, f.G = f.G,
      B = B, A = A, beta.data = beta.data, typetest = typetest)


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
  if(object@typetest == 0 || object@typetest == 4 || object@typetest == 5){
    print("This slot is not available for the required typetest")
  }else{
    print(object@cova.h)
  }
  cat("\n")
  cat("Slot 'cov.u':")
  cat("\n")
  if(object@typetest == 0 || object@typetest == 4){
    print("This slot is not available for the required typetest")
  }else{
    print(object@cova.u)
  }
  cat("\n")
  cat("Slot 'f.G':")
  cat("\n")
  if(object@typetest == 0){
    print("This slot is not available for the required typetest")
  }else{
    print(object@f.G)
  }
  cat("\n")
  cat("Slot 'B':")
  cat("\n")
  if(object@typetest == 0){
    print("This slot is not available for the required typetest")
  }else{
    print(object@B)
  }
  cat("\n")
  cat("Slot 'A':")
  cat("\n")
  print(object@A)
  cat("\n")
  cat("Slot 'beta.data':")
  cat("\n")
  if(object@typetest != 5){
    print("This slot is not available for the required typetest")
  }else{
    print(object@beta.data)
  }
  cat("\n")
  cat("Slot 'typetest':")
  cat("\n")
  print(object@typetest)
}
)
