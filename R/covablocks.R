#' Class "covablocks"
#'
#' A class for the sample spatio-temporal covariances to be computed for each
#'block of data and for the selected spatial and temporal lags fixed
#' in \code{stpairs} (output from \code{couples}). Depending on the type of
#' test the empirical variance, the sample spatial and temporal marginal
#' covariances for each block of data are also computed. Moreover, the sample
#' covariances between the spatio-temporal covariances at the specified spatial
#' and temporal lags are determined.
#'
#' @slot mat.cova matrix of sample spatio-temporal covariances for each block,
#' computed for the spatial and temporal lags given in \code{stpairs}
#' (object of class \code{couples})
#' @slot mat.cova.h matrix of sample spatial marginal covariances
#' for the specified lags
#' @slot mat.cova.u matrix of sample temporal marginal covariances
#' for the specified lags
#' @slot mat.cova.cova matrix of sample covariances between space-time
#' covariances for each block, computed for the spatial and temporal lags given
#' in \code{stpairs} (object of class \code{couples})
#' @slot typetest character; contains the code of the test to be performed
#'
#' @rdname covablocks-class
#' @exportClass covablocks
setClass("covablocks", slots = c(mat.cova = "matrix",
                                 mat.cova.h = "matrix",
                                 mat.cova.u = "matrix",
                                 mat.cova.cova = "matrix",
                                 typetest = "character"))

#' @param stblocks object of class \code{blocks}
#'
#' @param stpairs object of class \code{couples}, containing the spatial
#' points and the corresponding temporal lags to be analyzed
#'
#' @param typetest character; set \code{typetest ="sym"} for symmetry test
#' (default choice), \code{typetest ="sep"} for separability test, \code{typetest ="tnSep"}
#' for type of non separability test,\code{typetest ="productSum"} for the test
#' on the product-sum class of models, \code{typetest ="intProduct"} for the test
#' on the integrated product class of models, \code{typetest ="gneiting"} for the
#' test on the Gneiting class of models
#'
#' @note {
#' \itemize{
#' \item If \code{typetest} is equal to \code{"sym"} (symmetry test) or \code{"intProduct"}
#' (test on the integrated product class of models) \code{mat.cova.h} and \code{mat.cova.u}
#' are not available
#'
#' \item If \code{typetest} is equal to \code{"gneiting"} (test on the Gneiting
#' class of models) \code{mat.cova.h} is not available
#'
#' \item If temporal lags in \code{stpairs} are not consistent with block length
#' (\code{lb}) in \code{stblocks} an error message will be returned
#'
#' \item If the proportion between the maximum temporal lag in \code{stpairs} and
#' the block length (\code{lb}) in \code{stblocks} is greater than 0.25 a warning
#' message will be returned since the covariance estimation might not be reliable
#' }
#' }
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
#' block.sym <- blocks(lb = 40, ls = 10, matdata = rr_13, pardata1 = 1, pardata2 = 1,
#' stpairs = couples.sym)
#'
#' covabl.sym <- covablocks(stblocks = block.sym, stpairs = couples.sym, typetest = "sym")
#'
#' ### method for covablock
#' #1. show
#' covabl.sym
#'
#' @seealso \linkS4class{blocks}
#' @seealso \linkS4class{couples}
#'
#' @references
#' Cappello, C., De Iaco, S., Posa, D., 2018, Testing the type of
#' non-separability and some classes of space-time covariance function models.
#' Stochastic Environmental Research and Risk Assessment,
#' \bold{32} 17--35
#'
#' Cappello, C., De Iaco, S., Posa, D., 2020, {covatest}: An {R} Package for
#' Selecting a Class of Space-Time Covariance Functions.
#' Journal of Statistical Software, \bold{94(1)} 1--42.
#'
#' De Iaco, S., Palma, M., Posa, D., 2016. A general procedure for selecting a
#' class of fully symmetric space-time covariance functions.
#' Environmentrics, \bold{27(4)} 212--224.
#'
#' Li, B., Genton, M.G., Sherman, M., 2007, A nonparametric assessment
#' of properties of spacetime covariance functions.
#' Journal of the American Statistical Association, \bold{102} 736--744.
#'
#'
#' @rdname covablocks-class
#' @export
covablocks <- function(stblocks, stpairs, typetest = "sym") {

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) <
    tol

  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}

  ### SOME CHECKS ON THE ARGUMENTS ###

  if (!inherits(stpairs, "couples")){
    message("Start error message. stpairs argument has to be of class couples.")
    stop("End error message. Stop running.")
  }

  if (!inherits(stblocks, "blocks")){
    message("Start error message. stblocks argument has to be of class blocks.")
    stop("End error message. Stop running.")
  }

  if(identical(stpairs@sel.staz,stblocks@sel.staz) == FALSE){
    message("Start error message. The objects stpairs and stblocks have not been defined for the same set of spatial points.")
    stop("End error message. Stop running.")
  }

  if(stpairs@typetest != typetest){
    message("Warning message: the argument typetest is different from the one defined in stpairs")
  }



  matblock <- stblocks@mat.block
  selstaz <- stpairs@sel.staz
  couples <- stpairs@couples.st
  nstaz <- length(selstaz)
  block.ncol <- as.integer(ncol(matblock)/nstaz)  # number of blocks for each station
  block.nrow <- nrow(matblock)  # number of elements for each block


  #=================================================================#
  #= type of tests: 0=symmetry (default choice), 1=separability,   =#
  #= 2= type of non separability, 3=type of variability            =#
  #= 4 up to 7 = type of model                                     =#
  #=================================================================#

  if (is.character(typetest) == FALSE) {
    message("Start error message. The argument for typetest is not admissible.")
    stop("End error message. Stop running.")
  }

  if (typetest != "sym" && typetest != "sep"  && typetest != "tnSep" && typetest != "productSum"
      && typetest != "intProduct" && typetest != "gneiting") {
    message("Start error message. The argument for typetest is not admissible.")
    stop("End error message. Stop running.")
  }

  if (typetest == "sym") {
    type.test <- 0
  }else{if (typetest == "sep"){
    type.test <- 1
  }else{if (typetest == "tnSep"){
    type.test <- 2
  }else{if (typetest == "productSum"){
    type.test <- 4
  }else{if (typetest == "intProduct"){
    type.test <- 5
  }else{type.test <- 6 #Gneiting
  }
  }
  }
  }
  }


  if(stpairs@typetest == "sym" && type.test >= 4){
    message("Start error message. The argument typetest is not consistent with respect to the one defined in stpairs. Please change typetest or define a new stpairs.")
    stop("End error message. Stop running.")
  }

  if(any(stpairs@tl.couples < 0) && type.test>=4){
    message("Start error message. The argument typetest is not consistent with respect to the one defined in stpairs. Please change typetest or define a new stpairs.")
    stop("End error message. Stop running.")
  }

  maxt <- max(abs(max(couples[, -c(1:2)])), abs(min(couples[, -c(1:2)])))



  if (maxt > nrow(matblock)) {
    message("Start error message. Temporal lags defined in function 'couples' are not consistent with block length (lb) fixed in function 'blocks'.")
    stop("End error message. Stop running.")
  }

  if (maxt > (nrow(matblock)/4)) {
    ratio.max.t <- round(maxt, 3)
    message("Warning message: the proportion ", ratio.max.t, " between the maximum temporal lag defined in function 'couples' and block length defined in function 'blocks' is greater than 0.25. The covariance estimation might not be reliable.")
  }


  #=================================================================#
  #= start test on type.test=0, 1, 2, 3 and 4 up to 7               =#
  #=================================================================#

  #= Compute spatio-temporal covariance (with hs and ht different from zero)  =#
  #= for each block =#
  couples.nrow <- nrow(couples)
  pdomain <- 0
  if (pdomain == 0) {
    nct <- 0

    for (i in 1:nrow(couples)) {
      for (j in 3:ncol(couples)) {
        if (couples[i, j] != 0) {
          nct <- nct + 1
        }
      }
    }

    couples.ncol <- as.integer((ncol(couples) - 2))

    matcov.ncol <- as.integer(nct)
    block.array <- array(matblock, c(block.nrow, block.ncol, nstaz))
    mat.cova <- matrix(0, block.ncol, matcov.ncol)
    vec.na <- matrix(NA, block.nrow, 1)

    couples.ncol.r <- 0
    info.na <- matrix(NA, 1, 5)
    for (i in 1:couples.nrow) {
      couples.nrow.r <- 0
      nf <- 0
      for (j in 1:couples.ncol) {

        if (couples[i, j + 2] != 0) {

          couples.nrow.r <- couples.nrow.r + 1
          matcov.n <- as.integer(couples.nrow.r + couples.ncol.r)

          flag.block <- 0
          for (l in 1:block.ncol) {

            if (couples[i, j + 2] > 0) {
              mat.cova[l, matcov.n] <- cov(block.array[-(nrow(matblock) -
                                                           couples[i, j + 2] + 1:nrow(matblock)), l, couples[i,
                                                                                                             1]], block.array[-(1:couples[i, j + 2]), l, couples[i,
                                                                                                                                                                 2]], use = "pairwise.complete.obs")

              if(is.na(mat.cova[l, matcov.n]) == TRUE){
              if(flag.block == 0){
              if(is.na(info.na[1,1]) == TRUE){
                info.na[1,1] <- couples[i,1]
                info.na[1,2] <- couples[i,2]
                info.na[1,5] <- couples[i,j +2]
                if(identical(vec.na,block.array[-(nrow(matblock) -
                         couples[i, j + 2] + 1:nrow(matblock)), l, couples[i, 1]]) == TRUE){info.na[1,3] <- 1}
                if(identical(vec.na,block.array[-(1:couples[i, j + 2]), l, couples[i, 2]]) == TRUE){info.na[1,4] <- 2}
              }else{
                info.na <- rbind(info.na, c(couples[i,1:2],NA,NA,NA))
                info.na[1,5] <- couples[i,j +2]
                if(identical(vec.na,block.array[-(nrow(matblock) -
                                                  couples[i, j + 2] + 1:nrow(matblock)), l, couples[i, 1]]) == FALSE){info.na[1,3] <- 1}
                if(identical(vec.na,block.array[-(1:couples[i, j + 2]), l, couples[i, 2]]) == FALSE){info.na[1,4] <- 2}


              }
              flag.block <- 1
                }
              }

              }

            if (couples[i, j + 2] < 0) {
              mat.cova[l, matcov.n] <- cov(block.array[-(1:(-couples[i,
                                                                     j + 2])), l, couples[i, 1]], block.array[-(nrow(matblock) +
                                                                                                                  couples[i, j + 2] + 1:nrow(matblock)), l, couples[i,
                                                                                                                                                                    2]], use = "pairwise.complete.obs")


              if(is.na(mat.cova[l, matcov.n]) == TRUE){
                if(flag.block == 0){
                if(is.na(info.na[1,1]) == TRUE){
                  info.na[1,1] <- couples[i,1]
                  info.na[1,2] <- couples[i,2]
                  info.na[1,5] <- couples[i,j +2]
                  if(identical(vec.na,block.array[-(1:(-couples[i,
                                                                j + 2])), l, couples[i, 1]]) == TRUE){info.na[1,3] <- 1}
                  if(identical(vec.na,block.array[-(nrow(matblock) +
                                                    couples[i, j + 2] + 1:nrow(matblock)), l, couples[i,
                                                                                                      2]]) == TRUE){info.na[1,4] <- 2}
                }else{
                  info.na <- rbind(info.na, c(couples[i,1:2],NA,NA,NA))
                  info.na[1,5] <- couples[i,j +2]
                  if(identical(vec.na,block.array[-(1:(-couples[i,
                                                                j + 2])), l, couples[i, 1]]) == FALSE){info.na[1,3] <- 1}
                  if(identical(vec.na,block.array[-(nrow(matblock) +
                                                    couples[i, j + 2] + 1:nrow(matblock)), l, couples[i,
                                                                                                      2]]) == FALSE){info.na[1,4] <- 2}


                }
                  flag.block <- 1
              }
              }
              }


          }

        }

      }
      if (nf == 0) {
        couples.ncol.r <- couples.ncol.r + couples.nrow.r
      }
      nf <- 1


    }

    #= Check on columns and rows with non-zero values =#

    if (type.test >= 4) {
      ii <- -2
      jj <- 0
      count_nozero_col <- matrix(0, nrow = (nrow(couples)/3), ncol = ((ncol(couples) -
                                                                         2)/2))
      count_nozero_row <- matrix(0, nrow = (nrow(couples)/3), ncol = ((ncol(couples) -
                                                                         2)/2))
      iflag <- matrix(0, nrow = (nrow(couples)/3), ncol = (((ncol(couples) -
                                                               2)/3)/2))
      for (i in 1:(nrow(couples)/3)) {
        kk <- -5
        ii <- ii + 3

        for (k in 1:(((ncol(couples) - 2)/3)/2)) {
          kk <- kk + 6
          col_tlag <- c((kk + 2), (kk + 4), (kk + 6))

          count_nozero_col[i, (1 + (k - 1) * 3):(3 + (k - 1) * 3)] <- colSums(couples[ii:(ii +
                                                                                            2), col_tlag] != 0)
          count_nozero_row[i, (1 + (k - 1) * 3):(3 + (k - 1) * 3)] <- rowSums(couples[ii:(ii +
                                                                                            2), col_tlag] != 0)

          if (match(3, count_nozero_col[i, (1 + (k - 1) * 3):(3 + (k -
                                                                   1) * 3)], nomatch = 0) != 0 && match(3, count_nozero_row[i,
                                                                                                                            (1 + (k - 1) * 3):(3 + (k - 1) * 3)], nomatch = 0) != 0) {
            iflag[i, k] <- 1
            jj <- jj + 1
          }
        }
      }
      if (jj == 0) {
        message("Start error message. No temporal lags have been specified.")
        stop("End error message. Stop running.")
      }

    }
    #= End check on columns and rows with non-zero values =#


      if(is.na(info.na[1,1]) == FALSE){
        message("There are no enough data for computing the covariance: spatial couples, #point in the couple non-valid, #point in the couple non-valid, temporal lag non-valid.")
        for (i in 1:length(info.na)){
          print(info.na[i,])
        }
        message("Start error message. Please exclude/change the non-valid spatial couples/points/temporal lag from the selection.")
        stop("End error message. Stop running.")
      }

    #==========================================================================#
    #= Start if on type of test (1=separability, 2=type of non separability,  =#
    #= 3=variability, 4-7=type of model)                                      =#
    #==========================================================================#

    if (type.test == 1 || type.test == 2 || type.test == 3 || type.test >= 4) {

      #= Compute C00

      #= C00 of blocks as variance###
      #= vec.cova00<-matrix(0,nrow=block.ncol, ncol=1)
      #= for(l in 1:block.ncol){ vec.cova00[l, ]<-var(as.vector( block.array[,
      #= l,1:nstaz]), na.rm = TRUE) }

      #== Compute the mean of the covariances C00 ==#
      block.array <- array(matblock, c(nrow(matblock), block.ncol, length(selstaz)))
      mat.cova00 <- matrix(0, nrow = block.ncol, ncol = length(selstaz))
      vec.cova00 <- matrix(0, nrow = block.ncol, ncol = 1)
      for (l in 1:block.ncol) {
        for (i in 1:nstaz) {
          mat.cova00[l, i] <- cov(block.array[, l, i], block.array[,
                                                                   l, i], use = "pairwise.complete.obs")


        }
        vec.cova00[l, ] <- mean(mat.cova00[l, ], na.rm = TRUE)
      }

      #=================================================================#
      #= Compute the spatial and temporal marginal covariance          =#
      #=================================================================#


      if (type.test == 1 || type.test == 2 || type.test == 3) {
        #= Compute the spatial marginal covariance for each block=#
        mat.cova.h <- matrix(0, block.ncol, nrow(couples))
        for (i in 1:nrow(couples)) {
          for (l in 1:block.ncol) {

            mat.cova.h[l, i] <- cov(block.array[, l, couples[i, 1]],
                                    block.array[, l, couples[i, 2]], use = "pairwise.complete.obs")


          }
        }

        #= Compute the temporal marginal covariance for each block =#
        mat.cova.u.ncol <- as.integer(couples.ncol/2)
        mat.cova.u <- matrix(0, block.ncol, mat.cova.u.ncol)
        mat.cova.ui <- matrix(0, block.ncol, nstaz)
        jj <- -1
        for (j in 1:(couples.ncol/2)) {
          jj <- jj + 2
          i <- 1
          while (couples[i, jj + 2] == 0) {
            i <- i + 1
          }
          #= temporal marginal covariance C(0,u) is computed as a mean of the
          #= covariance C_i(u), i=1,..nstaz


          for (z in 1:nstaz) {

            for (l in 1:block.ncol) {



              mat.cova.ui[l, z] <- cov(block.array[-(nrow(matblock) -
                                                       couples[i, jj + 2] + 1:nrow(matblock)), l, z], block.array[-(1:couples[i,
                                                                                                                              jj + 2]), l, z], use = "pairwise.complete.obs")



            }


          }

          for (l in 1:block.ncol) {
            mat.cova.u[l, j] <- mean(mat.cova.ui[l, ], na.rm = TRUE)
          }

        }

      }

      if (type.test >= 4) {


        #= Detect spatial points non used for comparisons (corresponding to rows of
        #= zeros in count_nozero_col)
        lstaz_zero <- 0
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

        #= Compute the spatial marginal covariance for each block =#

        mat.cova.h <- matrix(0, block.ncol, nrow(couples))
        mat.cova.hsel <- matrix(0, block.ncol, nspaz * 3)
        j <- 0
        for (i in 1:nrow(couples)) {
          mat.cova.h[, i] <- 0

          if (sum(count_nozero_col[as.integer(1 + ((i - 1)/3)), ] !=
                  0) != 0) {
            j <- j + 1

            for (l in 1:block.ncol) {


              mat.cova.hsel[l, j] <- cov(block.array[, l, couples[i,
                                                                  1]], block.array[, l, couples[i, 2]], use = "pairwise.complete.obs")

              mat.cova.h[l, i] <- mat.cova.hsel[l, j]
            }

          }
        }



        #= Compute the temporal marginal covariance for each block =#

        nstaz <- length(selstaz)

        ntemp <- 0

        for (i in 1:(couples.ncol/2)) {

          if (sum(count_nozero_col[, i] != 0) != 0) {
            ntemp <- ntemp + 1
          }
        }
        mat.cova.u.ncol <- as.integer(couples.ncol/2)
        mat.cova.u <- matrix(0, block.ncol, mat.cova.u.ncol)
        mat.cova.ui <- matrix(0, block.ncol, nstaz)

        mat.cova.u.ncolsel <- ntemp
        mat.cova.usel <- matrix(0, block.ncol, mat.cova.u.ncolsel)
        mat.cova.uisel <- matrix(0, block.ncol, length(setdiff(selstaz,
                                                               selstaz[staz_zero])))

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

            if (length(staz_zero) == 0) {
              for (z in 1:nstaz) {
                for (l in 1:block.ncol) {

                  mat.cova.uisel[l, z] <- cov(block.array[-(nrow(matblock) -
                                                              couples[i, jj + 2] + 1:nrow(matblock)), l, z], block.array[-(1:couples[i,
                                                                                                                                     jj + 2]), l, z], use = "pairwise.complete.obs")

                }
              }
            } else {
              zz <- 0

              for (z in 1:nstaz) {
                if (length(intersect(selstaz[z], selstaz[staz_zero])) ==
                    0) {
                  zz <- zz + 1
                  for (l in 1:block.ncol) {
                    mat.cova.uisel[l, zz] <- cov(block.array[-(nrow(matblock) -
                                                                 couples[i, jj + 2] + 1:nrow(matblock)), l, z],
                                                 block.array[-(1:couples[i, jj + 2]), l, z], use = "pairwise.complete.obs")
                  }

                }
              }
            }



          }

          if (sum(count_nozero_col[, j] != 0) != 0) {
            jjj <- jjj + 1
            for (l in 1:block.ncol) {
              mat.cova.usel[l, jjj] <- mean(mat.cova.uisel[l, ], na.rm = TRUE)
              mat.cova.u[l, j] <- mean(mat.cova.uisel[l, ], na.rm = TRUE)
            }
          }

        }


      }

      #==========================================#
      #= Compute the matrix of covariances      =#
      #==========================================#

      if (type.test == 1 || type.test == 2) {
        mat.cova <- cbind(vec.cova00, mat.cova, mat.cova.h, mat.cova.u)
      }

      if (type.test == 3) {
        mat.cova <- cbind(mat.cova, mat.cova.h, mat.cova.u)
      }
      if (type.test >= 4) {



        #= Start permutation of each column of mat.cova:store sequentally the
        #= covariance for each set of comparisons (at least 3 spatial lags and 3
        #= temporal lags)

        mat.cova2 <- mat.cova
        mat.cova <- matrix(0, nrow = nrow(mat.cova), ncol = ncol(mat.cova))

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
                  for (zz in 1:((ncol(couples) - 2)/2/3)) {
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
                mat.cova[, z2:(z2 + count_nozero_row[ii, j + (jj - 1) *
                                                       3] - 1)] <- mat.cova2[, (count.skip + z3):(count.skip +
                                                                                                    z3 + count_nozero_row[ii, j + (jj - 1) * 3] - 1)]
                z2 <- z2 + count_nozero_row[ii, j + (jj - 1) * 3]
                z3 <- z3 + count_nozero_row[ii, j + (jj - 1) * 3] +
                  count.skip


              }


            }

          }
          z <- z2

          #= End permutation of each column of mat.cova:store sequentally the
          #= covariance for each set of comparisons (at least 3 spatial lags and 3
          #= temporal lags)


        }




        if (type.test == 4) {
          mat.cova <- cbind(mat.cova, mat.cova.hsel, mat.cova.usel)
        }
        if (type.test == 5) {
          mat.cova <- mat.cova
        }
        if (type.test == 6 || type.test == 7) {
          mat.cova <- cbind(mat.cova, mat.cova.usel)
        }

      }
      #= Write the matrix_cova after permutation


      #= end if on type of test
    }
    #= covariance matrix of covariances =#
    mat.cova.cova <- matrix(0, ncol(mat.cova), ncol(mat.cova))
    for (t in 1:ncol(mat.cova)) {
      mat.cova.cova[t, ] <- cov(mat.cova[, t], mat.cova, use = "pairwise.complete.obs")

    }



  }
  #==============================================#
  #= end test on type.test=0, 1, 2, 3, 4-7       =#
  #==============================================#


  #= Start defining the class covblocks =#

  if (type.test == 0) {

    #= mat.cova <- mat.cova
    mat.cova.h <- matrix(NA, 1, 1)
    mat.cova.u <- matrix(NA, 1, 1)
    #= mat.cova.cova <- mat.cova.cova
  }

  if (type.test == 1 || type.test == 2 || type.test == 3) {
    #= mat.cova = mat.cova mat.cova.h = mat.cova.h mat.cova.u = mat.cova.u
    #= mat.cova.cova = mat.cova.cova
  }

  if (type.test == 4) {

    #= mat.cova = mat.cova
    mat.cova.h <- mat.cova.hsel
    mat.cova.u <- mat.cova.usel
    #= mat.cova.cova = mat.cova.cova
  }
  if (type.test == 5) {

    #= mat.cova = mat.cova
    mat.cova.h <- matrix(NA, 1, 1)
    mat.cova.u <- matrix(NA, 1, 1)
    #= mat.cova.cova = mat.cova.cova
  }
  if (type.test == 6 || type.test == 7) {

    #= mat.cova = mat.cova,
    mat.cova.h <- matrix(NA, 1, 1)
    mat.cova.u <- mat.cova.usel
    #= mat.cova.cova = mat.cova.cova
  }



  new("covablocks", mat.cova = mat.cova, mat.cova.h = mat.cova.h, mat.cova.u = mat.cova.u,
      mat.cova.cova = mat.cova.cova, typetest = typetest)

  #= End defining the class covastat  =#

}
#' @include sepindex.R couples.R blocks.R
NULL
#' @param object object of class \code{covablocks} for method \code{show}
#' @rdname covablocks-class
#' @aliases covablocks-method
#' @export
setMethod(f="show", signature="covablocks", definition=function(object) {
  nb <- nrow(object@mat.cova)

  cat("An object of class 'covablocks', with", "\n")
  cat("number of blocks=", nb, "\n")
  cat("Slot 'mat.cova':")
  cat("\n")
  print(object@mat.cova)
  cat("\n")
  cat("Slot 'mat.cova.h':")
  cat("\n")
  # if(is.na(object@mat.cova.h) == TRUE){
  #   print("This slot is not available for the required typetest")
  # }else{
  #   print(object@mat.cova.h)
  # }
  if(object@typetest == "sym" || object@typetest == "intProduct" || object@typetest == "gneiting"){
    print("This slot is not available for the required typetest")
  }else{
    print(object@mat.cova.h)
  }
  cat("\n")
  cat("Slot 'mat.cova.u':")
  cat("\n")
  # if(!is.na(object@mat.cova.u)){
  #   print(object@mat.cova.u)
  # }else{
  #   print("This slot is not available for the required typetest")
  # }
  if(object@typetest == "sym" || object@typetest == "intProduct"){
    print("This slot is not available for the required typetest")
  }else{
    print(object@mat.cova.u)
  }
  cat("\n")
  cat("Slot 'mat.cova.cova':")
  cat("\n")
  print(object@mat.cova.cova)
  cat("\n")
  cat("Slot 'typetest':")
  cat("\n")
  print(object@typetest)
}

)
