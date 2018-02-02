#' Class "covaprop"
#'
#' A function for testing some properties (symmetry, separability, type of
#' non-separability) of spatio-temporal covariance functions and some classes
#' of space-time covariance models
#'
#' @slot test.statistics numeric; the value of the test statistic
#' @slot p.value numeric; the lower tail p-value of the test statistic
#' @slot df numeric; the degrees of freedom, if available
#' @slot typetest numeric; contains the code of the test to be performed
#'
#' @rdname covaprop-class
#' @exportClass covaprop
setClass("covaprop", slots = c(test.statistics = "matrix",
                               p.value = "matrix",
                               df = "ANY",
                               typetest = "numeric"))

#' @param cblock object of class \code{covablocks}
#' @param cstat object of class \code{covastat}
#' @param nonseptype integer; this argumet is required only for \code{typetest=2},
#' otherwise it has to be set equal to \code{NULL} (default choice).
#' Set \code{nonseptype=0} for testing the null hypothesis that the non
#' separability is non positive; set \code{nonseptype=1} for testing the null
#' hypothesis that the non separability is non negative
#' @param sign.level level of significance
#'
#' @note {
#' A message helps to decide for either reject the null hypothesis in favor of
#' the alternative or not reject it, at a specific level of significance
#' \itemize{
#' \item The test on full symmetry (\code{typetest=0}) represents the first step
#' for the selection of a suitable class of spatio-temporal covariance functions.
#' According to the definition of full symmetry, the null hypothesis to be
#' tested is \eqn{H_0: C(h,u) - C(h,-u)=0}
#'
#' \item The test of separability (\code{typetest=1}) represents the second
#' step of the testing procedure. According to the definition of separability,
#' the null hypothesis to be tested is
#' \eqn{H_0: C(h, u)/C(h, 0) - C(0, u)/C(0,0)=0}
#'
#' \item The test on the type of non separability (\code{typetest=2})
#' represents the third step for the selection of a suitable class of space-time
#' covariance functions.
#' According to the definition of type of non separability, the null hypothesis
#' to be tested is that the non separability is non negative
#' \eqn{H_0: C(h,u)/C(h,0) - C(0,u)/C(0,0) > 0}
#' or
#' \eqn{H_0: C(h,u)/C(h,0) - C(0,u)/C(0,0) < 0}
#' if the null hypothesis to test is that the non separability is non positive
#'
#' \item \code{typetest} from 3 to 5 are useful to test the goodness of a
#' specific class of space-time covariance function. For this testing procedure
#' the generic null hypothesis is: \eqn{H_0: Af(G)=0}
#' }
#' For the analytic expression of each test statistic and its probability
#' distribution see Cappello et al. (2017). In the same contribution the
#' different \code{f(G)} are given for each test to be computed.
#'
#' A stop occurs if the type of test set in \code{cblock} is not consistent with
#' the type of test set in \code{cstat}.
#' }
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
#' @seealso \linkS4class{couples}
#' @seealso \linkS4class{blocks}
#' @seealso \linkS4class{covablocks}
#' @seealso \linkS4class{covastat}
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
#' block.sym <- blocks(lb=40, ls=10, matdata = rr_13, pardata1 = 1, pardata2 = 1,
#' stpairs = couples.sym)
#'
#' covabl.sym <- covablocks(stblocks = block.sym, stpairs = couples.sym, typetest = 0)
#'
#' covast.sym <- covastat(matdata = rr_13, pardata1 = 1, pardata2 = 1,
#' stpairs = couples.sym, typetest = 0, beta.data = NULL)
#'
#' test.sym <- covaprop(cblock = covabl.sym, cstat = covast.sym, nonseptype = NULL,
#' sign.level = 0.05)
#'
#' ###method for covaprop
#' #1. show
#' test.sym
#'
#' @rdname covaprop-class
#' @export
covaprop <- function(cblock, cstat, nonseptype = NULL, sign.level = 0.05) {

  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) abs(x -
                                                                round(x)) < tol

  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}

  ### SOME CHECKS ON THE ARGUMENTS ###


  if (!inherits(cblock, "covablocks")){
    stop("cblock argument has to be of class covablocks")
  }

  if (!inherits(cstat, "covastat")){
    stop("cstat argument has to be of class covastat")
  }

 if(cblock@typetest != cstat@typetest){
    stop("The arguments cblocks and cstat are referred to different typetest. Please verify the consistency of the arguments")
  }

  typetest <- cblock@typetest

  if (is.null(nonseptype) == FALSE && typetest != 2) {
    message("The nonseptype will be ignored since it is not required for the
            selected type of test.")
  }
  if (typetest == 2) {
  if (is.scalar(nonseptype) == FALSE || nonseptype < 0 || nonseptype >
      1) {stop("The argument for nonseptype is not admissible.")}
  }

  if (is.numeric(sign.level) == FALSE) {
    stop("The argument for sign.level is not admissible.")
    }

  if (sign.level <= 0 || sign.level >= 1) {
    message("The specified level of significance is not admissible.
            It must be between 0 and 1")
    stop("Stop running")
  }


  if (typetest >= 3) {
    typetest <- typetest + 1
  }

  if (typetest == 0) {
    matrix_cova_cova <- cblock@mat.cova.cova
    G_vect <- cstat@G
    matrixA1 <- cstat@A
    df <- nrow(cstat@A)

    x <- tcrossprod(tcrossprod(matrixA1, matrix_cova_cova), matrixA1)

    test.statistics <- (t(matrixA1 %*% G_vect)) %*% (solve(x)) %*% (matrixA1 %*% G_vect)


    message("The test statistic is equal to")
    print(test.statistics)

    critical.value <- pchisq(test.statistics, df = nrow(matrixA1),
                             lower.tail = FALSE)

    message("The p-value is equal to")
    print(critical.value)

    if (critical.value >= sign.level) {
      message("Don't reject the null hypothesis of symmetry at ",
              sign.level, " level of significance")
    } else {
      message("Reject the null hypothesis of symmetry at ", sign.level,
              " level of significance")
    }
  }

  if (typetest == 1) {
    matrix_cova_cova <- cblock@mat.cova.cova
    f_G <- cstat@f.G
    matrixA1 <- cstat@A
    mat_B <- cstat@B
    df <- nrow(cstat@A)

    y <- crossprod(matrix_cova_cova, mat_B)
    x <- crossprod(mat_B, y)

    test.statistics <- (t(matrixA1 %*% f_G)) %*% (solve(matrixA1 %*% x %*% t(matrixA1))) %*% (matrixA1 %*% f_G)
    critical.value <- pchisq(test.statistics, df = nrow(matrixA1),
                             lower.tail = FALSE)


    message("The test statistic is equal to")
    print(test.statistics)

    message("The p-value is equal to")
    print(critical.value)

    if (critical.value >= sign.level) {
      message("Don't reject the null hypothesis of separability at ",
              sign.level, " level of significance")
    } else {
      message("Reject the null hypothesis of separability at ", sign.level,
              " level of significance")
    }
  }

  if (typetest == 2) {
    if (is.null(nonseptype) == TRUE) {
      message("The argument nonseptype must to be set")
      stop("Stop running")

    }




    matrix_cova_cova <- cblock@mat.cova.cova
    f_G <- cstat@f.G
    matrixA2 <- cstat@A
    mat_B <- cstat@B
    ### test_type is equal to 0 for testing the null hypothesis that the non
    ### separability is non positive ### test_type is equal to 1 for testing
    ### the null hypothesis that the non separability is non negative ###

    onevec <- matrix(1, nrow = nrow(matrixA2), ncol = 1)

    y <- crossprod(matrix_cova_cova, mat_B)
    x <- crossprod(mat_B, y)

    test.statistics <- (t(onevec) %*% (matrixA2 %*% f_G)) / sqrt(t(onevec) %*%
                       (matrixA2 %*% x %*% t(matrixA2)) %*% (onevec))

    message("The test statistic is equal to")
    print(test.statistics)

    if (nonseptype == 0) {
      critical.value <- pnorm(test.statistics, mean = 0, sd = 1,
                              lower.tail = FALSE)
      critical.value2 <- pnorm(test.statistics, mean = 0.5, sd = 1,
                               lower.tail = TRUE)
      message("The p-value is equal to")
      print(critical.value)
      if (critical.value >= sign.level) {
        message("Don't reject the null hypothesis of non positive non separability at
          ", sign.level, " level of significance")
      } else {
        message("Reject the null hypothesis of non positive non separability at
          ", sign.level, " level of significance")
      }
    }


    if (nonseptype == 1) {
      critical.value <- pnorm(test.statistics, mean = 0, sd = 1,
                              lower.tail = TRUE)
      critical.value2 <- pnorm(test.statistics, mean = -0.5, sd = 1,
                               lower.tail = FALSE)

      message("The p-value is equal to")
      print(critical.value)

      if (critical.value >= sign.level) {
        message("Don't reject the null hypothesis of non negative non separability at
          ", sign.level, " level of significance")
      } else {
        message("Reject the null hypothesis of non negative non separability at
          ", sign.level, " level of significance")
      }
    }

    df <- NA
  }

  if (typetest >= 4) {
    matrix_cova_cova <- cblock@mat.cova.cova
    f_G <- cstat@f.G
    matrixA.m <- cstat@A
    df <- nrow(cstat@A)
    mat_B <- cstat@B


    # the test statistic is equal to the one used for the separability test
    if (typetest == 4 || typetest == 5 || typetest == 7) {

      y <- crossprod(matrix_cova_cova, mat_B)
      x <- crossprod(mat_B, y)
      test.statistics <- (t(matrixA.m %*% f_G)) %*% (solve(matrixA.m %*%
                          x %*% t(matrixA.m))) %*% (matrixA.m %*% f_G)

      critical.value <- pchisq(test.statistics, df = nrow(matrixA.m),
                               lower.tail = FALSE)

      message("The test statistic is equal to")
      print(test.statistics)

      message("The p-value is equal to")
      print(critical.value)


      if (critical.value >= sign.level) {
        message("Don't reject the null hypothesis on the type of the model at ",
                sign.level, " level of significance")
      } else {
        message("Reject the null hypothesis on the type of the model at ",
                sign.level, " level of significance")
      }
    }


    if (typetest == 6) {

      nbeta <- dim(f_G)[3]
      test.statistics <- matrix(0, nrow = nbeta, ncol = 2)
      critical.value <- matrix(0, nrow = nbeta, ncol = 2)

      colnames(test.statistics) <- c("Test Statistics", "Beta")
      colnames(critical.value) <- c("p-value", "Beta")

      for (i in 1:nbeta) {
        y <- crossprod(matrix_cova_cova, mat_B[, , i])
        x <- crossprod(mat_B[, , i], y)
        test.sing <- matrixcalc::is.singular.matrix((matrixA.m %*% x %*%
                                          t(matrixA.m)), tol = 1e-08)

        test.statistics[i, 1] <- (t(matrixA.m %*% f_G[, , i])) %*%
          (solve(matrixA.m %*% x %*% t(matrixA.m))) %*% (matrixA.m %*%
                                                          f_G[, , i])
        test.statistics[i, 2] <- cstat@beta.data[i]

    critical.value[i, 1] <- pchisq(test.statistics[i, 1], df = nrow(matrixA.m),
                                       lower.tail = FALSE)
    critical.value[i, 2] <- cstat@beta.data[i]

        message("The test statistic is equal to")
        print(test.statistics)

        message("The p-value is equal to")
        print(critical.value)

        if (critical.value[i, 1] >= sign.level) {
        message("Don't reject the null hypothesis on the type of the model at ",
                  sign.level, " level of significance for beta", i)
        } else {
          message("Reject the null hypothesis on the type of the model at ",
                  sign.level, " level of significance for beta", i)
        }

      }
    }
  }

  if (typetest >= 4) {
    typetest <- typetest -1
  }

  new("covaprop", test.statistics = test.statistics, p.value = critical.value,
      df = df, typetest = typetest)

}

#' @include sepindex.R couples.R blocks.R covablocks.R covastat.R
NULL
#' @param object object of class \code{covaprop} for method \code{show}
#' @rdname covaprop-class
#' @aliases covaprop-class
#' @aliases covaprop-method
#' @aliases show
#' @export
setMethod(f="show", signature="covaprop", definition=function(object) {
  if(object@typetest == 0){
    test <- "symmetry test"
  }

  if(object@typetest == 1){
    test <- "separability test"
  }

  if(object@typetest == 2){
    test <- "type of non-separability test"
  }

  if(object@typetest == 3){
    test <- "test on the product-sum class of models"
  }

  if(object@typetest == 4){
    test <- "test on the integrated class of models"
  }

  if(object@typetest == 5){
    test <- "test on the Gneiting class of models"
  }

  cat("Results of the ",test, "\n")

  if(object@typetest == 0 || object@typetest == 1 || object@typetest == 3 || object@typetest == 4){
    cat("X-squared=", object@test.statistics, "\n")
    cat("df=", object@df, "\n")
    cat("p-value=", object@p.value, "\n")
  }

  if(object@typetest == 2){
    cat("Z=", object@test.statistics, "\n")
    cat("p-value=", object@p.value, "\n")
  }

  if(object@typetest == 5){
    cat("X-squared=", object@test.statistics[,1], "\n")
    cat("Beta=", object@test.statistics[,2], "\n")
    cat("df=", object@df, "\n")
    cat("p-value=", object@p.value, "\n")
  }


}
)
