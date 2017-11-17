#' Tests properties of S-T covariance functions
#'
#' A function for testing some properties (symmetry, separability, type of
#' non-separability) of spatio-temporal covariance functions and some classes
#' of space-time covariance models
#'
#' @param cblock object of class \code{covablocks}
#'
#' @param cstat object of class \code{covastat}
#'
#' @param typetest integer; set \code{typetest=0} for symmetry test (default
#' choice), \code{typetest=1} for separability test, \code{typetest=2} for type
#' of non separability test, \code{typetest=3} for the test on the product-sum
#' class of models, \code{typetest=4} for the test on the integrated product
#' class of models, \code{typetest=5} for the test on the Gneiting class of
#' models
#'
#' @param nonseptype integer; this argumet is required only for \code{typetest}=2,
#' otherwise it has to be set equal to \code{NULL} (default choice).
#' Set \code{nonseptype=0} for testing the null hypothesis that the non
#' separability is non positive; set \code{nonseptype=1} for testing the null
#' hypothesis that the non separability is non negative
#'
#' @param sign.level level of significance
#'
#' @return The value of the test statistic and the lower tail p-value of the
#' test statistic. Moreover a message helps to decide for either reject the null
#' hypothesis in favor of the alternative or not reject it, at a specific level
#' of significance
#'
#' @note {
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
#' distribution see references. In the same contribution the
#' different \code{f(G)} are given for each test to be computed.
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
#' Cappello, C., De Iaco, S., Posa, D., 2016, Testing the type of
#' non-separability and some classes of covariance models for space-time data.
#' Stochastic Environmental Research and Risk Assessment,
#' doi 10.1007/s00477-017-1472-2
#'
#' @examples
#' # Before running this function, it is necessary to execute couples, blocks,
#' # covablocks and covastat, as specified in the examples of the corresponding
#' # help pages.
#' #
#' # To run the example, paste and copy the following lines
#' # (without the symbol '#') in the console
#' #
#' # coupl_sim <- couples(typetest = 0, typecode = character())
#' # blocks_sim <- blocks(lb = 40, ls = 10, matdata = rr_13, stpairs = coupl_sim)
#' # covabl_sim <- covablocks(stblocks = blocks_sim, stpairs = coupl_sim, typetest = 0)
#' # covast_sim <- covastat(matdata = rr_13, stpairs = coupl_sim, typetest = 0)
#' # covaprop(cblock = covabl_sim, cstat = covast_sim, typetest = 0, nonseptype = NULL,
#' # sign.level = 0.05)
#'
#' @seealso \linkS4class{couples}
#' @seealso \linkS4class{blocks}
#' @seealso \linkS4class{covablocks}
#' @seealso \linkS4class{covastat}
#'
#' @export
covaprop <- function(cblock, cstat, typetest, nonseptype = NULL, sign.level = 0.05) {

  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) abs(x -
                                                                round(x)) < tol

  if (is.wholenumber(typetest) == FALSE || typetest < 0 || typetest > 5) {
    stop("The argument for typetest is not admissible.")
  }

  if (is.null(nonseptype) == FALSE && typetest != 2) {
    message("The nonseptype will be ignored since it is not required for the
            selected type of test.")
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


    test.statistics <- (t(matrixA1 %*% G_vect)) %*% (solve(matrixA1 %*%
                                                             matrix_cova_cova %*% t(matrixA1))) %*% (matrixA1 %*% G_vect)

    message("The test statistic is equal to")
    print(test.statistics)

    critical.value <- pchisq(test.statistics, df = nrow(matrixA1),
                             lower.tail = FALSE)

    message("The critical value is equal to")
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

    test.statistics <- (t(matrixA1 %*% f_G)) %*% (solve(matrixA1 %*%
                                                          t(mat_B) %*% matrix_cova_cova %*% mat_B %*% t(matrixA1))) %*%
      (matrixA1 %*% f_G)

    critical.value <- pchisq(test.statistics, df = nrow(matrixA1),
                             lower.tail = FALSE)


    message("The test statistic is equal to")
    print(test.statistics)

    message("The critical value is equal to")
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
    if (is.wholenumber(nonseptype) == FALSE || nonseptype < 0 || nonseptype >
        1) {
      stop("The argument for nonseptype is not admissible.")
    }

    matrix_cova_cova <- cblock@mat.cova.cova
    f_G <- cstat@f.G
    matrixA2 <- cstat@A
    mat_B <- cstat@B
    ### test_type is equal to 0 for testing the null hypothesis that the non
    ### separability is non positive ### test_type is equal to 1 for testing
    ### the null hypothesis that the non separability is non negative ###

    onevec <- matrix(1, nrow = nrow(matrixA2), ncol = 1)

    test.statistics <- (t(onevec) %*% (matrixA2 %*% f_G)) / sqrt(t(onevec) %*%
                       (matrixA2 %*% t(mat_B) %*% matrix_cova_cova %*% mat_B %*%
                           t(matrixA2)) %*% (onevec))

    message("The test statistic is equal to")
    print(test.statistics)

    if (nonseptype == 0) {
      critical.value <- pnorm(test.statistics, mean = 0, sd = 1,
                              lower.tail = FALSE)
      critical.value2 <- pnorm(test.statistics, mean = 0.5, sd = 1,
                               lower.tail = TRUE)
      message("The critical value is equal to")
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

      message("The critical value is equal to")
      print(critical.value)

      if (critical.value >= sign.level) {
        message("Don't reject the null hypothesis of non negative non separability at
          ", sign.level, " level of significance")
      } else {
        message("Reject the null hypothesis of non negative non separability at
          ", sign.level, " level of significance")
      }
    }


    #= print(critical.value2)
  }

  if (typetest >= 4) {
    matrix_cova_cova <- cblock@mat.cova.cova
    f_G <- cstat@f.G
    matrixA.m <- cstat@A
    mat_B <- cstat@B
    # the test statistic is equal to the one used for the separability test
    if (typetest == 4 || typetest == 5 || typetest == 7) {
      test.statistics <- (t(matrixA.m %*% f_G)) %*% (solve(matrixA.m %*%
                          t(mat_B) %*% matrix_cova_cova %*% mat_B %*%
                            t(matrixA.m))) %*% (matrixA.m %*% f_G)

      critical.value <- pchisq(test.statistics, df = nrow(matrixA.m),
                               lower.tail = FALSE)

      message("The test statistic is equal to")
      print(test.statistics)

      message("The critical value is equal to")
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
      test.statistics <- matrix(0, nrow = nbeta, ncol = 1)
      critical.value <- matrix(0, nrow = nbeta, ncol = 1)
      for (i in 1:nbeta) {
        test.sing <- matrixcalc::is.singular.matrix((matrixA.m %*% t(mat_B[, , i]) %*%
                                        matrix_cova_cova %*% mat_B[, , i] %*%
                                          t(matrixA.m)), tol = 1e-08)

        test.statistics[i, 1] <- (t(matrixA.m %*% f_G[, , i])) %*%
          (solve(matrixA.m %*% t(mat_B[, , i]) %*% matrix_cova_cova %*%
                   mat_B[, , i] %*% t(matrixA.m))) %*% (matrixA.m %*%
                                                          f_G[, , i])

    critical.value[i, 1] <- pchisq(test.statistics[i, 1], df = nrow(matrixA.m),
                                       lower.tail = FALSE)

        message("The test statistic is equal to")
        print(test.statistics)

        message("The critical value is equal to")
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

}

#' @include sepindex.R couples.R blocks.R covablocks.R covastat.R
NULL
