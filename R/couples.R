#' Class "couples"
#'
#' A class for spatial points and the corresponding temporal lags to be
#' analyzed in order to test some covariance properties and some well known
#' classes of space-time covariance functions models
#'
#' @slot couples.st matrix; in which the first two columns contain the
#' couples of spatial points (denoted with order numbers) to be analyzed
#' and the other columns the temporal lags associated with each spatial couples
#' @slot sel.staz numeric or character; contains the ID codes of the
#' selected spatial points
#' @slot sp.couples data.frame; contains the couples of order numbers
#' associated with the spatial points to be analyzed and the couples of the
#' ID codes
#' @slot tl.couples numeric; contains the temporal lags associated to the
#' couples of the selected spatial points
#'
#' @rdname couples-class
#' @exportClass couples
#' @export
setClass("couples", slots = c(couples.st = "matrix",
                                         sel.staz = "ANY",
                                         sp.couples = "data.frame",
                                         tl.couples = "numeric"))

#' @param typetest integer; set typetest=0 for simmetry test (default choice),
#' \code{typetest=1} for separability test, \code{typetest=2} for type of non
#' separability test,
#' \code{typetest=3} for the test on the product-sum class of models,
#' \code{typetest=4} for the test on the integrated product class of models,
#' \code{typetest=5} for the test on the Gneiting class of models
#'
#' @param typecode numeric or character; specifies the codification of the
#' spatial points in the \code{data frame} or in the STFDF/STSDF, typically
#' output from \code{\link{dataprep}}
#'
#' @details The function requires the user to set some external arguments:
#' \enumerate{
#' \item the number of spatial points to be analyzed
#' \item the sequence of ID codes which denote the spatial points to be analyzed
#' \item the number of spatial couples to be compared
#' \item the couples of different spatial points
#' \item the couples of positive and negative temporal lags
#' }
#' If some temporal lags, corresponding to some couples of spatial
#' points, are not required for the specific test, might be set equal to zero.
#'
#' @note It is important to note that for:
#' \itemize{
#' \item symmetry test (\code{typetest=0}), both positive and negative
#' temporal lags have to be considered;
#'
#' \item separability and type of non separability tests (\code{typetest=1} and
#' \code{typetest=2}, respectively), both positive and negative
#' temporal lags have to be considered. Hovewer, if the simmetry hyphotesis
#' has not been rejected, only positive temporal lags should be set.
#' Moreover for \code{typetest=2} the temporal lags should be chosen according
#' to the results of the sample non separability ratios, plotted through
#' a boxplot classified for temporal lags (see \linkS4class{sepindex} for more
#' details);
#'
#' \item model tests (\code{typetest} from 3 to 5), the number of analyzed spatial
#' points must be used to create at least 3 spatial couples or multiple of 3,
#' such that each triplet satisfies the condition h1-h2=h2-h3. The number
#' of positive temporal lags must be at least 3, or multiple of 3, too. The
#' condition u1-u2=u2-u3 must be satisfied for each triplet.
#'
#' }
#'
#'
#' Moreover, errors occur if
#' \itemize{
#' \item some spatial points, given in the sequence at the beginning of the
#' function, have not been used to generate the couples of spatial points
#'
#' \item there is at least one spatial couple with no specification of
#' temporal lags
#'
#' \item no temporal lags have been specified
#'
#' \item the number of spatial points fixed in \code{stpairs} (object of class
#' \code{couples}) is less than 2
#' }
#'
#' @examples
#' ## The function requires to set some external arguments.
#' # In the example regarding the simmetry test (typetest = 0),
#' # 12 spatial points, with ID codes: DERP016, ..., DESN049, have been selected,
#' # then 6 spatial couples ([DERP016, DENW065], .., [DEBY047, DESN049])
#' # have been formed for comparison. Moreover, 4 positive and negative temporal
#' # lags have been considered (i.e. +1, -1, +2, -2). Finally, no (N) temporal
#' # lag has been set equal to zero. Hence, for the simmetry test 24 (6*4)
#' # spatio-temporal comparisons have been fixed.
#' #
#' # To run the example, paste and copy the following lines
#' # (without the symbol '#') in the console
#' #
#' #coupl_sim <- couples(typetest = 0, typecode = character())
#' #12
#' #DERP016
#' #DENW065
#' #DENW063
#' #DENI019
#' #DENW068
#' #DEHE046
#' #DEHE051
#' #DETH026
#' #DEUB029
#' #DETH061
#' #DEBY047
#' #DESN049
#' #6
#' #DERP016
#' #DENW065
#' #DENW063
#' #DENI019
#' #DENW068
#' #DEHE046
#' #DEHE051
#' #DETH026
#' #DEUB029
#' #DETH061
#' #DEBY047
#' #DESN049
#' #4
#' #1
#' #-1
#' #2
#' #-2
#' #N
#' @rdname couples-class
#' @export
couples <- function(typetest = 0, typecode = numeric()) {


  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) abs(x -
                                                                round(x)) < tol


  if (is.wholenumber(typetest) == FALSE || typetest < 0 || typetest >
      5) {
    stop("The argument for typetest is not admissible.")
  }

  if (typetest >= 3) {
    typetest <- typetest + 1
  }

  if (typeof(typecode) != "double" && typeof(typecode) != "character") {
    stop("The argument for typecode is not admissible.")
  }


  ### SPATIAL COUPLES ###

  ns <- readline(prompt = "Enter the number (at least 2) of spatial points to be analyzed: ")
  ns <- as.integer(ns)

  if (ns < 2) {
    message("The number of spatial points to be analyzed must be at least equal to two")
    stop("Stop running")
  }


  if (typetest >= 4) {
    message("Note that the number of spatial points must be used to create at least 3 spatial couples ")
    message("or multiple of 3 spatial couples")


  }

  message("Enter a sequence of ID codes which denote the spatial points to be analyzed: ")
  sel.staz <- scan(file = "", n = ns, what = typecode)
  j <- 0

  while (anyDuplicated(sel.staz) >= 1) {
    j <- j + 1
    if (as.integer(j / 5) == (j / 5)) {
      message("The last sequences of points were not admissible")
      ans_YN <- readline(prompt = "Would you like to continue? (Y/N)")
      if (ans_YN == "N" || ans_YN == "n") {
        stop("Stop running")
      }

    }
    message("The spatial points in the sequence must be different ")
    message("Enter a new sequence of spatial points to be analyzed: ")
    sel.staz <- scan(file = "", n = ns, what = typecode)
  }


  message("********************************************************************************")
  message("*** One of the slots of couples class has been created.                      ***")
  message("*** This is called @sel.staz and contains the spatial points to be analyzed. ***")
  message("*** This first output is an argument for further functions.                  ***")
  message("********************************************************************************")

  nc <- readline(prompt = "Enter the number of spatial couples to be compared: ")
  nc <- as.integer(nc)

  if (typetest >= 4) {

    ns_multiple <- as.integer(nc / 3)
    if ((ns_multiple * 3) != nc) {
      message("The total number of spatial couples is not consistent with the number of couples required by the test")
      stop("Stop running")
    }
  }

  for (i in 1:nc) {
    message("Enter the couple (# ", i, ") of different spatial points")

    s.couples <- scan(file = "", n = 2, what = typecode)


    check.na <- match(s.couples, sel.staz)
    while (s.couples[2] == s.couples[1] || anyNA(check.na) == TRUE) {
      message("This couple of points is not admissible: the points must be different")
      message("and must be chosen among the spatial points previously selected")
      message("Enter again the couple (# ", i, ") of different spatial points")

      s.couples <- scan(file = "", n = 2, what = typecode)
      check.na <- match(s.couples, sel.staz)
    }



    if (i == 1) {
      sp.couples <- matrix(s.couples, nrow = 1, ncol = 2, byrow = TRUE)
    }
    if (i > 1) {
      s.couples.perm <- c(s.couples[2], s.couples[1])

      sp.couples2 <- matrix(s.couples, nrow = 1, ncol = 2, byrow = TRUE)
      sp.couples <- rbind(sp.couples, sp.couples2)
      dupli.couple <- duplicated(sp.couples, fromLast = T)

      sp.couples2.perm <- matrix(s.couples.perm, nrow = 1, ncol = 2,
                                 byrow = TRUE)
      sp.couples.perm <- rbind(sp.couples, sp.couples2.perm)
      dupli.couple.perm <- duplicated(sp.couples.perm, fromLast = T)

      j <- 0
      while (dupli.couple[1] == TRUE || dupli.couple.perm[1] == TRUE) {

        if (dupli.couple[1] == TRUE) {
          j <- j + 1

          sp.couples <- sp.couples[-i, ]
          if (as.integer(j / 5) == (j / 5)) {
            message("The last couples of points were not admissible")
            ans_YN <- readline(prompt = "Would you like to print the couples indicated up to now? (Y/N)")
            if (ans_YN == "Y" || ans_YN == "y") {
              print(sp.couples)
            }
            ans_YN <- readline(prompt = "Would you like to continue? (Y/N)")
            if (ans_YN == "N" || ans_YN == "n") {
              stop("Stop running")
            }

          }


          message("The previous couple of points already exists")
          message("Enter another couple (# ", i, ") of different spatial points")

          s.couples <- scan(file = "", n = 2, what = typecode)

          check.na <- match(s.couples, sel.staz)
          while (s.couples[2] == s.couples[1] || anyNA(check.na) ==
                 TRUE) {
            message("This couple of points is not admissible: the points must be different")
            message("and must be chosen among the spatial points previously selected")
            message("Enter again the couple (# ", i, ") of different spatial points")

            s.couples <- scan(file = "", n = 2, what = typecode)
            check.na <- match(s.couples, sel.staz)
          }

          s.couples.perm <- c(s.couples[2], s.couples[1])

          sp.couples2 <- matrix(s.couples, nrow = 1, ncol = 2,
                                byrow = TRUE)
          sp.couples <- rbind(sp.couples, sp.couples2)
          dupli.couple <- duplicated(sp.couples, fromLast = T)

          sp.couples2.perm <- matrix(s.couples.perm, nrow = 1,
                                     ncol = 2, byrow = TRUE)
          sp.couples.perm <- rbind(sp.couples, sp.couples2.perm)
          dupli.couple.perm <- duplicated(sp.couples.perm, fromLast = T)
        }

        if (dupli.couple.perm[1] == TRUE) {
          message("This couple of points, indicated in different order, already exists")
          ans_YN <- readline(prompt = "Would you like to change this couple of points? (Y/N) ")
          dupli.couple.perm[1] <- FALSE
          if (ans_YN == "Y" || ans_YN == "y") {

            sp.couples <- sp.couples[-i, ]

            message("Enter another couple (# ", i, ") of different spatial points")

            s.couples <- scan(file = "", n = 2, what = typecode)

            check.na <- match(s.couples, sel.staz)
            while (s.couples[2] == s.couples[1] || anyNA(check.na) ==
                   TRUE) {
              message("This couple of points is not admissible: the points must be different")
              message("and must be chosen among the spatial points previously selected")
              message("Enter again the couple (# ", i, ") of different spatial points")

              s.couples <- scan(file = "", n = 2, what = typecode)
              check.na <- match(s.couples, sel.staz)
            }

            s.couples.perm <- c(s.couples[2], s.couples[1])

            sp.couples2 <- matrix(s.couples, nrow = 1, ncol = 2,
                                  byrow = TRUE)
            sp.couples <- rbind(sp.couples, sp.couples2)
            dupli.couple <- duplicated(sp.couples, fromLast = T)

            sp.couples2.perm <- matrix(s.couples.perm, nrow = 1,
                                       ncol = 2, byrow = TRUE)
            sp.couples.perm <- rbind(sp.couples, sp.couples2.perm)
            dupli.couple.perm <- duplicated(sp.couples.perm, fromLast = T)

          }
        }


      }  # endwhile on duplicates


    }


    if (typetest >= 4 && (as.integer(i / 3)) == (i / 3)) {
      k <- (as.integer(i / 3))
      message("The  triplet (# ", k, ") of different spatial points will be used for the spatial comparison")
      message("The condition h1-h2=h2-h3 must be satisfied for each triplet")
      ans_YN <- readline(prompt = "Does this spatial triplet satisfy this condition? (Y/N)")
      if (ans_YN == "N" || ans_YN == "n") {
        stop("Stop running")
      }
    }


  }
  if (length(setdiff(sel.staz, sp.couples)) >= 1) {
    message("Error: the following spatial points given in the sequence at the beginning, have not been used to generate the couples of spatial points")
    print(setdiff(sel.staz, sp.couples))
    stop("Stop running")
  }
  sp.couples2 <- sp.couples
  sp.couples <- matrix(match(sp.couples, sel.staz), nrow = nc, ncol = 2)



  ### TEMPORAL COUPLES ###

  #===================================================================#
  #= Start defining temporal couples for                             =#
  #= test of symmetry (typetest=0), separability (typetest=1) and    =#
  #= type of non separability (typetest=2)                           =#
  #===================================================================#

  if (typetest == 1 || typetest == 0 || typetest == 2) {
    n.temp <- readline(prompt = "Enter the maximum number of positive and negative temporal lags to be considered: ")
    n.temp <- as.integer(n.temp)

    while (as.integer(n.temp / 2) != (n.temp / 2)) {
      message("The maximum number of temporal lags must be an even number")
      message("The number of positive temporal lags is equal to the number of negative temporal lags")
      message("Before the end of this process, temporal lags not used can be set equal to 0 ")

      n.temp <- readline(prompt = "Enter again the maximum number of positive and negative temporal lags to be considered: ")
      n.temp <- as.integer(n.temp)
    }
    for (i in 1:(n.temp / 2)) {
      message("Enter the couple (# ", i, ") of positive and negative temporal lags (press enter after each input)")

      t.couples <- scan(file = "", n = 2)

      while ((t.couples[2] + t.couples[1] != 0) || (t.couples[1] == 0)) {
        message("This couple of temporal lags is not admissible: the temporal lags must be equal in absolut value")
        message(" and their absolut values must be greater than zero")
        message("Enter again the couple (# ", i, ") of positive and negative temporal lags")

        t.couples <- scan(file = "", n = 2)
      }


      if (i == 1) {
        tl.couples <- matrix(t.couples, nrow = 1, ncol = 2, byrow = TRUE)
      }
      if (i > 1) {
        t.couples.perm <- c(t.couples[2], t.couples[1])

        tl.couples2 <- matrix(t.couples, nrow = 1, ncol = 2, byrow = TRUE)
        tl.couples <- rbind(tl.couples, tl.couples2)
        dupli.couple <- duplicated(tl.couples, fromLast = T)

        tl.couples2.perm <- matrix(t.couples.perm, nrow = 1, ncol = 2,
                                   byrow = TRUE)
        tl.couples.perm <- rbind(tl.couples, tl.couples2.perm)
        dupli.couple.perm <- duplicated(tl.couples.perm, fromLast = T)

        j <- 0
        while (dupli.couple[1] == TRUE || dupli.couple.perm[1] ==
               TRUE) {

          j <- j + 1

          tl.couples <- tl.couples[-i, ]
          if (as.integer(j / 5) == (j / 5)) {
            message("The last couples of temporal lags were not admissible")
            ans_YN <- readline(prompt = "Would you like to print the couples of temporal lags indicated up to now? (Y/N)")
            if (ans_YN == "Y" || ans_YN == "y") {
              print(tl.couples)
            }
            ans_YN <- readline(prompt = "Would you like to continue? (Y/N)")
            if (ans_YN == "N" || ans_YN == "n") {
              stop("Stop running")
            }

          }


          message("The previous couple of temporal lags already exists or have been included in a different order")
          message("Enter another couple (# ", i, ") of positive and negative temporal lags")

          t.couples <- scan(file = "", n = 2)

          while ((t.couples[2] + t.couples[1] != 0) || (t.couples[1] == 0)) {

            message("This couple of temporal lags is not admissible: the temporal lags must be equal in absolut value")
            message("Enter again the couple (# ", i, ") of positive and negative temporal lags")

            t.couples <- scan(file = "", n = 2)
          }

          t.couples.perm <- c(t.couples[2], t.couples[1])

          tl.couples2 <- matrix(t.couples, nrow = 1, ncol = 2,
                                byrow = TRUE)
          tl.couples <- rbind(tl.couples, tl.couples2)
          dupli.couple <- duplicated(tl.couples, fromLast = T)

          tl.couples2.perm <- matrix(s.couples.perm, nrow = 1,
                                     ncol = 2, byrow = TRUE)
          tl.couples.perm <- rbind(tl.couples, tl.couples2.perm)
          dupli.couple.perm <- duplicated(tl.couples.perm, fromLast = T)


        }

      }
    }


    for (i in 1:(n.temp / 2)) {
      if (i == 1) {
        tl.couples.mat <- tl.couples[i, ]
      }
      if (i > 1) {
        tl.couples.mat2 <- tl.couples[i, ]
        tl.couples.mat <- cbind(tl.couples.mat, tl.couples.mat2)
      }
    }

    temp.lags <- matrix(tl.couples.mat, nrow = nrow(sp.couples), ncol = n.temp,
                        byrow = TRUE)

    couples.st <- cbind(sp.couples, temp.lags)
    message("This is a preview of the couples of the spatial points and the temporal lags to be analyzed.")
    print(cbind(sp.couples2, temp.lags))

    message("If some temporal lags, corresponding to some couples of spatial points, are not required,")
    message("they can be set equal to 0. In this case, please answer YES to the following question")
    edit_YN <- readline(prompt = "Would you like to set some temporal lags equal to 0? (Y/N)")
    if (edit_YN == "Y" || edit_YN == "y") {
      couples.st2 <- edit(couples.st)


      ### Check on nonadmissible changes on the spatial points

      iflag <- 0
      for (i in 1:nrow(sp.couples)) {
        if ((couples.st[i, 1]) != (couples.st2[i, 1]) || (couples.st[i,
                                                   2]) != (couples.st2[i, 2])) {
          if (iflag == 0) {
            message("The editing of spatial points is not admissible, the initial input has been restored")
            iflag <- 1
          }
          couples.st2[i, 1] <- couples.st[i, 1]
          couples.st2[i, 2] <- couples.st[i, 2]
        }

        #== Check on nonadmissible changes on the temporparal lags (typetest=1
        #== and 2)


        if (typetest == 1 || typetest == 2) {
          for (k in 1:n.temp) {

            if ((couples.st[i, k + 2]) != (couples.st2[i, k + 2]) &&
                (couples.st2[i, k + 2]) != 0) {
              message("Regarding the spatial couple (", sp.couples2[i, 1], " ; ", sp.couples2[i, 2], " the temporal lag ", couples.st2[i, jj + 2], " is not admissible")
              message("You can replace it with the previous one i.e. ", couples.st[i, jj + 2])
              message("Otherwise, you can substitute 0, just to esclude this temporal lag for this specific spatial couple")
              edit.couple_YN <- readline(prompt = "Do you want to a) replace it or b) substitute it with 0? (a/b)")

              j <- 0

              while (edit.couple_YN != "b" && edit.couple_YN !=
                     "B" && edit.couple_YN != "a" && edit.couple_YN !=
                     "A") {
                j <- j + 1
                if (j == 5)
                  stop("The number of attempts has been exceeded")
                edit.couple_YN <- readline(prompt = "Please press a or b ")
              }


              if (edit.couple_YN == "b" || edit.couple_YN == "B") {
                couples.st2[i, jj + 2] <- 0
              } else {
                couples.st2[i, jj + 2] <- couples.st[i, jj + 2]
              }

            }
          }
        }

        #== Check on nonadmissible changes on the temporparal lags (typetest=0)

        if (typetest == 0) {
          jj <- -1
          for (j in 1:(n.temp / 2)) {
            jj <- jj + 2
            if (((couples.st[i, jj + 2] != couples.st2[i, jj + 2]) ||
                 (couples.st[i, jj + 3] != couples.st2[i, jj + 3])) &&
                 (couples.st2[i, jj + 3] != 0 ||
                  couples.st2[i, jj + 2] != 0)) {

              message("Regarding the spatial couple (", sp.couples2[i, 1], " ; ", sp.couples2[i, 2], ") the temporal lags (", couples.st2[i, jj + 2], " ; ", couples.st2[i, jj + 3], ") are not admissible")
              message("You can replace them with the previous ones i.e. (", couples.st[i, jj + 2], " ; ", couples.st[i, jj + 3], ")")
              message("Otherwise, you can substitute (0 ; 0), just to esclude these temporal lags for this specific spatial couple")
              edit.couple_YN <- readline(prompt = "Do you want to a) replace them with the previous one or b) substitute them with (0 ; 0)? (a/b)")

              j <- 0

              while (edit.couple_YN != "b" && edit.couple_YN !=
                     "B" && edit.couple_YN != "a" && edit.couple_YN !=
                     "A") {
                j <- j + 1
                if (j == 5)
                  stop("The number of attempts has been exceeded")
                edit.couple_YN <- readline(prompt = "Please press a or b ")
              }

              if (edit.couple_YN == "b" || edit.couple_YN == "B") {
                couples.st2[i, jj + 2] <- 0
                couples.st2[i, jj + 3] <- 0
              } else {
                couples.st2[i, jj + 2] <- couples.st[i, jj + 2]
                couples.st2[i, jj + 3] <- couples.st[i, jj + 3]
              }
            }
          }
        }

        if (identical(as.integer(matrix(0, nrow = 1, ncol = (ncol(couples.st) -
                                                             2))), as.integer(couples.st2[i, 3:ncol(couples.st)])) ==
            TRUE) {
          stop("There is at least one spatial couple with no specification of temporal lags.")
        }
      }
      couples.st <- couples.st2

    }


  }
  #==============================================================#
  #= end test on type=0 (simmetry), type=1 (separability) and   =#
  #= type=2(type of non separability)                           =#
  #==============================================================#


  #=================================================#
  #== Start defining temporal couples for          =#
  #== test on the type of model (typetest>=4)      =#
  #=================================================#

  if (typetest >= 4) {
    n.temp <- readline(prompt = "Enter the maximum number of positive temporal lags to be considered: ")
    n.temp <- as.integer(n.temp) * 2
    t.couples <- matrix(0, nrow = 2, ncol = 1)

    message("Note that the number of positive temporal lags must be at least 3")
    message("or multiple of 3 ")

    nt_multiple <- as.integer(n.temp / 3)
    if ((nt_multiple * 3) != n.temp) {
      message("The total number of temporal lags is not consistent with the number of lags required by the test.")
      stop("Stop running")
    }



    for (i in 1:(n.temp / 2)) {
      message("Enter the positive temporal lag (# ", i, ") ")


      t.couples[1] <- scan(file = "", n = 1)
      t.couples[2] <- 0


      while (t.couples[1] <= 0) {
        message("This temporal lag is not admissible: the temporal lags must be greater than zero")
        message("Enter again the temporal lag (# ", i, ") ")

        t.couples[1] <- scan(file = "", n = 1)
        t.couples[2] <- 0
      }




      if (i == 1) {
        tl.couples <- matrix(t.couples, nrow = 1, ncol = 2, byrow = TRUE)
      }
      if (i > 1) {

        tl.couples2 <- matrix(t.couples, nrow = 1, ncol = 2, byrow = TRUE)
        tl.couples <- rbind(tl.couples, tl.couples2)
        dupli.couple <- duplicated(tl.couples, fromLast = T)


        j <- 0
        while (dupli.couple[1] == TRUE) {

          j <- j + 1

          tl.couples <- tl.couples[-i, ]
          if (as.integer(j / 5) == (j / 5)) {
            message("The last temporal lags were not admissible")
            ans_YN <- readline(prompt = "Would you like to print the temporal lags indicated up to now? (Y/N)")
            if (ans_YN == "Y" || ans_YN == "y") {
              print(tl.couples[1, ])
            }
            ans_YN <- readline(prompt = "Would you like to continue? (Y/N)")
            if (ans_YN == "N" || ans_YN == "n") {
              stop("Stop running")
            }

          }


          message("The previous temporal lag already exists")
          message("Enter another positive temporal lag (# ", i,
                  ") ")

          t.couples[1] <- scan(file = "", n = 1)
          t.couples[2] <- 0

          while (t.couples[1] <= 0) {
            message("This temporal lag is not admissible: the temporal lags must be greater than zero")
            message("Enter again the temporal lag (# ", i, ") ")

            t.couples[1] <- scan(file = "", n = 1)
            t.couples[2] <- 0
          }

          tl.couples2 <- matrix(t.couples, nrow = 1, ncol = 2, byrow = TRUE)
          tl.couples <- rbind(tl.couples, tl.couples2)
          dupli.couple <- duplicated(tl.couples, fromLast = T)



        }

      }


      if ((as.integer(i / 3)) == (i / 3)) {
        k <- (as.integer(i / 3))
        message("The  triplet (# ", k, ") of different temporal lags will be used for the temporal comparison")
        message("The condition u1-u2=u2-u3 must be satisfied for each triplet")
        ans_YN <- readline(prompt = "Does the temporal triplet satisfy this condition? (N to stop; any key to continue)")
        if (ans_YN == "N" || ans_YN == "n") {
          stop("Stop running")
        }

      }

    }


    for (i in 1:(n.temp / 2)) {
      if (i == 1) {
        tl.couples.mat <- tl.couples[i, ]
      }
      if (i > 1) {
        tl.couples.mat2 <- tl.couples[i, ]
        tl.couples.mat <- cbind(tl.couples.mat, tl.couples.mat2)
      }
    }

    temp.lags <- matrix(tl.couples.mat, nrow = nrow(sp.couples), ncol = n.temp,
                        byrow = TRUE)

    couples.st <- cbind(sp.couples, temp.lags)
    message("This is a preview of the couples of the spatial points and the temporal lags to be analyzed.")
    print(cbind(sp.couples2, temp.lags))

    message("If some temporal lags, corresponding to some couples of spatial points, are not required,")
    message("they can be set equal to 0. Note that al least one row and one column must contain positive values.")
    message("For changes in the temporal lags setting, please answer YES to the following question")
    edit_YN <- readline(prompt = "Would you like to set some temporal lags equal to 0? (Y/N)")
    couples.st2 <- couples.st
    if (edit_YN == "Y" || edit_YN == "y") {
      couples.st2 <- edit(couples.st)


      #= Check on nonadmissible changes on the spatial points

      iflag <- 0
      for (i in 1:nrow(sp.couples)) {
        if ((couples.st[i, 1]) != (couples.st2[i, 1]) || (couples.st[i,
                                                                     2]) != (couples.st2[i, 2])) {
          if (iflag == 0) {
            message("The editing of spatial points is not admissible, the initial input has been restored")
            iflag <- 1
          }
          couples.st2[i, 1] <- couples.st[i, 1]
          couples.st2[i, 2] <- couples.st[i, 2]
        }

        #= Check on nonadmissible changes on the temporparal lags (typetest>=4)


        for (k in 1:n.temp) {

          if ((couples.st[i, k + 2]) != (couples.st2[i, k + 2]) &&
              (couples.st2[i, k + 2]) != 0) {

            message("Regarding the spatial couple (", sel.staz[couples[i,
                                                                       1]], "", sel.staz[couples[i, 2]], " the temporal lag ",
                    couples.st2[i, k + 2], " is not admissible")
            message("You can replace it with the previous one i.e. ",
                    couples.st[i, k + 2])
            message("Otherwise, you can substitute 0, just to esclude this temporal lag for this specific spatial couple")
            message("Do you want to a) replace it with ", couples.st[i, k + 2], " or b) substitute it with 0?")
            edit.couple_YN <- readline(prompt = "(a/b)")

            j <- 0

            while (edit.couple_YN != "b" && edit.couple_YN != "B" &&
                   edit.couple_YN != "a" && edit.couple_YN != "A") {
              j <- j + 1
              if (j == 5)
                stop("The number of attempts has been exceeded")
              edit.couple_YN <- readline(prompt = "Please press a or b ")
            }

            if (edit.couple_YN == "b" || edit.couple_YN == "B") {
              couples.st2[i, k + 2] <- 0
            } else {
              couples.st2[i, k + 2] <- couples.st[i, k + 2]
            }
          }

          # end cicle on k
        }

        if (identical(matrix(0, nrow = 1, ncol = (ncol(couples.st) -
                                                  2)), as.integer(couples.st2[i, 3:ncol(couples.st)])) ==
            TRUE) {
          stop("There is at least one spatial couple with no specification of temporal lags.")

        }


      }
      #= end cicle on i

      #= Check on columns and rows with non-zero values

      ii <- -2
      jj <- 0
      iflag <- matrix(0, nrow = (nrow(sp.couples) / 3), ncol = ((n.temp / 3) / 2))
      for (i in 1:(nrow(sp.couples) / 3)) {
        kk <- -5
        ii <- ii + 3
        for (k in 1:((n.temp / 3) / 2)) {
          kk <- kk + 6
          col_tlag <- c((kk + 2), (kk + 4), (kk + 6))
          count_nozero_col <- colSums(couples.st2[ii:(ii + 2),
                                                  col_tlag] != 0)
          count_nozero_row <- rowSums(couples.st2[ii:(ii + 2),
                                                  col_tlag] != 0)


          if (match(3, count_nozero_col, nomatch = 0) != 0 && match(3,
                                                                    count_nozero_row, nomatch = 0) != 0) {
            iflag[i, k] <- 1
            jj <- jj + 1
          }

        }
      }
      if (jj == 0) {
        stop("No temporal lags have been specified.")
      }

      ii <- -2
      for (i in 1:(nrow(sp.couples) / 3)) {
        kk <- -5
        ii <- ii + 3
        for (k in 1:((n.temp / 3) / 2)) {
          kk <- kk + 6
          if (iflag[i, k] == 0) {
            couples.st2[ii:(ii + 2), (kk + 2):(kk + 7)] <- 0

          }
        }
      }
      # End check on columns and rows with non-zero values

    }
    couples.st <- couples.st2



  }
  #=================================================#
  #= End test on the type of model (typetest>=4)   =#
  #=================================================#


  message("*****************************************************************")
  message("*** One of the slots of couples class has been created.       ***")
  message("*** This is called @couples.st and contains the couples of    ***")
  message("*** the spatial points and the temporal lags to be analyzed.  ***")
  message("*** This output is an object for further functions.           ***")
  message("*****************************************************************")

  tl.couples <- unique(as.vector(couples.st[, -(1 : 2)]))
  tl.couples <- tl.couples[tl.couples != 0]
  sp.couples.nm <- data.frame(matrix(NA, nrow = nrow(sp.couples), ncol = 2))
  for(i in 1:(nrow(sp.couples))){
    sp.couples.nm[i,1] <- sel.staz[sp.couples[i,1]]
    sp.couples.nm[i,2] <- sel.staz[sp.couples[i,2]]
  }
  sp.couples.nm <- cbind(sp.couples,sp.couples.nm)
  colnames(sp.couples.nm) <- c("id.1", "id.2", "point.1", "point.2")

  new("couples", couples.st = couples.st,
      sel.staz = sel.staz,
      sp.couples = sp.couples.nm,
      tl.couples = tl.couples)
}
#' @include sepindex.R
NULL
