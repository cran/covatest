#' Class "couples"
#'
#' A class for spatial points and the corresponding temporal lags to be
#' analyzed in order to test some covariance properties and some well known
#' classes of space-time covariance functions models
#'
#' @slot couples.st matrix, in which the first two columns contain the
#' couples of spatial points (denoted with order numbers) to be analyzed
#' and the other columns the temporal lags associated with each spatial couples
#' @slot sel.staz numeric or character, contains the ID codes of the
#' selected spatial points
#' @slot sp.couples data.frame, contains the couples of order numbers
#' associated with the spatial points to be analyzed and the couples of the
#' ID codes
#' @slot tl.couples numeric, contains the temporal lags associated to the
#' couples of the selected spatial points
#' @slot typetest character; contains the code of the test to be performed
#'
#' @rdname couples-class
#' @exportClass couples
#' @export
setClass("couples", slots = c(couples.st = "matrix",
                                         sel.staz = "ANY",
                                         sp.couples = "data.frame",
                                         tl.couples = "numeric",
                                         typetest = "character"))

#' @param sel.staz vector, the sequence of ID codes which denote the spatial
#' points to be analyzed
#' @param sp.couples.in two-column matrix: rows corresponding to the couples
#' of different spatial points, chosen among the ones fixed in \code{sel.staz}
#' argument, to be compared
#' @param t.couples.in vector of only positive (negative) temporal lags to be
#' analyzed. The corresponding negative (positive) temporal lags are included
#' authomatically for \code{typetest = "sym", "sep", "tnSep"}. If some temporal
#' lags, corresponding to some couples of spatial points, are not required for
#' the specific test, they can be set equal to zero, through the specific
#' \code{\link{setzero}} method
#' @param typetest character, set \code{typetest = "sym"} for symmetry test
#' (default choice), \code{typetest = "sep"} for separability test,
#' \code{typetest = "tnSep"} for type of non separability test,
#' \code{typetest = "productSum"} for the test on the product-sum class of models,
#' \code{typetest = "intProduct"} for the test on the integrated product class
#' of models, \code{typetest = "gneiting"} for the test on the Gneiting class
#' of models
#' @param typecode type of object, i.e. numeric() or character(), specifies the
#' type of codification of the spatial points in the \code{data frame} or in the
#' STFDF/STSDF
#'
#' @details
#' It is important to point out that:
#' \itemize{
#' \item both positive and negative temporal lags are automatically considered in
#' the slot \code{@couples.st} and \code{@tl.couples} for symmetry test
#' (\code{typetest = "sym"}), separability test (\code{typetest = "sep"}) and
#' type of non separability tests
#' (\code{typetest = "tnSep"}). If the symmetry hyphotesis has not been rejected, only
#' positive temporal lags might be considered for the test on separability and type
#' of non separability (\code{typetest = "sep"} and \code{typetest = "tnSep"}), hence the
#' specific \code{\link{setzero}} method must be used to set the negative temporal
#' lags equal to zero
#'
#' \item for \code{typetest = "tnSep"} the temporal lags should be chosen according to
#' the results of the sample non separability ratios, plotted through a boxplot
#' classified for temporal lags (see \linkS4class{sepindex} for more details)
#'
#' \item for model tests (\code{typetest} equal to \code{"productSum"},
#' \code{"intProduct"} and \code{"gneiting"}), the number of analyzed spatial
#' points must be used to create at least 3 spatial couples or multiple of 3,
#' such that each triplet satisfies the condition
#' \deqn{||\mathbf{h}_{1}||^{2\gamma}- ||\mathbf{h}_{2}||^{2\gamma} = ||\mathbf{h}_{2}||^{2\gamma}-||\mathbf{h}_{3}||^{2\gamma}}{||h_1||^{2\gamma} - ||h_2||^{2\gamma} = ||h_2||^{2\gamma} - ||h_3||^{2\gamma}}
#' where \eqn{\gamma \in ]0,1]} only for \code{typetest = "intProduct"}
#' and \code{"gneiting"}.
#' The number of positive temporal lags must be at least 3, or multiple
#' of 3, too. The condition \deqn{u_{1}^{2\alpha}-u_{2}^{2\alpha}=u_{2}^{2\alpha}-u_{3}^{2\alpha}}{u_1^{2\alpha} - u_2^{2\alpha} = u_2^{2\alpha} - u_3^{2\alpha}}
#' where \eqn{\alpha \in ]0,1]} must be satisfied for each triplet
#' (only for \code{typetest = "intProduct"} and \code{"gneiting"}), as clarified
#' in Cappello et al., 2018. The values of \eqn{\gamma} and \eqn{\alpha} are usually fixed
#' equal to 0.5 or 1 according that the behavior near the origin of the spatial
#' and temporal marginal covariograms is linear or quadratic, respectively.
#' Note that for each spatial triplet and each temporal triplet, 6 contrasts can
#' be defined. However, for \code{typetest = "intProduct"} (test on the integrated
#' model) the user has to set arbitrarily one temporal lag equal to zero for each
#' spatial triplet in order to delete redundant contrasts, through the specific
#' \code{\link{setzero}} method
#' }
#'
#' @note{
#' Errors occur if
#' \itemize{
#' \item some spatial points, given in the sequence at the beginning of the
#' function, have not been used to generate the couples of spatial points
#'
#' \item there is at least one spatial couple with no specification of
#' temporal lags
#'
#' \item no temporal lags have been specified
#'
#' \item the number of spatial points fixed in \code{sel.staz} is less than 2
#'
#' \item the construction of the \code{sp.couples.in} is not consistent with the
#' test to be performed
#' }
#' }
#'
#' @seealso \code{\link{setzero}}
#'
#' @references
#' Cappello, C., De Iaco, S., Posa, D., 2020, {covatest}: An {R} Package for
#' Selecting a Class of Space-Time Covariance Functions.
#' Journal of Statistical Software, \bold{94(1)} 1--42.
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
#' t.couples.in = t.couples.in.sym, typetest = "sym", typecode = character())
#'
#' ### methods for couples
#' #1. show
#' couples.sym
#'
#' #2. [ extract
#' couples.sym[3, by.row = FALSE]
#' couples.sym[3, by.row = TRUE]
#'
#' #3. summary
#' summary(couples.sym)
#'
#' @rdname couples-class
#' @export
couples <- function(sel.staz, sp.couples.in, t.couples.in, typetest = "sym", typecode = numeric()) {


  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) abs(x -
                                                                round(x)) < tol

  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}

  ### SOME CHECKS ON THE ARGUMENTS ###

  if (is.character(typetest) == FALSE) {
    message("Start error message. The argument for typetest is not admissible.")
    stop("End error message. Stop running.")
  }

  if (typetest != "sym" && typetest != "sep"  && typetest != "tnSep" && typetest != "productSum"
      && typetest != "intProduct" && typetest != "gneiting") {
    message("Start error message. The argument for typetest is not admissible.")
    stop("End error message. Stop running.")
  }

  if (typeof(typecode) != "double" && typeof(typecode) != "character") {
    message("Start error message. The argument for typecode is not admissible.")
    stop("End error message. Stop running.")
  }

  if (is.vector(sel.staz) == FALSE || length(sel.staz) <2) {
    message("Start error message. The argument sel.staz must be a vector of at least 2 elements.")
    stop("End error message. Stop running.")
  }

  if (typeof(typecode) !=  typeof(sel.staz)) {
    message("Start error message. The type of data in sel.staz are not consistent with the declared typecode.")
    stop("End error message. Stop running.")
  }

  if (is.matrix(sp.couples.in) == FALSE || ncol(sp.couples.in) != 2) {
    message("Start error message. The argument sp.couples.in must be a matrix of 2 column. Please revise appropriately the argument sp.couples.in and run the function again.")
    stop("End error message. Stop running.")
  }

  if (is.vector(t.couples.in) == FALSE || is.numeric(t.couples.in) == FALSE
      || match(0, t.couples.in, nomatch = 0) != 0) {
    message("Start error message. The argument t.couples.in must be a numeric vector, with no zeros. Please revise appropriately the argument t.couples.in and run the function again.")
    stop("End error message. Stop running.")
  }

  t.couples.in <- abs(t.couples.in)
  if (match("TRUE", duplicated(t.couples.in, fromLast = TRUE), nomatch = 0) != 0) {
    message("Start error message. The argument t.couples.in must contains different temporal lags in absolut value. Please revise appropriately the argument t.couples.in and run the function again.")
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


  t.couples.in <- matrix(data = t.couples.in, ncol=1)

  if (type.test <=3){
    t.couples.in <- cbind(t.couples.in, -t.couples.in)
    colnames(t.couples.in) <- NULL
  }

  if (type.test >=4){
    t.couples.in <- cbind(t.couples.in, 0)
    colnames(t.couples.in) <- NULL
  }


  ### SPATIAL COUPLES ###

  ns <- length(sel.staz)

  if (ns < 2) {
    message("Start error message. The number of spatial points to be analyzed must be at least equal to two.")
    stop("End error message. Stop running.")
  }



  if (type.test >= 4) {
    message("***********************************************************************************************************")
    message("For the test on the type of class of models, note that:")
    message("- the spatial points must be used to create at least 3 spatial couples or multiple of 3 spatial couples ")
    message("- each triplet of different spatial points will be used for the spatial comparison. ")
    message("- the number of positive temporal lags must be at least 3 or multiple of 3 ")
    message("- each triplet of different temporal lags will be used for the temporal comparison.")
    message("A condition on the considered spatial and temporal lags must be satisfied for each triplet. Please see the ")
    message("  manual for more details and check this condition. ")
    message("*************************************************************************************************************")

  }

 # message("Enter a sequence of ID codes which denote the spatial points to be analyzed: ")
 if(mode(sel.staz) != mode(typecode)) {
   message("Start error message. The ID codes are not consistent with the typecode.")
   stop("End error message. Stop running.")}

  j <- 0

  while (anyDuplicated(sel.staz) >= 1) {
    message("Start error message. There are duplicates in the sequence of ID codes. Please revise appropriately the argument sel.staz and run the function again.")
    stop("End error message. Stop running.")
  }


  #message("********************************************************************************")
  #message("*** One of the slots of couples class has been created.                      ***")
  #message("*** This is called @sel.staz and contains the spatial points to be analyzed. ***")
  #message("*** This first output is an argument for further functions.                  ***")
  #message("********************************************************************************")

  nc <- nrow(sp.couples.in)
  if(mode(sp.couples.in) != mode(typecode)) {
    message("Start error message. The ID codes in sp.couples.in are not consistent with the typecode. Please revise appropriately the argument sp.couples.in and run the function again.")
    stop("End error message. Stop running.")
    }

  for (i in 1:nc) {
    if (sp.couples.in[i,1] == sp.couples.in[i,2]){
      message("Start error message. Some couples of points are not admissible: the points must be different and must be chosen among the selected spatial points. Please revise appropriately the argument sp.couples.in and run the function again.")
      stop("End error message. Stop running.")
    }
  }
  # check on spatial couples for symmetry, non-separability and type of nonseparability
  if (type.test < 4) {
    if (match(TRUE,duplicated(sp.couples.in),nomatch = 0) != 0){
    message("Start error message. Some couples of points are duplicated. Please revise appropriately the argument sp.couples.in and run the function again.")
      stop("End error message. Stop running.")
    }
    for (i in 1:(nc-1)) {
      for (j in (i+1):nc) {
          sp.couples.perm <- c(sp.couples.in[j,2],sp.couples.in[j,1])
          sp.couples.perm <- rbind(sp.couples.perm,sp.couples.in[i,])
          if (match(TRUE,duplicated(sp.couples.perm),nomatch = 0) != 0){
            message("Start error message. Some couples are only permuted and then are considered duplicated. Please revise appropriately the argument sp.couples.in and run the function again.")
            stop("End error message. Stop running.")
        }
      }
    }
  }

  # check on spatial couples for the tests on models
  if (type.test >= 4) {
    ns_multiple <- as.integer(nc / 3)
    if ((ns_multiple * 3) != nc) {
      message("Start error message. The total number of spatial couples is not consistent with the number of couples required by the test. Please revise appropriately the arguments and run the function again.")
      stop("End error message. Stop running.")
    }
  }

  # check WITHIN the triplets where three couples are replicated.
  if (type.test >= 4) {
    k <- 1
    for (j in 1:ns_multiple) {
      sp.couples.transf1 <- sp.couples.in[k,]
      sp.couples.transf2 <- sp.couples.in[k+1,]
      sp.couples.transf3 <- sp.couples.in[k+2,]

      if (length(setdiff(sp.couples.transf1, sp.couples.transf2)) == 0 &&
          length(setdiff(sp.couples.transf1, sp.couples.transf3)) == 0 &&
          length(setdiff(sp.couples.transf2, sp.couples.transf1)) == 0 &&
          length(setdiff(sp.couples.transf3, sp.couples.transf1)) == 0 ){
        message("Start error message. Some triplets of couples are not valid. There are triplets where the three couples are replicated.")
        message("Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }
      k <- k + 3
    }
  }

  # check WITHIN the triplets where two couples are replicated.
  if (type.test >= 4) {
    k <- 1
    for (j in 1:ns_multiple) {
      sp.couples.transf1 <- sp.couples.in[k:(k+1),]
      sp.couples.transf2 <- sp.couples.in[(k+1):(k+2),]
      sp.couples.transf3 <- rbind(sp.couples.in[k,],sp.couples.in[k+2,])

      sp.couples.in.triplet1 <-matrix(sapply(1:nrow(sp.couples.transf1), function(i)(sp.couples.transf1[i])), ncol = 4, byrow = T)
      sp.couples.in.triplet2 <-matrix(sapply(1:nrow(sp.couples.transf2), function(i)(sp.couples.transf2[i])), ncol = 4, byrow = T)
      sp.couples.in.triplet3 <-matrix(sapply(1:nrow(sp.couples.transf3), function(i)(sp.couples.transf3[i])), ncol = 4, byrow = T)

      if ((length(setdiff(sp.couples.in.triplet1, sp.couples.in.triplet2)) == 0) && (length(setdiff(sp.couples.in.triplet2, sp.couples.in.triplet1)) == 0)){
        message("Start error message. Some triplets of couples are not valid. There are triplets where two couples are replicated.")
        message("Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }
      if ((length(setdiff(sp.couples.in.triplet1, sp.couples.in.triplet3)) == 0) && (length(setdiff(sp.couples.in.triplet3, sp.couples.in.triplet1)) == 0)){
        message("Start error message. Some triplets of couples are not valid. There are triplets where two couples are replicated.")
        message("Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }
      if ((length(setdiff(sp.couples.in.triplet2, sp.couples.in.triplet3)) == 0) && (length(setdiff(sp.couples.in.triplet3, sp.couples.in.triplet2)) == 0)){
        message("Start error message. Some triplets of couples are not valid. There are triplets where two couples are replicated.")
        message("Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }
      k <- k + 3
    }
  }

  # check AMONG triplets with the same couples in the same order.
  if (type.test >= 4) {
    if(ns_multiple > 1){
    sp.couples.in.triplet <-matrix(sapply(1:nrow(sp.couples.in), function(i)(sp.couples.in[i,])), ncol = 6, byrow = T)
    if (match(TRUE,duplicated(sp.couples.in.triplet),nomatch = 0) != 0){
      message("Start error message. Some triplets of couples are duplicated. Please revise appropriately the argument sp.couples.in and run the function again.")
      stop("End error message. Stop running.")
    }
    }
    }

    # check AMONG triplets with 1) exchanges the first and third couple of the triplet, or 2) permutations in the triplet or 3) both
    # 1) exchanges the first and third couple of the triplet
     if (type.test >= 4) {
       if(ns_multiple > 1){
    k <- 1
    for (j in 2:ns_multiple) {
          sp.couples.transf <-sp.couples.in
          sp.couples.transf[k+2,] <- sp.couples.in[k,]
          sp.couples.transf[k,] <- sp.couples.in[k+2,]
          sp.couples.in.triplet <- matrix(sapply(1:nrow(sp.couples.transf), function(i)(sp.couples.transf[i,])), ncol = 6, byrow = T)

          if (match(TRUE,duplicated(sp.couples.in.triplet),nomatch = 0) != 0){
            message("Start error message. Some triplets of couples can be considered duplicated, since at least two triples differ just for the position of the first and the third couple. Please revise appropriately the argument sp.couples.in and run the function again.")
            stop("End error message. Stop running.")
          }
      k <- k + 3
    }
    # 2) permutation in the couples of the triplet
       k <- 1
    for (j in 2:ns_multiple) {
      sp.couples.transf <-sp.couples.in
      sp.couples.transf[k,] <- c(sp.couples.in[k,2],sp.couples.in[k,1])
      sp.couples.transf[k+1,] <- c(sp.couples.in[k+1,2],sp.couples.in[k+1,1])
      sp.couples.transf[k+2,] <- c(sp.couples.in[k+2,2],sp.couples.in[k+2,1])
      sp.couples.in.triplet <-matrix(sapply(1:nrow(sp.couples.transf), function(i)(sp.couples.transf[i,])), ncol = 6, byrow = T)

      if (match(TRUE,duplicated(sp.couples.in.triplet),nomatch = 0) != 0){
        message("Start error message. Some triplets of couples can be considered duplicated, since at least two triples are characterized by the same couples but taken in a different order. Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }
      k <- k + 3
    }
       # 3) both exchange the first and third couple of the triplet and permutation in the triplet
       k <- 1
       for (j in 2:ns_multiple) {
         sp.couples.transf <-sp.couples.in
         sp.couples.transf[k,] <- c(sp.couples.in[k+2,2],sp.couples.in[k+2,1])
         sp.couples.transf[k+1,] <- c(sp.couples.in[k+1,2],sp.couples.in[k+1,1])
         sp.couples.transf[k+2,] <- c(sp.couples.in[k,2],sp.couples.in[k,1])
         sp.couples.in.triplet <-matrix(sapply(1:nrow(sp.couples.transf), function(i)(sp.couples.transf[i,])), ncol = 6, byrow = T)

         if (match(TRUE,duplicated(sp.couples.in.triplet),nomatch = 0) != 0){
           message("Start error message. Some triplets of couples can be considered duplicated, since:")
           message(" at least two triples differ just for the position of the first and the third couple;")
           message(" or at least two triples are characterized by the same couples but taken in a different order.")
           message(" Please revise appropriately the argument sp.couples.in and run the function again.")
           stop("End error message. Stop running.")
         }
         k <- k + 3
       }
     }
  }

  # check AMONG triplets where three couples are replicated.
  # check on triplets that are exactly equal or are equal except for the order in which the spatial points are coupled.
  if (type.test >= 5) {
    if(ns_multiple > 1){
    sp.couples.in.triplet <-matrix(sapply(1:nrow(sp.couples.in), function(i)(sp.couples.in[i,])), ncol = 6, byrow = T)
    for (j in 1:(ns_multiple-1)) {
      jj <- (j-1)*3 + 1
      for (l in (j+1):ns_multiple) {
        ll <- (l-1)*3 + 1
        rbind(sp.couples.in[jj,],sp.couples.in[jj+1,],sp.couples.in[jj+2,],sp.couples.in[ll,])
      if ((length(setdiff(sp.couples.in.triplet[j,], sp.couples.in.triplet[l,])) == 0) &&
          (length(setdiff(sp.couples.in.triplet[l,], sp.couples.in.triplet[j,])) == 0) &&
          ((length(setdiff(sp.couples.in[jj,], sp.couples.in[ll,])) == 0) ||
           (length(setdiff(sp.couples.in[jj,], sp.couples.in[ll+1,])) == 0) ||
           (length(setdiff(sp.couples.in[jj,], sp.couples.in[ll+2,])) == 0)) &&
          ((length(setdiff(sp.couples.in[jj+1,], sp.couples.in[ll,])) == 0) ||
           (length(setdiff(sp.couples.in[jj+1,], sp.couples.in[ll+1,])) == 0) ||
           (length(setdiff(sp.couples.in[jj+1,], sp.couples.in[ll+2,])) == 0)) &&
          ((length(setdiff(sp.couples.in[jj+2,], sp.couples.in[ll,])) == 0) ||
           (length(setdiff(sp.couples.in[jj+2,], sp.couples.in[ll+1,])) == 0) ||
           (length(setdiff(sp.couples.in[jj+2,], sp.couples.in[ll+2,])) == 0))
          ){
        message("Start error message. Some triplets of couples are not valid.")
        message("There might be 1) at least two triplets with the same spatial couples (in different order) or ")
        message("2) the same spatial points are coupled in such a way that the condition on the lags is not satisfied.")
        message("Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")

      }
    }
    }
    }
    }
    # if (type.test == 4) {
    #   if(ns_multiple > 1){
    # sp.couples.in.triplet <-matrix(sapply(1:nrow(sp.couples.in), function(i)(sp.couples.in[i,])), ncol = 6, byrow = T)
    #  for (j in 1:(ns_multiple-1)) {
    #   for (l in (j+1):ns_multiple) {
    #     if ((length(setdiff(sp.couples.in.triplet[j,], sp.couples.in.triplet[l,])) == 0) && (length(setdiff(sp.couples.in.triplet[l,], sp.couples.in.triplet[j,])) == 0)){
    #       message("Start error message. Some triplets of couples are not valid.")
    #       message("There might be 1) at least two triplets with the same spatial couples (in different order) or ")
    #       message("2) at least two triplets that are equal except for the order in which the spatial points are coupled.")
    #       message("Please revise appropriately the argument sp.couples.in and run the function again.")
    #       stop("End error message. Stop running.")
    #     }
    #    }
    #   }
    #  }
    # }

  sp.couples <- sp.couples.in
  if (length(setdiff(sel.staz, sp.couples)) >= 1) {
    message("Start error message. The following spatial points have not been used to generate the couples of spatial points.")
    print(setdiff(sel.staz, sp.couples))
    message("Please revise appropriately the argument sp.couples.in and run the function again.")
    stop("End error message. Stop running.")
  }
  sp.couples2 <- sp.couples
  sp.couples <- matrix(match(sp.couples, sel.staz), nrow = nc, ncol = 2)


  ### TEMPORAL COUPLES ###

  #===================================================================#
  #= Start defining temporal couples for                             =#
  #= test of symmetry (type.test=0), separability (type.test=1) and  =#
  #= type of non separability (type.test=2)                          =#
  #===================================================================#

  if (type.test == 1 || type.test == 0 || type.test == 2) {
    n.temp <- 2 * nrow(t.couples.in)

    for (i in 1:(n.temp/2)) {
      t.couples <- t.couples.in[i,]
      if((t.couples[2] + t.couples[1] != 0) || (t.couples[1] == 0)) {
        message("Start error message. This couple of temporal lags is not admissible: the temporal lags must be equal in absolut value and their absolut values must be greater than zero. Please revise appropriately the argument t.couples.in and run the function again.")
        stop("End error message. Stop running.")
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

        if(dupli.couple[1] == TRUE || dupli.couple.perm[1] == TRUE) {
          message("Start error message. The couple (# ", i, ") of temporal lags already exists or have been included in a different order. Please revise appropriately the argument t.couples.in and run the function again.")
          stop("End error message. Stop running.")
        }
      }
    }

    for (i in 1:(n.temp/2)) {
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

  }
  #==============================================================#
  #= end test on type=0 (symmetry), type=1 (separability) and   =#
  #= type=2(type of non separability)                           =#
  #==============================================================#


  #=================================================#
  #== Start defining temporal couples for          =#
  #== test on the type of model (type.test>=4)     =#
  #=================================================#

  if (type.test >= 4) {
    n.temp <- 2 * nrow(t.couples.in)
    t.couples <- matrix(0, nrow = 2, ncol = 1)

    nt_multiple <- as.integer(n.temp/3)
    if ((nt_multiple * 3) != n.temp) {
     message("Start error message. The total number of temporal lags is not consistent with the number of lags required by the test. Please revise appropriately the argument t.couples.in and run the function again.")
     stop("End error message. Stop running.")
      }

    for (i in 1:(n.temp/2)) {

      t.couples[1] <- t.couples.in[i,1]
      t.couples[2] <- 0

      while (t.couples[1] <= 0) {
        message("Start error message. This temporal lag is not admissible: the temporal lags must be greater than zero. Please revise appropriately the argument t.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }

      if (i == 1) {
        tl.couples <- matrix(t.couples, nrow = 1, ncol = 2, byrow = TRUE)
      }
      if (i > 1) {

        tl.couples2 <- matrix(t.couples, nrow = 1, ncol = 2, byrow = TRUE)
        tl.couples <- rbind(tl.couples, tl.couples2)
        dupli.couple <- duplicated(tl.couples, fromLast = T)


        while (dupli.couple[1] == TRUE) {
          message("Start error message. Some temporal lags are duplicated. Please revise appropriately the argument t.couples.in and run the function again.")
          stop("End error message. Stop running.")
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

    message("********************************************************************")
    message("*** This is the slot called @couples.st and contains the couples ***")
    message("*** of the spatial points and the temporal lags to be analyzed.  ***")
    message("********************************************************************")
    message("*** Some temporal lags not used can be set equal to 0 through    ***")
    message("*** the specific setzero method                                  ***")
    message("********************************************************************")

    # check if there are at least two spatial triplets, where two couples are replicated
    if(ns_multiple > 1){
      index.couples1 <- matrix(NA,ncol = 2)
      index.couples2 <- matrix(NA,ncol = 2)
      icheck <- 0
      kk <- 0
      for (j in 1:(ns_multiple-1)) {
      jj <- ((j-1)*3)+1
      for (j2 in 1:2) {
        for (j3 in (j2+1):3) {
          sp.couples.transf1.1 <- sp.couples.in[jj+j2-1,]
          sp.couples.transf1.2 <- sp.couples.in[jj+j3-1,]

          for (l in (j+1):ns_multiple) {
            ll <- ((l-1)*3)+1
            for (l2 in 1:2) {
              for (l3 in (j2+1):3) {
                sp.couples.transf2.1 <- sp.couples.in[ll+l2-1,]
                sp.couples.transf2.2 <- sp.couples.in[ll+l3-1,]
                if (length(setdiff(sp.couples.transf1.1, sp.couples.transf2.1)) == 0 &&
                    length(setdiff(sp.couples.transf1.2, sp.couples.transf2.2)) == 0){
                  icheck <- icheck +1
                  kk <- kk+1
                  if(kk == 1){
                    index.couples1[1,] <- c(j,l)
                  }else{
                    index.couples1 <- rbind(index.couples1,c(j,l))
                  }
                }
                if (length(setdiff(sp.couples.transf1.1, sp.couples.transf2.2)) == 0 &&
                    length(setdiff(sp.couples.transf1.2, sp.couples.transf2.1)) == 0){
                  icheck <- icheck +1
                  kk <- kk+1
                  if(kk == 1){
                    index.couples1[1,] <- c(j,l)
                  }else{
                    index.couples1 <- rbind(index.couples1,c(j,l))
                  }
                }
              } #end cicle on l3
            } #end cicle on l2
          } #end cicle on j3
        } #end cicle on j2
      } #end cicle on l
    } #end cicle on j


    # check if there are at least two spatial triplets, where one couple is replicated
      icheck2 <- 0
      icheck3 <- 0
      kk <- 0
      for (j in 1:(ns_multiple-1)) {
        for (l in (j+1):ns_multiple) {
          for(jj in (((j-1)*3)+1):(((j-1)*3)+3)){
            for(ll in (((l-1)*3)+1):(((l-1)*3)+3)){
              sp.couples.transf1 <- sp.couples.in[jj,]
              sp.couples.transf2 <- sp.couples.in[ll,]
              if (length(setdiff(sp.couples.transf1, sp.couples.transf2)) == 0 &&
                  length(setdiff(sp.couples.transf2, sp.couples.transf1)) == 0){
                  icheck3 <- icheck3 +1
                  if(length(setdiff(c(j,l),index.couples1)) != 0){
                   icheck2 <- icheck2 +1
                   kk <- kk+1
                   if(kk == 1){
                     index.couples2[1,] <- c(j,l)
                   }else{
                     index.couples2 <- rbind(index.couples2,c(j,l))
                   }
                   }
                }
            }
          }

        }
      }
      # index.couples <- rbind(index.couples1,index.couples2)
      # if (length(setdiff(c(1:ns_multiple), as.vector(index.couples))) == 0 ){
      #    message("Start error message. In all the triplets, there are at least two spatial triplets, where two couples are replicated. Please revise appropriately the argument sp.couples.in and run the function again.")
      #    stop("End error message. Stop running.")
      # }
      if(icheck == choose(ns_multiple,2)){
        message("Start error message. There are two couples replicated in all the triplets. Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }
      if(icheck3 == choose(ns_multiple,2)){
        message("Start error message. There is  one couple replicated in all the triplets. Please revise appropriately the argument sp.couples.in and run the function again.")
        stop("End error message. Stop running.")
      }
      if(icheck >= 1){
        message("Note that there are at least two spatial triplets, where TWO couples are replicated.")
        message("Replications have been found in the following triplets: ")
        for(j in 1: nrow(index.couples1)){
        message(index.couples1[j,1], " ", index.couples1[j,2])
          }
        message("Remember to use the method setzero in order to avoid singularity problems in performing the tests.")
        message("See the manual for further details.")
        message("  ")
      }
      if(icheck2 >= 1){
        message("Note that there are at least two spatial triplets, where ONE couple is replicated.")
        message("Replications have been found in the following triplets: ")
        for(j in 1: nrow(index.couples2)){
          message(index.couples2[j,1], " ", index.couples2[j,2])
        }
        message("Remember to use the method setzero in order to avoid singularity problems in performing the tests.")
        message("See the manual for further details.")
      }
    }
    }
  #=================================================#
  #= End test on the type of model (type.test>=4)   =#
  #=================================================#


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
      tl.couples = tl.couples,
      typetest = typetest)
}
#' @include sepindex.R
NULL
#' @param object object of class \code{couples} for methods \code{show} and \code{summary}
#' @rdname couples-class
#' @aliases couples-method
#' @aliases show
#' @export
setMethod(f="show", signature="couples", definition=function(object) {
cat("Spatio-temporal lags defined throught the objects of the class 'couples'", "\n")
  couples.st2 <- matrix(object@couples.st, nrow = nrow(object@couples.st), ncol = ncol(object@couples.st))
  couples.st2[, 1:2] <- as.matrix(object@sp.couples[, 3:4], nrow = nrow(object@sp.couples))
  rownames(couples.st2) <- NULL
  colnames(couples.st2) <- NULL
  print(couples.st2)
  cat("\n")
  cat("Slot 'sel.staz':")
  cat("\n")
  print(object@sel.staz)
  cat("\n")
  cat("Slot 'sp.couples':")
  cat("\n")
  print(object@sp.couples)
  cat("\n")
  cat("Slot 'tl.couples':")
  cat("\n")
  print(object@tl.couples)
  cat("\n")
  cat("Slot 'typetest':")
  cat("\n")
  print(object@typetest)
}
)
#' @param x object of class \code{couples} for method \code{extract}
#' @param i index specifing rows or columns of the slot \code{@couples.st}.
#' Rows or columns depending on the logical parameter \code{by.row} to be set
#' @param by.row logical, if \code{TRUE} rows of the slot \code{@couples.st} are
#' selected (the temporal lags associated to the i-th spatial couple are given).
#' If \code{FALSE} (the default) columns of the slot \code{@couples.st} are
#' selected. In particular, the spatial couples associated to the i-th temporal
#' lag (i >= 3, temporal lags are stored from the third column) are given
#'
#' @rdname couples-class
#' @aliases couples-method
#' @aliases select
#' @export
setMethod(f="[", signature="couples", definition=function(x, i, by.row = FALSE) {
  if(by.row == FALSE){
    if(i <= 2){
      message("Start error message. The column selected does not contain the temporal lags. Please select a column greater than 2.")
      stop("End error message. Stop running.")
    }

    y <-as.data.frame(cbind(x@sp.couples[, 3:4], x@couples.st[, i]))
    names(y) <- NULL
    return(y)
  }
  if(by.row == TRUE){
    y <- as.data.frame(c(x@sp.couples[i, 3:4], x@couples.st[i, -(1:2)]))
    names(y) <- NULL
    return(y)
  }

}
)
#' @rdname couples-class
#' @aliases couples-method
#' @aliases summary
#' @export
#' @export
setMethod(f = "summary", signature = "couples",
          definition = function(object) {
            cat("Number of temporal lags = " , length(object@tl.couples) , "\n")
            cat("Number of spatial points involved = " , length(object@sel.staz) , "\n")
            cat("Number of spatial couples = " , nrow(object@couples.st), "\n")

          })
#' setzero
#'
#' Through the function {\link{couples}}, \code{m} spatial couples and \code{n}
#' temporal lags are provided, hence a set of \code{m x n} spatio-temporal lags
#' are defined. If some of these lags are not required for the specific test, they
#' can be set equal to zero by using the \code{setzero} method for object of class
#' \code{couples}
#'
#' @param x object of class \code{couples}
#' @param zero logical, if \code{TRUE} (the default) all negative temporal
#' lags are replaced with zero. If \code{x@typetest} is equal to
#' \code{"sym"} (symmetry test) the argument \code{setzero} is ignored because
#' both positive and negative temporal lags are required for computing the test
#' @param index two column matrix. Each row of the matrix \code{index} contains
#' the specific row and column, of the slot \code{@couples.st}, for which the
#' spatio-temporal covariance is not required
#' @param value numeric, the value to be replaced. Note that this method is reasonable
#' to be used only to replace a value equal to zero
#' @seealso \code{\link{couples}}
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
#' t.couples.in = t.couples.in.sym, typetest = "sym", typecode = character())
#'
#' zero.index <- matrix(data=c(1,3,1,4,2,5,2,6), ncol=2, byrow = TRUE)
#'
#' setzero(couples.sym, zero = FALSE, index = zero.index, value = 0)
#'
#' @export
#' @docType methods
#' @rdname setzero-methods
setGeneric("setzero", function(x, zero = TRUE, index = NULL, value) {
  standardGeneric("setzero")
})
#' @rdname setzero-methods
#' @aliases setzero,ANY,ANY-method
#' @export
setMethod(f="setzero", signature(x="couples"),
          definition=function(x, zero, index, value) {

            setzero <- zero

            if(setzero == TRUE && is.null(index) == FALSE){
              message("Start error message. Arguments not valid.")
              stop("End error message. Stop running.")
            }

            if(setzero == FALSE && is.null(index) == TRUE){
              message("Start error message. Arguments not valid.")
              stop("End error message. Stop running.")
            }

            if(setzero == TRUE && x@typetest == "sym"){
              message("Start error message. Argument setzero is ignored for symmetry test. See manual for details.")
              stop("End error message. Stop running.")
              setzero == FALSE
            }

            if(value != 0){
              message("Start error message. The value to be replaced has to be only zero.")
              stop("End error message. Stop running.")
            }

            couples.st2 <- x@couples.st
            t.lag <- (ncol(x@couples.st)-2)

            if(setzero == TRUE){
              for(i in 1:nrow(x@couples.st)){
                jj <- 1
                for(j in 1:(t.lag/2)){
                  jj <- jj + 2
                  couples.st2[i, jj+1] <- value
                }
              }
            }
            if(setzero == FALSE){
              for(i in 1:nrow(index)){
                if(index[i,1] <= 2 && index[i,2] <= 2){
                  message("Start error message. The editing of spatial points is not admissible.")
                  stop("End error message. Stop running.")
                }
                couples.st2[index[i,1],index[i,2]] <- value
              }





              ### Check on nonadmissible changes on the spatial points

              for (i in 1:nrow(x@sp.couples)) {

                #== Check on nonadmissible changes on the temporparal lags (typetest = 1
                #== and 2)
                if (x@typetest == "sep" || x@typetest == "tnSep") {
                  for (k in 1:t.lag) {

                    if ((x@couples.st[i, k + 2]) != (couples.st2[i, k + 2]) &&
                        (couples.st2[i, k + 2]) != 0) {
                      message("Start error message. Regarding the spatial couple (", x@sp.couples[i, 3], " ; ", x@sp.couples[i, 4], " the temporal lag ", couples.st2[i, jj + 2], " is not admissible.")
                      message("You can only substitute zeros to the existing lags just to esclude these temporal lags for this specific spatial couple. Please revise the arguments and run again.")
                      stop("End error message. Stop running.")

                    }
                  }
                }



                # #== Check on nonadmissible changes on the temporparal lags (typetest=0)
                #


                if (x@typetest == "sym") {
                  jj <- -1
                  for (j in 1:(t.lag / 2)) {
                    jj <- jj + 2
                    if (couples.st2[i, jj + 2] != abs(couples.st2[i, jj + 3])) {
                      message("Start error message. Regarding the spatial couple (", x@sp.couples[i, 3], " ; ", x@sp.couples[i, 4], ") the temporal lags (", couples.st2[i, jj + 2], " ; ", couples.st2[i, jj + 3], ") are not admissible")
                      message("For symmetry test, you need to substitute zero to both negative and positive temporal lags. Please revise the argument 'index' and run again.")
                      stop("End error message. Stop running.")
                    }
                  }
                }


                if (identical(as.integer(matrix(0, nrow = 1, ncol = (ncol(x@couples.st) -
                                                                     2))), as.integer(couples.st2[i, 3:ncol(x@couples.st)])) ==
                    TRUE) {
                  message("Start error message. There is at least one spatial couple with no specification of temporal lags.")
                  stop("End error message. Stop running.")
                }
              }


              n.temp <- (ncol(x@couples.st) - 2)
              if (x@typetest == "productSum" || x@typetest == "intProduct" ||x@typetest == "gneiting" ) {
                #= Check on nonadmissible changes on the spatial points

                for (i in 1:nrow(x@sp.couples)) {

                  #= Check on nonadmissible changes on the temporparal lags (typetest>=4)


                  for (k in 1:t.lag) {

                    if ((x@couples.st[i, k + 2]) != (couples.st2[i, k + 2]) &&
                        (couples.st2[i, k + 2]) != 0) {

                      message("Start error message. Regarding the spatial couple (", x@sel.staz[couples.st2[i,1]], "", x@sel.staz[couples.st2[i, 2]], " the temporal lag ", couples.st2[i, k + 2], " is not admissible.")
                      message("You can only substitute zeros to the existing lags just to esclude these temporal lags for this specific spatial couple. Please revise the arguments and run again.")
                      stop("End error message. Stop running.")

                    }

                    # end cicle on k
                  }

                  if (identical(matrix(0, nrow = 1, ncol = (ncol(x@couples.st) -
                                                            2)), as.integer(couples.st2[i, 3:ncol(x@couples.st)])) ==
                      TRUE) {
                    message("Start error message. There is at least one spatial couple with no specification of temporal lags.")
                    stop("End error message. Stop running.")

                  }


                }
                #= end cicle on i

                #= Check on columns and rows with non-zero values


                ii <- -2
                jj <- 0
                iflag <- matrix(0, nrow = (nrow(x@sp.couples) / 3), ncol = ((n.temp / 3) / 2))
                for (i in 1:(nrow(x@sp.couples) / 3)) {
                  kk <- -5
                  ii <- ii + 3
                  for (k in 1:((t.lag / 3) / 2)) {
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
                  message("Start error message. No temporal lags have been specified.")
                  stop("End error message. Stop running.")
                }

                ii <- -2
                for (i in 1:(nrow(x@sp.couples) / 3)) {
                  kk <- -5
                  ii <- ii + 3
                  for (k in 1:((t.lag / 3) / 2)) {
                    kk <- kk + 6
                    if (iflag[i, k] == 0) {
                      couples.st2[ii:(ii + 2), (kk + 2):(kk + 7)] <- 0

                    }
                  }
                }
                # End check on columns and rows with non-zero values
              }# End check on couples for models
            }

            x@couples.st <- couples.st2
            tl.couples <- unique(as.vector(x@couples.st[, -(1 : 2)]))
            tl.couples <- tl.couples[tl.couples != 0]
            x@tl.couples <- tl.couples
            #validObject(x)
            return(x)
          }
)
