#' Imports a GSLIB file in R
#'
#' Imports a \code{GSLIB} file and convert it into a data frame or in a \code{STFDF},
#' according to the standard of the \code{spacetime} package.
#'
#' @details The function requires the user to set some external arguments, that
#' is:
#' \enumerate{
#'    \item the GSLIB file name and its extension
#'    \item the number of variables in the file
#'    \item the code of missing values
#'    \item the column in which the coordinates of the spatial points are stored
#'    \item the number of different spatial points
#'    \item the flag for using/not using an existing identification ID for the
#'    spatial locations
#'    \item the column in which the time coordinates are stored
#'    \item the number of temporal observations for each spatial point
#'    \item the column in which the values of the variable are stored
#'    }
#' Moreover it is allowed the user to choose between two options for saving the
#' data: 1. save the data as an STFDF (\code{spacetime} class) 2. save the
#' data as a \code{data frame}. If the first option is selected (save the data
#' as a STFDF), it is also necessary to specify
#'   \itemize{
#'    \item the \code{Start Date} of the STFDF
#'    \item the interval of time between two temporal observations
#'    }
#'
#' The user must provide the aforementioned arguments concerning the
#' structure of the \code{GSLIB} file in order to convert it in a \code{data frame}
#' or in a STFDF
#'
#' @return object of the STFDF-class or \code{data frame}, which contains
#' coordinates of the spatial points, the identification code of the spatial
#' locations, the indentification code of the temporal observations and the observed
#' values of the variable
#'
#' @note As described by Remy et al. (2009), the \code{GSLIB} file is a \code{ASCII}
#' file organized by lines:
#' \itemize{
#' \item the first line gives a title
#' \item the second lines is a single number \code{n} indicating the number of
#' properties in the object (i.e. the number of columns of data)
#' \item the \code{n} following lines contain the name of each property
#' \item each remaining line contains the values of each properties (\code{n} values
#' per line) separated by spaces or tabulations, for each spatio-temporal points.
#' The order of the property values along each line is given by the order in
#' which the property names were entered
#' }
#'
#'
#' @references
#' Bivand, R. S., Pebesma, E., Gomez-Rubio, V., 2013,
#' \emph{Applied spatial data analysis with R}, Second edition. New York: Springer.
#' \url{http://www.asdar-book.org/}
#'
#' Pebesma, E.J., 2012, spacetime: Spatio-Temporal Data in R.
#' Journal of Statistical Software, \bold{51(7)} 1--30.
#' \url{http://www.jstatsoft.org/v51/i07/}
#'
#' Remy, N., Boucher, A., Wu, J, 2009, \emph{Applied Geostatistics with SGeMS:
#' A User's Guide}. Cambridge
#'
#' @seealso \code{\link[spacetime]{STFDF}}
#'
#' @importFrom methods as
#' @importFrom utils edit read.table
#' @importFrom stats cov pchisq pnorm var
#' @import spacetime
#'
#' @examples
#' ## The function requires to set some external arguments.
#' # In the GSLIB file, used as example, the measurements of PM10 in 13 rural
#' # background monitoring stations, in the period 2005-2006  (730 days), are
#' # stored. The information required to load the dataset concern:
#' # - the file name and its extension: PM10data.txt
#' # - the number of variables in the file: 6
#' # - the code of missing values: -999
#' # - the column in which the x coordinates are stored: 1
#' # - the column in which the y coordinates are stored: 2
#' # - the number of different spatial points: 13
#' # - the availability of an identification id for all the spatial points: y
#' # - the column in which the time coordinates are stored: 5
#' # - the number of temporal observations for each spatial point: 730
#' # - the column in which the values of the variable are stored: 6
#' # The user could choose to save the data as STFDF (option 1) or dataframe
#' # (option 2). If the option 1 will be choosen, it is also necessary to specify:
#' # - the Start Date: 2005-01-01
#' # - the interval of time between two temporal observations: day
#'
#' # To run the example, paste and copy the following lines
#' # (without the symbol '#') in the console
#' #
#' #datafile <- dataprep()
#' #PM10data.txt
#' # 5
#' #-999
#' #1
#' #2
#' #13
#' #y
#' #4
#' #730
#' #5
#' #1
#' #2005-01-01
#' #day
#' @export
dataprep <- function() # file, nvar, iclx, icly, iclt, iclvr, nt, ns, missingvalue file= name of
  # file to be imported ('nov2006.dat') nvar= number of variables (nvar+2=#
  # rows of strings to skip in reading the GsLib file) nt=number of temporal
  # points ns=number of spatial points missing=code for missing values
{
  message("*****************************************************************************")
  message("* Some information about the GSLIB file to be prepared are required.        *")
  message("* Please answer to the following questions.                                 *")
  message("*****************************************************************************")

  file <- readline(prompt = "Enter the GSLIB file name and its extension (i.e. demo.txt): ")
  nvar <- as.integer(readline(prompt = "Enter the the number of variables in the file: "))
  missing.v <- as.integer(readline(prompt = "Enter the code of missing values: "))
  iclx <- as.integer(readline(prompt = "Enter the column in which the x coordinates are stored: "))
  icly <- as.integer(readline(prompt = "Enter the column in which the y coordinates are stored: "))
  n.stat <- as.integer(readline(prompt = "Enter the number of different spatial points: "))
  code_YN <- readline(prompt = "Are the spatial points coded with an identification id? (Y/N)")
  iclt <- as.integer(readline(prompt = "Enter the column in which the time coordinates are stored: "))
  n.time <- as.integer(readline(prompt = "Enter the number of temporal observations for each spatial point: "))
  iclvr <- as.integer(readline(prompt = "Enter the column in which the values of the variable are stored: "))
  n.stpoint <- n.time * n.stat

  importFl <- read.table(file, skip = (nvar + 2), na.strings = missing.v)

  importFl[importFl == missing.v] <- NA
  message("There are two options for saving the data: 1. save the data as an STFDF (spacetime class) 2. save the data as a data frame. Please choose 1 or 2")
  code_12 <- as.integer(readline(prompt = "Please choose 1 or 2 "))
  if (code_12 != 1 && code_12 != 2)
    stop("The digit entered is not allowed.")
  if (code_12 == 2) {

    ### IMPORT GSLIB FILE ###

    importFl1 <- importFl[order(importFl[, iclx], importFl[, icly], importFl[,
                                                                             iclt]), ]


    if (n.stpoint != nrow(importFl1)) {
      stop("Stop running: error in data file. Maybe the number of temporal observations is not the same for each spatial point.")
    }

    if (code_YN == "N" || code_YN == "n") {
      ID_points <- c(rep(1:n.stat, each = n.time))
      importFl1 <- cbind(ID_points, importFl1)
      message("The spatial points have been coded by using consecutive numbers starting from 1.")
      message("The identification IDs have been written in the first column. ")
    }


    return(importFl1)

  }


  #=========================================#
  #=       CREATE STFDF (gstat package)    =#
  #=========================================#
  if (code_12 == 1) {


    #== SPATIAL DB ==#

    importFl <- importFl[order(importFl[, iclx], importFl[, icly], importFl[,
                                                                            iclt]), ]
    sp <- cbind(importFl[, iclx], importFl[, icly])
    sp <- unique(sp)
    colnames(sp) <- c("x", "y")
    sp2 <- sp::SpatialPoints(sp)

    #== TEMPORAL DB ==#

    s.time <- readline(prompt = "Enter the Start Date in the format YYYY-MM-DD %H:%M:%S (i.e. ): ")
    by.time <- readline(prompt = "Enter the interval of time between two temporal observations (i.e. day/min/...): ")
    data.time <- seq(from = as.POSIXct(s.time), by = by.time, length = n.time)



    #== DATI ST ==#

    importFl <- importFl[order(importFl[, iclx], importFl[, icly], importFl[,
                                                                            iclt]), ]
    mydata <- importFl[, iclvr]

    mydata <- matrix(mydata, ncol = n.time, byrow = TRUE)



    #== STFDF ==#
    stfdf <- STFDF(sp2, data.time, data.frame(variable = as.vector(as.matrix(mydata))))

    return(stfdf)

  }

}
