#' Imports a text file in R
#'
#' A function for importing a text file containing spatio-temporal data.
#' In particular, it (a) generates the spatial and temporal IDs,
#' (b) converts the time series of each spatial point (with non-existing values
#' for some dates) into a regularly spaced object
#' within the observed time period, by filling the missing dates with ‘NA’
#' (c) converts the data into a \code{STFDF}, according to the standard of the
#' \code{spacetime} package, or into a data frame
#'
#' @param file the name of the data file and its extension. The file is searched
#' in the current working directory, otherwise the absolute path has to be included
#' in the file name. Note that data for each spatial point and
#' each temporal point are given by row; each row of the file contains at least
#' the x and y coordinates of a spatial point, the temporal code (or date) and
#' the measurement of the variable of interest.
#' @param header logical value indicating whether the file contains the names
#' of the variables in the first line. If this argument is missing, header is
#' set to \code{FALSE} (default choice)
#' @param dec character used to indicate decimal points
#' @param sep field separator character. If sep = "" (default choice)
#' columns of the file are separated by white space or tabs (see \code{\link[utils]{read.table}}
#' for more details)
#' @param iclx numeric; the column in which the x-coordinate of the spatial points are
#' stored
#' @param icly numeric; the column in which the y-coordinate of the spatial points are
#' stored
#' @param iclt numeric; the column in which numeric temporal codes are stored.
#' This argument is provided  only if the \code{icldate} argument is not available:
#' \code{iclt} and \code{icldate}  are mutually exclusive. This argument is set
#' equal to \code{0}, if not available
#' @param icldate numeric vector to set the columns in which the dates are stored.
#' The user has to set \code{icl.date} if the date is stored in a single column, otherwise
#' the user has to specify the colunn in which the years (\code{iclty}), the months
#' (\code{icltm}) or the days (\code{icltd}) are stored separately.
#' This argument is set equal to \code{0} (default choice) if not available
#' @param icltime numeric vector to set the columns in which the time component
#' (hour, minute, second) of a date (if available) is stored. The user has to
#' set \code{icl.time} if the time is stored in a single column, otherwise the
#' user has to specify the colunn in which the hours (\code{icltH}), the minutes
#' (\code{icltM}) or the seconds (\code{icltS}) are stored separately.
#' This argument is set equal to \code{0} (default choice) if not available
#' @param iclvr numeric; the column in which the values of the variable are stored
#' @param iclsp numeric; the column in which the identification codes (IDs) for the
#' spatial locations are stored. This argument is set equal to \code{0} (default
#' choice) if IDs for the spatial locations are not available
#' @param missing.v the code used to indicate the presence of missing values in
#' the imported data. By default this argument is set equal to \code{NA}
#' @param save.as character, indicating the class of the data to be returned.
#' It is allowed to choose between two options for saving the file (\code{"STFDF"}
#' or \code{"data.frame"})
#' @param date.format vector, whose first element \code{date.format[1]} denotes the class of
#' the temporal component to be imported and the second one \code{date.format[2]} represents
#' the corresponding format. Note that the supported class of dates are
#' \code{"yearmon"}, \code{"yearqtr"},
#' \code{"Date"}, \code{"POSIX"} (see \pkg{Base}, \pkg{lubridate}, \pkg{zoo});
#' moreover the personalized options \code{"year"} and \code{"code"} are also admissible and
#' are used if the temporal coordinate is given by year or as a numerical code,
#' respectively. By default, the argument \code{date.format} is set equal to \code{("code", format = NA)}.
#' If the temporal component, provided for example in year and month, is given in
#' separeted columns in the text file,
#' the required format in \code{date.format[2]} is of the type \code{"\%Y \%m"};
#' in general the format requires the use of white spaces between two consecutive time units
#' @param bytime character, which denotes the time disaggregation of interest,
#' set \code{NA} (default choice) for numeric temporal code, otherwise \code{"\%Y"}
#' or \code{"\%y"} if values are taken by year, \code{"\%m"} if values are taken by
#' month, \code{"\%d"} if values are taken by day, \code{"\%q"} if values
#' are taken by quarter, \code{"\%H"} if values are taken by hour, \code{"\%M"}
#' if values are taken by minute and \code{"\%S"} if values are taken by seconds
#' @param tlag numeric; time increment\/lag between two temporal observations
#' @param time.zone character; time zone for dates with time component
#'
#' @return object of the \code{STFDF}-class or \code{data.frame}, which contains
#' coordinates of the spatial points, the spatial IDs, the temporal IDs, the dates
#' (if available in the input file) and the observed values of the variable of interest
#'
#'
#' @note
#' \itemize{
#' \item Uncomplete time series, for each spatial point, are filled with NA.
#' \item Some checks on the admissibility of the supported classes of dates are
#' implemented.
#' \item Time indexes for temporal points are coded for data.frame output by using
#' consecutive numbers starting from 1 (column 'timeIndex').
#' \item The spatial points are coded by using the string 'id' and the consecutive numbers
#' starting from 1 (column 'spatialIndex').
#' }
#'
#'
#' @references
#' Bivand, R. S., Pebesma, E., Gomez-Rubio, V., 2013,
#' \emph{Applied spatial data analysis with R}, Second edition. New York: Springer.
#' \url{http://www.asdar-book.org/}
#'
#' Grolemund, G, Wickham, H., 2011, Dates and Times Made Easy with lubridate.
#' Journal of Statistical Software, \bold{40(3)} 1--25.
#' \url{http://www.jstatsoft.org/v40/i03/}
#'
#' Pebesma, E.J., 2012, spacetime: Spatio-Temporal Data in R.
#' Journal of Statistical Software, \bold{51(7)} 1--30.
#' \url{http://www.jstatsoft.org/v51/i07/}
#'
#' Zeileis, A., Grothendieck. G., 2005, zoo: S3 Infrastructure for Regular and
#' Irregular Time Series.
#' Journal of Statistical Software, \bold{14(6)} 1--27.
#'
#' @seealso \code{\link[spacetime]{STFDF}}
#' @seealso \code{\link[utils]{read.table}}
#' @seealso \code{\link[zoo]{yearmon}}
#' @seealso \code{\link[zoo]{yearqtr}}
#' @seealso \code{\link[base]{Dates}} for dates without times
#' @seealso \code{\link[base]{DateTimeClasses}}
#'
#'
#' @importFrom methods as
#' @importFrom utils edit
#' @importFrom utils read.table
#' @importFrom stats cov
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats var
#' @import lubridate
#' @import zoo
#' @import spacetime
#'
#'
#'
#' @examples
#' #example 1: import a text file, with dates stored in a single column (the 4th)
#' # and fill missing time points in monthly time series, with time lag equal to one
#'
#' file_date <- system.file("extdata", "file_date.txt", package = "covatest")
#' db.date <- read.STdata(file = file_date, header = TRUE, iclx = 2, icly = 3, iclt = 0,
#' icldate = c(icl.date = 4, iclty = 0, icltm = 0, icltd = 0),
#' icltime = c(icl.time = 0, icltH =0, icltM = 0, icltS = 0),
#' iclvr = 5, iclsp = 1, missing.v = -99999, save.as = "data.frame",
#' date.format = c("Date", "%d-%m-%Y"), bytime = "%m", tlag = 1)
#'
#'
#' #example 2: import a text file, with dates and times stored in different columns
#' # (from the 4th to the 9th) and fill missing time points in hourly time series,
#' # with time lag equal to three
#'
#' file_datetime <- system.file("extdata", "file_datetime.txt", package = "covatest")
#' db.datetime <- read.STdata(file = file_datetime, header = TRUE, iclx = 2, icly = 3, iclt = 0,
#' icldate = c(icl.date = 0, iclty = 6, icltm = 5, icltd = 4),
#' icltime = c(icl.time = 0, icltH = 7, icltM = 8, icltS = 9),
#' iclvr = 10, iclsp = 1, missing.v = -99999, save.as = "data.frame",
#' date.format = c("POSIX", "%Y %m %d %H %M %S"), bytime = "%H", tlag = 3)
#'
#'
#' #example 3: import a text file, with dates and times stored in different columns
#' # (from the 4th to the 9th) and fill missing time points in quarterly time series,
#' # with time lag equal to one
#'
#' file_yq <- system.file("extdata", "file_yq.txt", package = "covatest")
#' db.yq <- read.STdata(file = file_yq, header = TRUE, iclx = 2, icly = 3, iclt = 0,
#' icldate = c(icl.date = 4, iclty = 0, icltm = 0, icltd = 0),
#' icltime = c(icl.time = 0, icltH =0, icltM = 0, icltS = 0),
#' iclvr = 5, iclsp = 1, missing.v = -99999, save.as = "data.frame",
#' date.format = c("yearqtr", "%Y-Q%q"), bytime = "%q", tlag = 1)
#'
#'
#' @rdname read.STdata
#' @export
read.STdata <- function(file, header = FALSE, dec = ".", sep = "", iclx, icly, iclt,
                     icldate = c(icl.date = 0, iclty = 0, icltm = 0, icltd = 0),
                     icltime = c(icl.time = 0, icltH = 0, icltM = 0, icltS = 0),
                     iclvr, iclsp = 0, missing.v = NA,
                     save.as = "data.frame", date.format = c("code",format = NA),
                     bytime = NA, tlag, time.zone = ""){

    ### SOME CHECKS ON THE ARGUMENTS ###
  is.scalar <- function (x){length(x) == 1L && is.vector(x, mode = "numeric")}
  if(time.zone == ""){
    if(is.na(Sys.timezone())){
      message("Start error message. The current time zone is not set in your system.")
      message("Please set the specific time zone by using the parameter time.zone.")
      stop("End error message. Stop running.")
    }else{
      time.zone <- Sys.timezone()
    }
  }
  if (is.scalar(iclx) == FALSE || is.scalar(icly) == FALSE  ||
      is.scalar(iclvr) == FALSE || is.scalar(iclsp) == FALSE) {
    message("Start error message. Some of the arguments iclsp, iclx, icly, iclt, iclvr are not numeric.")
    stop("End error message. Stop running.")}
  print("icldate")
  if(is.vector(icldate) == FALSE || is.vector(icltime) == FALSE){
    message("Start error message. The arguments icldate or icltime are not vector.")
    stop("End error message. Stop running.")}
  print(icldate)
  if(iclx != as.integer(iclx) || icly != as.integer(icly) || iclt != as.integer(iclt)
     || iclvr != as.integer(iclvr) || icldate[1] != as.integer(icldate[1]) ||
     icltime[1] != as.integer(icltime[1]) || iclsp != as.integer(iclsp) 
     || icldate[2] != as.integer(icldate[2]) || icltime[2] != as.integer(icltime[2])
     || icldate[3] != as.integer(icldate[3]) || icltime[3] != as.integer(icltime[3])
     || icldate[4] != as.integer(icldate[4]) || icltime[4] != as.integer(icltime[4])){
    iclx <- as.integer(iclx)
    icly <- as.integer(icly)
    iclt <- as.integer(iclt)
    icldate <- as.integer(icldate)
    icltime <- as.integer(icltime)
    iclvr <- as.integer(iclvr)
    iclsp <- as.integer(iclsp)
    message("Warning message: the arguments iclsp, iclx, icly, iclvr are forced to be integer numbers.")
  }
  if(iclx == 0 || icly == 0 || iclvr == 0){
    message("Start error message. The columns iclx, icly, iclvr have to be defined.")
    stop("End error message. Stop running.")
  }
  if (save.as != "data.frame" && save.as != "STFDF"){
    message("Start error message. The class of data has to be data.frame or STFDF. The class entered is not allowed.")
    stop("End error message. Stop running.")
  }
  if(is.vector(date.format) == FALSE){
    message("Start error message. The argument date.format is not a vector.")
    stop("End error message. Stop running.")
  }
  iflagt <- 0
  if(date.format[1] == "code"){
    if(is.na(date.format[2]) == FALSE){
      message("Start error message. The argument date.format[2] is not valid: it has to be set equal to NA.")
      stop("End error message. Stop running.")
    }
    date_format <- 0
  }else{
    if(date.format[1] == "year"){
      date_format <- 0
      iflagt <- 1
    }else{
      if(date.format[1] == "yearmon"){
        date_format <- date.format[2]
      }else{
        if(date.format[1] == "yearqtr"){
          date_format <- date.format[2]
        }else{
          if(date.format[1] == "Date"){
            date_format <- date.format[2]
          }else{
            if(date.format[1] == "POSIX"){
              date_format <- date.format[2]
            }else{
              message("Start error message. The argument date.format is not admissible.")
              stop("End error message. Stop running.")
            }}}}}}
  if(bytime != "%Y" && bytime != "%y" && bytime != "%m" && bytime != "%q" &&
     bytime != "%d" && bytime != "%H" &&
     bytime != "%M" && bytime != "%S"  && is.na(bytime) != TRUE){
    message("Start error message. The parameter in bytime is not in the required format.")
    stop("End error message. Stop running.")}
  if(date_format == 0 && iclt == 0 && iflagt == 0){
    message("Start error message. The arguments iclt and date.format are not consistent: please set one of them.")
    stop("End error message. Stop running.")}
  if(date_format != 0 && iclt != 0){ # only one column has to be given because of the missing values inclusion
    message("Start error message. The arguments iclt and date.format are not consistent: please set only one of them.")
    stop("End error message. Stop running.")}
  if(date_format != 0 && length(setdiff(icldate, 0)) == 0){
    message("Start error message. The arguments icldate and date.format are not consistent: no column has been indicated for dates.")
    stop("End error message. Stop running.")}
  if(date_format != 0 && icldate[1] != 0 && length(setdiff(icldate[2:4], 0)) > 0){
    message("Start error message. The arguments icldate and date.format are not consistent: too many columns have been indicated for dates.")
    stop("End error message. Stop running.")}
  if(date.format[1] != "year" && date.format[1] != "yearmon" &&
     date.format[1] != "yearqtr"){
  if(date_format != 0 && icldate[1] == 0 && length(setdiff(icldate[2:4], 0)) < 3){
    message("Start error message. The arguments icldate and date.format are not consistent: no information about columns for year or month or day.")
    stop("End error message. Stop running.")}}else{
    if(date.format[1] == "yearmon" || date.format[1] == "yearqtr"){
      if(icldate[1] == 0 && length(setdiff(icldate[2:3], 0)) < 2){
        message("Start error message. The arguments icldate and date.format are not consistent: no information about columns for year or month/quarter.")
        stop("End error message. Stop running.")}}
      if(date.format[1] == "year"){
        if(icldate[1] == 0 && icldate[2] == 0){
          message("Start error message. The arguments icldate and date.format are not consistent: no information about column for year.")
          stop("End error message. Stop running.")}
        if(date.format[2] != "%Y" && date.format[2] != "%y"){
          message("Start error message. The arguments date.format[1] and date.format[2] are not consistent.")
          stop("End error message. Stop running.")}
        if(is.na(bytime) == FALSE){
          message("Start error message. The arguments bytime is not consistent. Please use the default choice NA.")
          stop("End error message. Stop running.")}
        date.format[1] <- "code"
        if(icldate[1] != 0){
        iclt <- icldate[1]
        }else{iclt <- icldate[2]
        icldate[2] <- 0}
    }
  }
  if(date.format[1] == "yearmon" || date.format[1] == "yearqtr"){
    if(icldate[1] == 0 && length(setdiff(icldate[2:4], 0)) == 3){
      icldate[4] <- 0
      message("Warning message: the argument in icldate[4] is forced to be zero for classes yearmon or yearqtr.")}
    if(icltime[1] != 0 || length(setdiff(icltime[2:4], 0)) != 0){
      icltime <- 0
      message("Warning message: the arguments in icltime are forced to be zero for classes yearmon or yearqtr.")}
  }else{
  if(date_format != 0 && icltime[1] != 0 && length(setdiff(icltime[2:4], 0)) > 0){
    message("Start error message. The arguments icltime and date.format are not consistent: too many columns have been indicated for hours.")
    stop("End error message. Stop running.")}
  if(regexpr('%H', date_format)[1] > 0  && icltime[1] == 0 && icltime[2] == 0){
    message("Start error message. The argument date.format includes hours but no column in icltime provides information about hour.")
    stop("End error message. Stop running.")}
  if(regexpr('%H', date_format)[1] > 0 && regexpr('%M', date_format)[1] > 0 &&
     icltime[1] == 0 && (icltime[2] == 0 || icltime[3] == 0)){
    message("Start error message. The argument date.format includes hours and minutes but some information about hours or minutes are missing in icltime.")
    stop("End error message. Stop running.")}
  if((regexpr('%H', date_format)[1] > 0 && regexpr('%M', date_format)[1] > 0 &&
     regexpr('%S', date_format)[1] > 0) &&
     icltime[1] == 0 && (icltime[2] == 0 || icltime[3] == 0 || icltime[4] == 0)){
    message("Start error message. The argument date.format includes hours, minutes and seconds but some information are missing in icltime.")
    stop("End error message. Stop running.")}}
  if(is.scalar(tlag) == FALSE || tlag <= 0){
    message("Start error message. The argument tlag has to be integer and greater than zero.")
    stop("End error message. Stop running.")}
  iclw1 <- c(iclx,icly,iclt,iclvr,iclsp)
  iclw1 <- iclw1[iclw1 != 0]
  iclw2 <- icldate[icldate != 0]
  iclw3 <- icltime[icltime != 0]
  if(sum(match(duplicated(c(iclw1,iclw2,iclw3)),TRUE,nomatch = 0)) > 0){
    message("Start error message. The columns (iclx, icly, ...) indicated have to be different.")
    stop("End error message. Stop running.")
  }
    if(is.scalar(tlag) == FALSE){
    tlag <- as.integer(tlag)
    message("Warning message: the argument tlag is forced to be integer.")
    }
  if(is.na(bytime) == FALSE){
  if(bytime == "%S"){
    if(regexpr('%S', date_format)[1] < 0){
      message("Start error message. The arguments bytime = %S and date.format are not consistent.")
      stop("End error message. Stop running.")}
  }else{
    if(bytime == "%M"){
      if(regexpr('%M', date_format)[1] < 0){
        message("Start error message. The arguments bytime = %M and date.format are not consistent.")
        stop("End error message. Stop running.")
      }}else{
      if(bytime == "%H"){
        if(regexpr('%H', date_format)[1] < 0){
          message("Start error message. The arguments bytime = %H and date.format are not consistent.")
          stop("End error message. Stop running.")
        }}else{
        if(bytime == "%d"){
          if(regexpr('%d', date_format)[1] < 0){
            message("Start error message. The arguments bytime = %d and date.format are not consistent.")
            stop("End error message. Stop running.")
          }}else{
          if(bytime == "%m"){
            if((regexpr('%m', date_format)[1] < 0)){
              message("Start error message. The arguments bytime = %m and date.format are not consistent.")
              stop("End error message. Stop running.")
            }}else{
            if(bytime == "%q"){
              if((regexpr('%q', date_format)[1] < 0)){
                message("Start error message. The arguments bytime = %q and date.format are not consistent.")
                stop("End error message. Stop running.")
              }}else{
              if(bytime == "%Y" || bytime == "%y"){
                if((regexpr('%Y', date_format)[1] < 0) && (regexpr('%y', date_format)[1] < 0)){
                  message("Start error message. The arguments bytime = %q and date.format are not consistent.")
                  stop("End error message. Stop running.")
                }}else{
                message("Start error message. The argument bytime is not admissible.")
                stop("End error message. Stop running.")
              }}}}}}}}
   nskip <- 0
   if(header == TRUE){
       nskip <- 1
   }
  importFl <- utils::read.table(file, dec = dec, sep = sep, skip = nskip, na.strings = missing.v)
  importFl[importFl == missing.v] <- NA
  ndata <- dim(importFl)[1]
  if(date_format != 0){
    i <- 1
    while(i <= ndata){
     if(icldate[1] == 0){
        if(icldate[2] != 0){
          i.date <-  importFl[i,icldate[2]]
          if(icldate[3] != 0){
           i.date <- paste(i.date,importFl[i,icldate[3]], sep=" ")
           if(icldate[4] != 0){
             i.date <- paste(i.date,importFl[i,icldate[4]], sep=" ")
           }
           }else{
            if(icldate[4] != 0){
              message("Start error message. The date cannot be defined: it is specified the column for the day but not the column for the month.")
              stop("End error message. Stop running.")}}}else{
          message("Start error message. The date cannot be defined: the column for the year is not specified.")
          stop("End error message. Stop running.")}
     }else{i.date <- importFl[i,icldate[1]]}
      if(length(setdiff(icltime, 0)) != 0){ # not all zeros
        if(icltime[1] == 0){
        if(icltime[2] != 0){
          i.time <-  importFl[i,icltime[2]]
          if(icltime[3] != 0){
            i.time <- paste(i.time,importFl[i,icltime[3]], sep=" ")
            if(icltime[4] != 0){
              i.time <- paste(i.time,importFl[i,icltime[4]], sep=" ")
            }
          }else{
            if(icltime[4] != 0){
              message("Start error message. The time cannot be defined: it is specified the column for the second but not the column for the minute.")
              stop("End error message. Stop running.")}}
        }else{
          message("Start error message. The time cannot be defined: the column for the hour is not specified.")
          stop("End error message. Stop running.")}
          }else{i.time <-  importFl[i,icltime[1]]}
        importFl[i, "icldate"] <- paste(i.date, i.time)
      }else{importFl[i, "icldate"] <- i.date}
      if(date.format[1] != "yearmon" && date.format[1] != "yearqtr"){
      if(is.na(as.POSIXlt((importFl[i,"icldate"]), date_format, tz = time.zone)) == TRUE){
        message("Start error message. One of the possible causes is that:")
        message("- the dates in the file are not in the format specified in date.format or")
        message("- the data format is not valid or")
        message("- the dates in the file are not consistent with the given timezone (check for errors in dates involved in daylight saving time).")
        stop("End error message. Stop running.")}
        }else{
        if(date.format[1] == "yearmon"){
          if(is.na(as.yearmon(as.character(importFl[i,"icldate"]), date_format)) == TRUE){
          message("Start error message. The dates in the file are not in the format specified in date.format or the data format is not valid.")
          message("It might be hepful to provide dates with numeric elements.")
          stop("End error message. Stop running.")}}
        if(date.format[1] == "yearqtr"){
          if(is.na(as.yearqtr(as.character(importFl[i,"icldate"]), date_format)) == TRUE){
            message("Start error message. The dates in the file are not in the format specified in date.format or the data format is not valid.")
            message("It might be hepful to provide dates with numeric elements.")
            stop("End error message. Stop running.")}}
        }
      i <- i+1
    }
    if(date.format[1] != "yearmon" && date.format[1] != "yearqtr"){
          i.date <- as.POSIXlt((importFl[,"icldate"]), format = date_format, tz = time.zone)
          if(regexpr('%Y', date_format)[1] < 0){
            date_format <- "%y-%m-%d %H:%M:%S"}else{
              date_format <- "%Y-%m-%d %H:%M:%S"}
    }else{
      if(date.format[1] == "yearmon"){
      i.date <- as.POSIXlt(as.yearmon(as.character(importFl[,"icldate"]), format = date_format))
      }else{i.date <- as.POSIXlt(as.yearqtr(as.character(importFl[,"icldate"]), format = date_format))}
    }
    if(bytime == "%m"){
      i.date <- as.POSIXlt(as.yearmon(i.date))
      date.format[1] <- "yearmon"
      if(regexpr('%Y', date_format)[1] < 0){
        date_format <- "%y-%m-%d"}else{
          date_format <- "%Y-%m-%d"}
    }
    if(bytime == "%q"){
      i.date <- as.POSIXlt(as.yearqtr(i.date))
      date.format[1] <- "yearqtr"
      if(regexpr('%Y', date_format)[1] < 0){
        date_format <- "%y-%m-%d"}else{
          date_format <- "%Y-%m-%d"}
    }
    if(bytime == "%d"){
      i.date <- as.POSIXlt(as.Date(i.date))
      if(regexpr('%Y', date_format)[1] < 0){
        date_format <- "%y-%m-%d"}else{
          date_format <- "%Y-%m-%d"}
    }
    if(bytime == "%H"){
      minute(i.date) <- 0
      second(i.date) <- 0
    }
    if(bytime == "%M"){
      second(i.date) <- 0
    }
    if(bytime == "%Y" || bytime == "%y"){
      i.date <- year(i.date)
      date.format[1] <- "code"
      date_format <- 0
      bytime <- NA
      if(regexpr('%Y', date.format)[2] < 0){
        date.format[2] <- "%y"}else{
          date.format[2] <- "%Y"}
      iclt <-  ncol(importFl)
      iflagt <- 1
    }
          importFl <- importFl[,-c(ncol(importFl))]
          importFl <- cbind(importFl,i.date)
          names(importFl)[names(importFl) == "i.date"] <- "iclt"
          iclt <- ncol(importFl)
  }else{
    if(is.numeric(importFl[1,iclt]) == FALSE){
      message("Start error message. The temporal codes in the file are not numeric.")
      stop("End error message. Stop running.")} # check just on the first element since it is a data.frame
    names(importFl)[iclt] <- "iclt"
  }

  #=========================================#
  #=       IMPORT DATA FILE                =#
  #=========================================#
  importFl1 <- importFl[order(importFl[, iclx], importFl[, icly], importFl[, iclt]), ]
  print("importFl1")
  print(importFl1)
    code.time <- unique(importFl1[, iclt])
    tpar1 <- min(code.time)
    tpar2 <- max(code.time)
    vec.date <- c(tpar1)
    if(date.format[1] == "Date" || date.format[1] == "POSIX"){
      if(bytime == "%d"){
      tpar1 <- as.POSIXlt(as.Date(tpar1))
      tpar2 <- as.POSIXlt(as.Date(tpar2))
      vec.date.w <- as.POSIXlt(as.Date(tpar1))
      vec.date <- as.character(vec.date)
      }else{
        tpar1 <- as.POSIXlt(tpar1)
        tpar2 <- as.POSIXlt(tpar2)
        vec.date.w <- as.POSIXlt(tpar1)
        vec.date <- format(round(as.POSIXlt(vec.date, format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
        }
    }else{
      vec.date <- as.character(vec.date)
      if(date.format[1] == "yearmon"){
        tpar1 <- as.POSIXlt(as.yearmon(tpar1))
        tpar2 <- as.POSIXlt(as.yearmon(tpar2))
        vec.date.w <- as.POSIXlt(as.yearmon(tpar1))
      }
      if(date.format[1] == "yearqtr"){
        tpar1 <- as.POSIXlt(as.yearqtr(tpar1))
        tpar2 <- as.POSIXlt(as.yearqtr(tpar2))
        vec.date.w <- as.POSIXlt(as.yearqtr(tpar1))
      }
      if(date.format[1] == "code"){
        vec.date.w <- tpar1
      }
    }
    if(length(unique(c(tpar1,tpar2))) ==1){
      message("Start error message. The initial observed date or temporal point is coincident with the final one, thus the data structure is not spatio-temporal.")
      stop("End error message. Stop running.")
    }
    if(date_format != 0){
      if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
        if(bytime == "%d"){delta.time <- difftime(tpar2, tpar1, units = "days")}
        if(bytime == "%H"){delta.time <- difftime(tpar2, tpar1, units = "hours")}
        if(bytime == "%M"){delta.time <- difftime(tpar2, tpar1, units = "mins")}
        if(bytime == "%S"){delta.time <- difftime(tpar2, tpar1, units = "secs")}
      }else{
        if(bytime == "%m"){
          delta.time <- length(seq(from=as.Date(tpar1), to=as.Date(tpar2), by='month')) - 1
        }else{
          if(bytime == "%q"){
            delta.time <- length(seq(from=as.Date(tpar1), to=as.Date(tpar2), by='quarter')) - 1
          }}}}else{delta.time <- tpar2-tpar1}
    if(tlag>(delta.time)){
      message("Start error message. The argument tlag is not consistent: it has to be less than the temporal observed interval.")
      stop("End error message. Stop running.")}

    i <- 1
#     vec.date.w <- vec.date[1]
#     vec.date <- as.character(vec.date)
    while(delta.time >= tlag){
      i <- i + 1
      if(date_format != 0){
        if(bytime == "%Y" || bytime == "%y"){
          year(vec.date.w) <- year(vec.date.w) + tlag
        }
        if(bytime == "%d"){
          day(vec.date.w) <- day(vec.date.w) + tlag
        }
        if(bytime == "%m"){
          month(vec.date.w) <- month(vec.date.w) + tlag
        }
        if(bytime == "%q"){
          month(vec.date.w) <- month(vec.date.w) + tlag*3
        }
        if(bytime == "%H"){
          delta.w <- vec.date.w
          hour(delta.w) <- hour(vec.date.w) + tlag
          if(is.na(delta.w) == TRUE){
            vec.date.w1 <-vec.date.w
            hour(vec.date.w1) <- hour(vec.date.w1) + tlag + 1
            if(difftime(vec.date.w1, vec.date.w, units = "hours") < tlag){
              hour(vec.date.w1) <- hour(vec.date.w1) + 1
            }
            hour(vec.date.w) <- hour(vec.date.w1)
          }else{
            hour(vec.date.w) <- hour(vec.date.w) + tlag}
        }
        if(bytime == "%M"){
          delta.w <- vec.date.w
          minute(delta.w) <- minute(vec.date.w) + tlag
          if(is.na(delta.w) == TRUE){
            vec.date.w1 <-vec.date.w
            minute(vec.date.w1) <- minute(vec.date.w1) + tlag + 60
            if(difftime(vec.date.w1, vec.date.w, units = "mins") < tlag){
              minute(vec.date.w1) <- minute(vec.date.w1) + 60
            }
            minute(vec.date.w) <- minute(vec.date.w1)
            }else{
            minute(vec.date.w) <- minute(vec.date.w) + tlag}
        }
        if(bytime == "%S"){
          delta.w <- vec.date.w
          second(delta.w) <- second(vec.date.w) + tlag
          if(is.na(delta.w) == TRUE){
            vec.date.w1 <-vec.date.w
            second(vec.date.w1) <- second(vec.date.w1) + tlag + 3600
            if(difftime(vec.date.w1, vec.date.w, units = "secs") < tlag){
              second(vec.date.w1) <- second(vec.date.w1) + 3600
            }
            second(vec.date.w) <- second(vec.date.w1)
            }else{
            second(vec.date.w) <- second(vec.date.w) + tlag}
        }
        if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
          if(bytime == "%d"){delta.time <- difftime(tpar2, vec.date.w, units = "days")}
          if(bytime == "%H"){delta.time <- difftime(tpar2, vec.date.w, units = "hours")}
          if(bytime == "%M"){delta.time <- difftime(tpar2, vec.date.w, units = "mins")}
          if(bytime == "%S"){delta.time <- difftime(tpar2, vec.date.w, units = "secs")}
        }else{
          if(bytime == "%m"){
            delta.time <- length(seq(from=as.Date(vec.date.w), to=as.Date(tpar2), by='month')) - 1
          }else{
            if(bytime == "%q"){
              delta.time <- length(seq(from=as.Date(vec.date.w), to=as.Date(tpar2), by='quarter')) - 1
            }}}
      }else{
        vec.date.w <- vec.date.w + tlag
        delta.time <- tpar2-vec.date.w
      }
      vec_date <- vec.date.w
      if(date.format[1] == "POSIX" && bytime != "%d"){
        if(hour(vec_date) == 0 && minute(vec_date) == 0 && second(vec_date) == 0){
      vec_date <- format(round(as.POSIXlt(vec.date.w, format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
        }}
      vec.date <- rbind(vec.date,as.character(vec_date))
    }
    n.time <- length(unique(importFl1[, iclt]))
    for(i in 1:n.time){
      if(date.format[1] == "POSIX" && bytime != "%d"){
        if(hour(importFl1[i, iclt]) == 0 && minute(importFl1[i, iclt]) == 0 && second(importFl1[i, iclt]) == 0){
          importFl1[i, iclt] <- format(round(as.POSIXlt(importFl1[i, iclt], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
        }}

    }
    importFl1[[iclt]] <- as.character(importFl1[,iclt])
    importFl1 <- subset(importFl1, iclt %in% vec.date)
    if(date_format != 0){
    if(date.format[1] == "Date" || date.format[1] == "POSIX"){
      if(bytime == "%d"){
        importFl1[[iclt]] <- as.POSIXlt(as.Date(importFl1[,iclt], format = date_format))
      }else{
        importFl1[[iclt]] <- (as.POSIXlt(importFl1[,iclt], format = date_format, tz = time.zone))
        }
    }else{
      if(date.format[1] == "yearmon"){
        importFl1[[iclt]] <- as.POSIXlt(as.yearmon((importFl1[,iclt]), format = date_format))
      }
      if(date.format[1] == "yearqtr"){
        importFl1[[iclt]] <- as.POSIXlt(as.yearqtr((importFl1[,iclt]), format = date_format))
      }
    }
    }else{importFl1[[iclt]]<-as.numeric(importFl1[,iclt])}
    sp <- cbind(importFl1[, iclx], importFl1[, icly])
    sp <- unique(sp)
    n.stat <- nrow(sp)
    n.time <- length(unique(importFl1[, iclt]))
    code.time <- unique(importFl1[, iclt])
    tpar2 <- max(code.time)
    if(date_format != 0){
    if(date.format[1] == "Date" || date.format[1] == "POSIX"){
      if(bytime == "%d"){
        tpar2 <- as.POSIXlt(as.Date(tpar2))
      }else{
        tpar2 <- (as.POSIXlt(tpar2))}
    }else{
      if(date.format[1] == "yearmon"){
        tpar2 <- as.POSIXlt(as.yearmon(tpar2))
      }
      if(date.format[1] == "yearqtr"){
        tpar2 <- as.POSIXlt(as.yearqtr(tpar2))
      }
    }}else{tpar2 <- as.numeric(tpar2)}
    if(length(unique(c(tpar1,tpar2))) == 1){
      message("Start error message. There is no temporal point for the selected tlag increment.")
      message("Thus the data structure is not spatio-temporal.")
      stop("End error message. Stop running.")
    }
    for(i in 1:n.stat){
      print("i")
      print(i)
      if(date_format != 0){
      if(sum(match(duplicated(importFl1[(importFl1[iclx] == sp[i,1] &
                              importFl1[icly] == sp[i,2]),iclt]), TRUE, nomatch = 0)) >= 1){
        message("Start error message. According to the argument 'bytime = '", bytime,", there are more than one observation for a fixed spatial point and a fixed time point.")
        stop("End error message. Stop running.")}
      }else{
        if(sum(match(duplicated(importFl1[(importFl1[iclx] == sp[i,1] &
                                           importFl1[icly] == sp[i,2]),iclt]), TRUE, nomatch = 0)) >= 1){
          message("Start error message. There are more than one observation for a fixed spatial point and a fixed time point.")
          stop("End error message. Stop running.")}}}
    if(iclsp == 0){
      importFl1 <- importFl1[,c(iclx,icly,iclt,iclvr)]
      colnames(importFl1) <- c("iclx","icly","iclt","iclvr")
      importFlw <- data.frame(iclx=NA,icly=NA,iclt=NA,iclvr=NA)
      importFlfull <- importFl1[1,]
      colnames(importFlfull) <- c("iclx","icly","iclt", "iclvr")
    }else{
      importFl1 <- importFl1[,c(iclsp,iclx,icly,iclt,iclvr)]
      colnames(importFl1) <- c("iclsp","iclx","icly","iclt","iclvr")
      importFlw <- data.frame(iclsp=NA,iclx=NA,icly=NA,iclt=NA,iclvr=NA)
      importFlfull <- importFl1[1,]
      colnames(importFlfull) <- c("iclsp","iclx","icly","iclt", "iclvr")
      iclsp <- "iclsp"
    }
    iclx <- "iclx"
    icly <- "icly"
    iclvr <- "iclvr"
    iclt <- "iclt"
    if(date.format[1] == "POSIX" && bytime != "%d"){
      if(hour(importFlfull[[iclt]]) == 0 && minute(importFlfull[[iclt]]) == 0 && second(importFlfull[[iclt]]) == 0){
        importFlfull[[iclt]] <- format(round(as.POSIXlt(importFlfull[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
      }}
    importFlfull[[iclt]] <- as.character(importFlfull[[iclt]])
    ### CHECK ON MISSING DATES ###
    t1 <- importFl1[1,iclt]
    if(date_format != 0){
    if(date.format[1] == "Date" || date.format[1] == "POSIX"){
      if(bytime == "%d"){
        t1 <- as.POSIXlt(as.Date(t1))
      }else{
      t1 <- as.POSIXlt(t1)}
    }else{
      if(date.format[1] == "yearmon"){
        t1 <- as.POSIXlt(as.yearmon(t1))
      }
      if(date.format[1] == "yearqtr"){
        t1 <- as.POSIXlt(as.yearqtr(t1))
      }
      }
    }else{t1 <- as.numeric(t1)}
    importFlw <- importFl1[1,]
    j <- 1
    if(date_format != 0){
      if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
        if(bytime == "%d"){
          delta.time <- difftime(t1,tpar1, units = "days")
        }
        if(bytime == "%H"){
          delta.time <- difftime(t1,tpar1, units = "hours")
        }
        if(bytime == "%M"){
          delta.time <- difftime(t1,tpar1, units = "mins")
        }
        if(bytime == "%S"){
          delta.time <- difftime(t1,tpar1, units = "secs")
        }
      }else{
        if(bytime == "%m"){
          delta.time <- length(seq(from=as.Date(tpar1), to=as.Date(t1), by='month')) - 1
        }else{
          if(bytime == "%q"){
            delta.time <- length(seq(from=as.Date(tpar1), to=as.Date(t1), by='quarter')) - 1
          }
        }
      }
    }else{
      delta.time <- t1-tpar1}
    while(delta.time >= tlag){
      j <- j + 1
      if(as.integer(j/500) == (j/500) ){
        message("The first ", j, " values have been included in the data file.")
      }
      if(date_format != 0){
        if(bytime == "%Y" || bytime == "%y"){
          year(t1) <- year(t1) - tlag
        }
        if(bytime == "%d"){
          day(t1) <- day(t1) - tlag
        }
        if(bytime == "%m"){
          month(t1) <- month(t1) - tlag
        }
        if(bytime == "%q"){
          month(t1) <- month(t1) - tlag*3
        }
        if(bytime == "%H"){
          delta.w <- t1
          hour(delta.w) <- hour(t1) - tlag
          if(is.na(delta.w) == TRUE){
            t11 <- t1
            hour(t11) <- hour(t11) - tlag - 1
            if(difftime(t1, t11, units = "hours") < tlag){
              hour(t11) <- hour(t11) - 1
            }
            hour(t1) <- hour(t11)
          }else{
          hour(t1) <- hour(t1) - tlag}
        }
        if(bytime == "%M"){
          delta.w <- t1
          minute(delta.w) <- minute(t1) - tlag
          if(is.na(delta.w) == TRUE){
            t11 <- t1
            minute(t11) <- minute(t11) - tlag - 60
            if(difftime(t1, t11, units = "mins") < tlag){
              minute(t11) <- minute(t11) - 60
            }
            minute(t1) <- minute(t11)
          }else{
            minute(t1) <- minute(t1) - tlag}
        }
        if(bytime == "%S"){
          delta.w <- t1
          second(delta.w) <- second(t1) - tlag
          if(is.na(delta.w) == TRUE){
            t11 <- t1
            second(t11) <- second(t11) - tlag - 3600
            if(difftime(t1, t11, units = "secs") < tlag){
              second(t11) <- second(t11) - 3600
            }
            second(t1) <- second(t11)
          }else{
            second(t1) <- second(t1) - tlag}
        }
        if(date.format[1] == "Date" || date.format[1] == "POSIX"){
          if(bytime == "%d"){
            importFlw[[iclt]]<- as.POSIXlt(as.Date(t1))
          }else{
            importFlw[[iclt]]<- (as.POSIXlt(t1))}
        }else{
          if(date.format[1] == "yearmon"){
            importFlw[[iclt]]<- as.POSIXlt(as.yearmon(t1))
          }
          if(date.format[1] == "yearqtr"){
            importFlw[[iclt]]<- as.POSIXlt(as.yearqtr(t1))
          }
        }
      }else{
        t1 <- t1 - tlag
       importFlw[[iclt]]<- importFlw[1,iclt] - tlag
      }
      importFlw[[iclvr]]<- NA
      if(date.format[1] == "POSIX" && bytime != "%d"){
        if(hour(importFlw[[iclt]]) == 0 && minute(importFlw[[iclt]]) == 0 && second(importFlw[[iclt]]) == 0){
          importFlw[[iclt]] <- format(round(as.POSIXlt(importFlw[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
        }}
      importFlw[[iclt]] <- as.character(importFlw[[iclt]])
      importFlfull <- rbind(importFlw,importFlfull)
      if(date_format != 0){
        if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
          if(bytime == "%d"){delta.time <- difftime(t1,tpar1, units = "days")}
          if(bytime == "%H"){delta.time <- difftime(t1,tpar1, units = "hours")}
          if(bytime == "%M"){delta.time <- difftime(t1,tpar1, units = "mins")}
          if(bytime == "%S"){delta.time <- difftime(t1,tpar1, units = "secs")}
        }else{
          if(bytime == "%m"){
            delta.time <- length(seq(from=as.Date(tpar1), to=as.Date(t1), by='month')) - 1
          }else{
            if(bytime == "%q"){
              delta.time <- length(seq(from=as.Date(tpar1), to=as.Date(t1), by='quarter')) - 1
            }
          }
        }
      }else{delta.time <- t1-tpar1}
    }
    i <- 2
    ndata <- dim(importFl1)[1]
    print("ndata")
    print(ndata)
    while (i<= ndata) {
      print("i second cicle")
      print(i)
     if((importFl1[i-1,iclx] == importFl1[i,iclx]) && (importFl1[i-1,icly] == importFl1[i,icly])){
       t1 <- importFl1[i-1,iclt]
       t2 <- importFl1[i,iclt]
       if(date_format != 0){
         if(date.format[1] == "Date" || date.format[1] == "POSIX"){
           if(bytime == "%d"){
             t1 <- as.POSIXlt(as.Date(t1))
             t2 <- as.POSIXlt(as.Date(t2))
           }else{
           t1 <- as.POSIXlt(t1)
           t2 <- as.POSIXlt(t2)
         }
         }else{
           if(date.format[1] == "yearmon"){
             t1 <- as.POSIXlt(as.yearmon(t1))
             t2 <- as.POSIXlt(as.yearmon(t2))
           }
           if(date.format[1] == "yearqtr"){
             t1 <- as.POSIXlt(as.yearqtr(t1))
             t2 <- as.POSIXlt(as.yearqtr(t2))
           }
         }
       }else{
         t1 <- as.numeric(t1)
         t2 <- as.numeric(t2)}
       if(date_format != 0){
         if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
           if(bytime == "%d"){
             delta.time <- difftime(t2, t1, units = "days")
           }
           if(bytime == "%H"){
             delta.time <- difftime(t2, t1, units = "hours")
           }
           if(bytime == "%M"){
             delta.time <- difftime(t2, t1, units = "mins")
           }
           if(bytime == "%S"){
             delta.time <- difftime(t2, t1, units = "secs")
           }
         }else{
           if(bytime == "%m"){
             delta.time <- length(seq(from=as.Date(t1), to=as.Date(t2), by='month')) - 1
           }else{
             if(bytime == "%q"){
               delta.time <- length(seq(from=as.Date(t1), to=as.Date(t2), by='quarter')) - 1
             }
           }
         }
       }else{delta.time <- t2-t1}
       while(delta.time > tlag){
         j <- j + 1
         if(as.integer(j/500) == (j/500) ){
           message("The first ", j, " values have been included in the data file.")
         }
         importFlw<-importFl1[i-1,]
         if(date_format != 0){
            if(bytime == "%Y" || bytime == "%y"){
              year(t1) <- year(t1) + tlag
            }
           if(bytime == "%d"){
             day(t1) <- day(t1) + tlag
           }
           if(bytime == "%m"){
             month(t1) <- month(t1) + tlag
           }
           if(bytime == "%q"){
             month(t1) <- month(t1) + tlag*3
           }
           if(bytime == "%H"){
             delta.w <- t1
             hour(delta.w) <- hour(t1) + tlag
             if(is.na(delta.w) == TRUE){
               t11 <- t1
               hour(t11) <- hour(t11) + tlag + 1
               if(difftime(t11, t1, units = "hours") < tlag){
                 hour(t11) <- hour(t11) + 1
               }
               hour(t1) <- hour(t11)
             }else{
               hour(t1) <- hour(t1) + tlag}
           }
           if(bytime == "%M"){
             delta.w <- t1
             minute(delta.w) <- minute(t1) + tlag
             if(is.na(delta.w) == TRUE){
               t11 <- t1
               minute(t11) <- minute(t11) + tlag + 60
               if(difftime(t11, t1, units = "mins") < tlag){
                 minute(t11) <- minute(t11) + 60
               }
               minute(t1) <- minute(t11)
             }else{
               minute(t1) <- minute(t1) + tlag}
           }
           if(bytime == "%S"){
             delta.w <- t1
             second(delta.w) <- second(t1) + tlag
             if(is.na(delta.w) == TRUE){
               t11 <- t1
               second(t11) <- second(t11) + tlag + 3600
               if(difftime(t11, t1, units = "secs") < tlag){
                 second(t11) <- second(t11) + 3600
               }
               second(t1) <- second(t11)
             }else{
               second(t1) <- second(t1) + tlag}
           }
           if(date.format[1] == "Date" || date.format[1] == "POSIX"){
             if(bytime == "%d"){
               importFlw[[iclt]]<- as.POSIXlt(as.Date(t1))
             }else{
               importFlw[[iclt]]<- as.POSIXlt(t1)}
           }else{
             if(date.format[1] == "yearmon"){
               importFlw[[iclt]]<- as.POSIXlt(as.yearmon(t1))
             }
             if(date.format[1] == "yearqtr"){
               importFlw[[iclt]]<- as.POSIXlt(as.yearqtr(t1))
             }
           }
         }else{
           t1 <- t1 + tlag
          importFlw[[iclt]]<- t1
         }
         importFlw[[iclvr]]<- NA
         if(date.format[1] == "POSIX" && bytime != "%d"){
           if(hour(importFlw[[iclt]]) == 0 && minute(importFlw[[iclt]]) == 0 && second(importFlw[[iclt]]) == 0){
             importFlw[[iclt]] <- format(round(as.POSIXlt(importFlw[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
           }}
         importFlw[[iclt]] <- as.character(importFlw[[iclt]])
         importFlfull <- rbind(importFlfull, importFlw)
         if(date_format != 0){
           if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
             if(bytime == "%d"){delta.time <- difftime(t2, t1, units = "days")}
             if(bytime == "%H"){delta.time <- difftime(t2, t1, units = "hours")}
             if(bytime == "%M"){delta.time <- difftime(t2, t1, units = "mins")}
             if(bytime == "%S"){delta.time <- difftime(t2, t1, units = "secs")}
           }else{
             if(bytime == "%m"){
               delta.time <- length(seq(from=as.Date(t1), to=as.Date(t2), by='month')) - 1
             }else{
               if(bytime == "%q"){
                 delta.time <- length(seq(from=as.Date(t1), to=as.Date(t2), by='quarter')) - 1
               }
             }
           }
         }else{delta.time <- t2 - t1}
       }
       importFlw <- importFl1[i,]
       if(date.format[1] == "POSIX" && bytime != "%d"){
         if(hour(importFlw[[iclt]]) == 0 && minute(importFlw[[iclt]]) == 0 && second(importFlw[[iclt]]) == 0){
           importFlw[[iclt]] <- format(round(as.POSIXlt(importFlw[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
         }}
       importFlw[[iclt]] <- as.character(importFlw[1, iclt])
       importFlfull <- rbind(importFlfull, importFlw)
         i <- i + 1
         j <- j +1
     }else{
       if(importFl1[i-1,iclt] != tpar2){
         t2 <- importFl1[i-1,iclt]
         if(date_format != 0){
           if(date.format[1] == "Date" || date.format[1] == "POSIX"){
             if(bytime == "%d"){
               t2 <- as.POSIXlt(as.Date(t2))
             }else{
           t2 <- as.POSIXlt(t2)}
           }else{
             if(date.format[1] == "yearmon"){
               t2 <- as.POSIXlt(as.yearmon(t2))
             }
             if(date.format[1] == "yearqtr"){
               t2 <- as.POSIXlt(as.yearqtr(t2))
             }
           }
         }else{t2 <- as.numeric(t2)}
         if(date_format != 0){
           if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
             if(bytime == "%d"){delta.time <- difftime(tpar2, t2, units = "days")}
             if(bytime == "%H"){delta.time <- difftime(tpar2, t2, units = "hours")}
             if(bytime == "%M"){delta.time <- difftime(tpar2, t2, units = "mins")}
             if(bytime == "%S"){delta.time <- difftime(tpar2, t2, units = "secs")}
           }else{
             if(bytime == "%m"){
               delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='month')) - 1
             }else{
               if(bytime == "%q"){
                 delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='quarter')) - 1
               }
             }
           }
         }else{delta.time <- tpar2-t2}
         while(delta.time >= tlag){
           j <- j +1
           if(as.integer(j/500) == (j/500) ){
             message("The first ", j, " values have been included in the data file.")
           }
           importFlw <- importFl1[i - 1,]
           if(date_format != 0){
             if(bytime == "%Y" || bytime == "%y"){
               year(t2) <- year(t2) + tlag
             }
             if(bytime == "%d"){
               day(t2) <- day(t2) + tlag
             }
             if(bytime == "%m"){
               month(t2) <- month(t2) + tlag
             }
             if(bytime == "%q"){
               month(t2) <- month(t2) + tlag*3
             }
             if(bytime == "%H"){
               delta.w <- t2
               hour(delta.w) <- hour(t2) + tlag
               if(is.na(delta.w) == TRUE){
                 t11 <- t2
                 hour(t11) <- hour(t11) + tlag + 1
                 if(difftime(t11, t2, units = "hours") < tlag){
                   hour(t11) <- hour(t11) + 1
                 }
                 hour(t2) <- hour(t11)
               }else{
                 hour(t2) <- hour(t2) + tlag}
             }
             if(bytime == "%M"){
               delta.w <- t2
               minute(delta.w) <- minute(t2) + tlag
               if(is.na(delta.w) == TRUE){
                 t11 <- t2
                 minute(t11) <- minute(t11) + tlag + 60
                 if(difftime(t11, t2, units = "mins") < tlag){
                   minute(t11) <- minute(t11) + 60
                 }
                 minute(t2) <- minute(t11)
               }else{
                   minute(t2) <- minute(t2) + tlag}
             }
             if(bytime == "%S"){
               delta.w <- t2
               second(delta.w) <- second(t2) + tlag
               if(is.na(delta.w) == TRUE){
                 t11 <- t2
                 second(t11) <- second(t11) + tlag + 3600
                 if(difftime(t11, t2, units = "secs") < tlag){
                   second(t11) <- second(t11) + 3600
                 }
                 second(t2) <- second(t11)
               }else{
                   second(t2) <- second(t2) + tlag}
             }
             if(date.format[1] == "Date" || date.format[1] == "POSIX"){
               if(bytime == "%d"){
                 t2 <- as.POSIXlt(as.Date(t2))
                 importFlw[[iclt]]<- as.POSIXlt(as.Date(t2))
               }else{
                 importFlw[[iclt]]<- (as.POSIXlt(t2))
              }
             }else{
               if(date.format[1] == "yearmon"){
                 importFlw[[iclt]]<- as.POSIXlt(as.yearmon(t2))
               }
               if(date.format[1] == "yearqtr"){
                 importFlw[[iclt]]<- as.POSIXlt(as.yearqtr(t2))
               }
             }
           }else{
             t2 <- t2 + tlag
            importFlw[[iclt]] <- t2
           }
           importFlw[[iclvr]] <- NA
           if(date.format[1] == "POSIX" && bytime != "%d"){
             if(hour(importFlw[[iclt]]) == 0 && minute(importFlw[[iclt]]) == 0 && second(importFlw[[iclt]]) == 0){
               importFlw[[iclt]] <- format(round(as.POSIXlt(importFlw[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
             }}
           importFlw[[iclt]] <- as.character(importFlw[[iclt]])
           importFlfull <- rbind(importFlfull, importFlw)
           if(date_format != 0){
             if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
               if(bytime == "%d"){delta.time <- difftime(tpar2, t2, units = "days")}
               if(bytime == "%H"){delta.time <- difftime(tpar2, t2, units = "hours")}
               if(bytime == "%M"){delta.time <- difftime(tpar2, t2, units = "mins")}
               if(bytime == "%S"){delta.time <- difftime(tpar2, t2, units = "secs")}
             }else{
               if(bytime == "%m"){
                 delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='month')) - 1
               }else{
                 if(bytime == "%q"){
                   delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='quarter')) - 1
                 }
               }
             }}else{delta.time <- tpar2-t2}
         }
       }
       t1 <- importFl1[i,iclt]
       t0 <- tpar1
       if(date_format != 0){
         if(date.format[1] == "Date" || date.format[1] == "POSIX"){
           if(bytime == "%d"){
             t1 <- as.POSIXlt(as.Date(t1))
             t0 <- as.POSIXlt(as.Date(t0))
           }else{
           t0 <- as.POSIXlt(t0)
           t1 <- as.POSIXlt(t1)}
         }else{
           if(date.format[1] == "yearmon"){
             t1 <- as.POSIXlt(as.yearmon(t1))
             t0 <- as.POSIXlt(as.yearmon(t0))
           }
           if(date.format[1] == "yearqtr"){
             t1 <- as.POSIXlt(as.yearqtr(t1))
             t0 <- as.POSIXlt(as.yearqtr(t0))
           }
         }
       }else{
         t0 <- as.numeric(t0)
         t1 <- as.numeric(t1)}
       importFlw <-importFl1[i,]
       if(date_format != 0){
         if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
           if(bytime == "%d"){delta.time <- difftime(t1, t0, units = "days")}
           if(bytime == "%H"){delta.time <- difftime(t1, t0, units = "hours")}
           if(bytime == "%M"){delta.time <- difftime(t1, t0, units = "mins")}
           if(bytime == "%S"){delta.time <- difftime(t1, t0, units = "secs")}
         }else{
           if(bytime == "%m"){
             delta.time <- length(seq(from=as.Date(t0), to=as.Date(t1), by='month')) - 1
           }else{
             if(bytime == "%q"){
               delta.time <- length(seq(from=as.Date(t0), to=as.Date(t1), by='quarter')) - 1
             }
           }
         }
       }else{
         delta.time <- t1-t0
       }
       while(delta.time >= tlag){
         j <- j +1
         if(as.integer(j/500) == (j/500) ){
           message("The first ", j, " values have been included in the data file.")
         }
         if(date_format != 0){
           if(date.format[1] == "Date" || date.format[1] == "POSIX"){
             if(bytime == "%d"){
               importFlw[[iclt]]<- as.POSIXlt(as.Date(t0))
             }else{
               importFlw[[iclt]]<- (as.POSIXlt(t0))}
           }else{
             if(date.format[1] == "yearmon"){
               importFlw[[iclt]]<- as.POSIXlt(as.yearmon(t0))
             }
             if(date.format[1] == "yearqtr"){
               importFlw[[iclt]]<- as.POSIXlt(as.yearqtr(t0))
               }
           }
           if(bytime == "%Y" || bytime == "%y"){
             year(t0) <-  year(t0) + tlag
           }
           if(bytime == "%d"){
             day(t0) <-  day(t0) + tlag
           }
           if(bytime == "%m"){
             month(t0) <-  month(t0) + tlag
           }
          if(bytime == "%q"){
            month(t0) <- month(t0) + tlag*3
          }
           if(bytime == "%H"){
             delta.w <- t0
             hour(delta.w) <- hour(t0) + tlag
             if(is.na(delta.w) == TRUE){
               t11 <- t0
               hour(t11) <- hour(t11) + tlag + 1
               if(difftime(t11, t0, units = "hours") < tlag){
                 hour(t11) <- hour(t11) + 1
               }
               hour(t0) <- hour(t11)
             }else{
               hour(t0) <- hour(t0) + tlag}
           }
           if(bytime == "%M"){
             delta.w <- t0
             minute(delta.w) <- minute(t0) + tlag
             if(is.na(delta.w) == TRUE){
               t11 <- t0
               minute(t11) <- minute(t11) + tlag + 60
               if(difftime(t11, t0, units = "mins") < tlag){
                 minute(t11) <- minute(t11) + 60
               }
               minute(t0) <- minute(t11)
             }else{
               minute(t0) <- minute(t0) + tlag}
           }
           if(bytime == "%S"){
             delta.w <- t0
             second(delta.w) <- second(t0) + tlag
             if(is.na(delta.w) == TRUE){
               t11 <- t0
               second(t11) <- second(t11) + tlag + 3600
               if(difftime(t11, t0, units = "secs") < tlag){
                 second(t11) <- second(t11) + 3600
               }
               second(t0) <- second(t11)
             }else{
               second(t0) <- second(t0) + tlag}
           }
         }else{
         importFlw[[iclt]] <-  t0
           t0 <- t0 + tlag
         }
         importFlw[[iclvr]]<- NA
         if(date.format[1] == "POSIX" && bytime != "%d"){
           if(hour(importFlw[[iclt]]) == 0 && minute(importFlw[[iclt]]) == 0 && second(importFlw[[iclt]]) == 0){
             importFlw[[iclt]] <- format(round(as.POSIXlt(importFlw[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
           }}
         importFlw[[iclt]] <- as.character(importFlw[[iclt]])
         importFlfull <- rbind(importFlfull,importFlw)
         if(date_format != 0){
           if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
             if(bytime == "%d"){delta.time <- difftime(t1, t0, units = "days")}
             if(bytime == "%H"){delta.time <- difftime(t1, t0, units = "hours")}
             if(bytime == "%M"){delta.time <- difftime(t1, t0, units = "mins")}
             if(bytime == "%S"){delta.time <- difftime(t1, t0, units = "secs")}
           }else{
             if(bytime == "%m"){
               delta.time <- length(seq(from=as.Date(t0), to=as.Date(t1), by='month')) - 1
             }else{
               if(bytime == "%q"){
                 delta.time <- length(seq(from=as.Date(t0), to=as.Date(t1), by='quarter')) - 1
               }
             }
           }
         }else{
           delta.time <- t1-t0
         }
       }
       importFlw <- importFl1[i,]
       if(date.format[1] == "POSIX" && bytime != "%d"){
         if(hour(importFlw[[iclt]]) == 0 && minute(importFlw[[iclt]]) == 0 && second(importFlw[[iclt]]) == 0){
           importFlw[[iclt]] <- format(round(as.POSIXlt(importFlw[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
         }}
       importFlw[[iclt]] <- as.character(importFlw[1, iclt])
       importFlfull <- rbind(importFlfull, importFlw)
       i <- i +1
       j <- j+1
     }
    }# END CICLE OVER THE DATA
    if(importFl1[ndata,iclt] != tpar2){
       t2 <- importFl1[ndata,iclt]
       if(date_format != 0){
       if(date.format[1] == "Date" || date.format[1] == "POSIX"){
         if(bytime == "%d"){
           t2 <- as.POSIXlt(as.Date(t2))
         }else{t2 <- as.POSIXlt(t2)}
       }else{
         if(date.format[1] == "yearmon"){
           t2 <- as.POSIXlt(as.yearmon(t2))
         }
         if(date.format[1] == "yearqtr"){t2 <- as.POSIXlt(as.yearqtr(t2))}
       }}else{t2 <- as.numeric(t2)}
       importFlw <-importFl1[ndata,]
       if(date_format != 0){
         if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
           if(bytime == "%d"){delta.time <- difftime(tpar2, t2, units = "days")}
           if(bytime == "%H"){delta.time <- difftime(tpar2, t2, units = "hours")}
           if(bytime == "%M"){delta.time <- difftime(tpar2, t2, units = "mins")}
           if(bytime == "%S"){delta.time <- difftime(tpar2, t2, units = "secs")}
         }else{
           if(bytime == "%m"){
             delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='month')) - 1
           }else{
             if(bytime == "%q"){
               delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='quarter')) - 1
             }
           }
         }
       }else{
         delta.time <- tpar2-t2
       }
       while(delta.time >= tlag){
        j <- j +1
        if(as.integer(j/500) == (j/500) ){
          message("The first ", j, " values have been included in the data file.")
        }
        if(date_format != 0){
          if(bytime == "%Y" || bytime == "%y"){
            year(t2) <- year(t2) + tlag
          }
          if(bytime == "%d"){
            day(t2) <- day(t2) + tlag
          }
          if(bytime == "%m"){
            month(t2) <- month(t2) + tlag
          }
          if(bytime == "%q"){
            month(t2) <- month(t2) + tlag*3
          }
          if(bytime == "%H"){
            delta.w <- t2
            hour(delta.w) <- hour(t2) + tlag
            if(is.na(delta.w) == TRUE){
              t11 <- t2
              hour(t11) <- hour(t11) + tlag + 1
              if(difftime(t11, t2, units = "hours") < tlag){
                hour(t11) <- hour(t11) + 1
              }
              hour(t2) <- hour(t11)
            }else{
              hour(t2) <- hour(t2) + tlag}
          }
          if(bytime == "%M"){
            delta.w <- t2
            minute(delta.w) <- minute(t2) + tlag
            if(is.na(delta.w) == TRUE){
              t11 <- t2
              minute(t11) <- minute(t11) + tlag + 60
              if(difftime(t11, t2, units = "mins") < tlag){
                minute(t11) <- minute(t11) + 60
              }
              minute(t2) <- minute(t11)
            }else{
              minute(t2) <- minute(t2) + tlag}
          }
          if(bytime == "%S"){
            delta.w <- t2
            second(delta.w) <- second(t2) + tlag
            if(is.na(delta.w) == TRUE){
              t11 <- t2
              second(t11) <- second(t11) + tlag + 3600
              if(difftime(t11, t2, units = "secs") < tlag){
                second(t11) <- second(t11) + 3600
              }
              second(t2) <- second(t11)
            }else{
              second(t2) <- second(t2) + tlag}
          }
          if(date.format[1] == "Date" || date.format[1] == "POSIX"){
            if(bytime == "%d"){
              importFlw[[iclt]] <- as.POSIXlt(as.Date(t2))
            }else{
              importFlw[[iclt]] <- as.POSIXlt(t2)}
          }else{
            if(date.format[1] == "yearmon"){
              importFlw[[iclt]] <- as.POSIXlt(as.yearmon(t2))
            }
            if(date.format[1] == "yearqtr"){
              importFlw[[iclt]] <- as.POSIXlt(as.yearqtr(t2))
            }
          }
        }else{
          t2 <- t2 + tlag
         importFlw[[iclt]] <- t2
        }
        importFlw[[iclvr]] <- NA
        if(date.format[1] == "POSIX" && bytime != "%d"){
          if(hour(importFlw[[iclt]]) == 0 && minute(importFlw[[iclt]]) == 0 && second(importFlw[[iclt]]) == 0){
            importFlw[[iclt]] <- format(round(as.POSIXlt(importFlw[[iclt]], format="%Y-%m-%d %H:%M:%S")), "%Y-%m-%d %H:%M:%S")
          }}
        importFlw[[iclt]] <- as.character(importFlw[[iclt]])
        importFlfull <- rbind(importFlfull, importFlw)
        if(date_format != 0){
          if(bytime == "%d" || bytime == "%H" || bytime == "%M" || bytime == "%S"){
            if(bytime == "%d"){delta.time <- difftime(tpar2, t2, units = "days")}
            if(bytime == "%H"){delta.time <- difftime(tpar2, t2, units = "hours")}
            if(bytime == "%M"){delta.time <- difftime(tpar2, t2, units = "mins")}
            if(bytime == "%S"){delta.time <- difftime(tpar2, t2, units = "secs")}
          }else{
            if(bytime == "%m"){
              delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='month')) - 1
            }else{
              if(bytime == "%q"){
                delta.time <- length(seq(from=as.Date(t2), to=as.Date(tpar2), by='quarter')) - 1
              }
            }
          }
        }else{delta.time <- tpar2-t2}
      }
    }
    ### END CHECK ON MISSING DATES ###
    print("Begin final conversion")
    if(date_format != 0){
      if(date.format[1] == "Date" || date.format[1] == "POSIX"){
        if(bytime == "%d"){
          importFlfull[[iclt]] <- as.Date(as.character(importFlfull[[iclt]]), format = date_format)
        }else{
          importFlfull[[iclt]] <- as.POSIXlt(importFlfull[[iclt]], format = date_format, tz = time.zone)
        }}else{
        if(date.format[1] == "yearmon"){
          importFlfull[[iclt]] <- as.yearmon(importFlfull[[iclt]], format = date_format)
        }
        if(date.format[1] == "yearqtr"){
          importFlfull[[iclt]] <- as.yearqtr(importFlfull[[iclt]], format = date_format)}
        }}else{importFlfull[[iclt]] <- as.numeric(importFlfull[[iclt]])}
    print("End final conversion")
    n.time <- length(unique(importFlfull[, iclt]))
    if (iclsp == 0) {
      ID_points <- c(rep(1:n.stat, each = n.time))
      ID_points <- paste("id_", ID_points , sep = "")
      importFlfull[,"iclsp"] <- ID_points
    }
    if (save.as == "data.frame") {
      print("begin new temporal index")
      if (date_format != 0 || iflagt == 1) {
        ID_times <- c(rep(1:n.time, times = n.stat))
        importFlfull[,"timeIndex"] <- ID_times
        importFlfull <- importFlfull[,c("iclsp", "iclx","icly","timeIndex","iclt","iclvr")]
        colnames(importFlfull) <- c("spatialIndex","x","y","timeIndex","date","variable")
      }else{
        colnames(importFlfull) <- c("spatialIndex","x","y","timeIndex","variable")}
      print("end new temporal index")
    return(importFlfull)
    }

  #=========================================#
  #=       CREATE STFDF (gstat package)    =#
  #=========================================#
  if (save.as == "STFDF") {
    if(date_format == 0 && iflagt == 0){
      message("Start error message. The data set cannot be saved as STFDF since time information is missing.")
      stop("End error message. Stop running.")
    }
    #== SPATIAL DB ==#
    importFl <- importFlfull
    sp <- cbind(importFl[, iclx], importFl[, icly])
    sp <- unique(sp)
    sp.names <- unique(importFl[, iclsp])
    colnames(sp) <- c("x", "y")
    sp2 <- sp::SpatialPoints(sp)
    row.names(sp2) <- sp.names
    n.time <- length(unique(importFl[, iclt]))
    #== TEMPORAL DB ==#
    if(iflagt == 1){
     importFl[[iclt]]<-as.character(importFl[,iclt])
     data.time <- (as.Date(unique(importFl[,iclt]), format = date.format[2]))
     month(data.time) <- 1
     day(data.time) <- 1
    }else{data.time <- unique(importFl[,iclt])}
    #== ST DATA  ==#
    mydata <- importFl[, iclvr]
    mydata <- matrix(mydata, ncol = n.time, byrow = TRUE)
    #== STFDF ==#
    stfdf <- STFDF(sp2, data.time, data.frame(variable = as.vector(as.matrix(mydata))))
    return(stfdf)
  }
}
