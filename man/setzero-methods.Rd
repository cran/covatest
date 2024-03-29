% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/couples.R
\docType{methods}
\name{setzero}
\alias{setzero}
\alias{setzero,couples-method}
\alias{setzero,ANY,ANY-method}
\title{setzero}
\usage{
setzero(x, zero = TRUE, index = NULL, value)

\S4method{setzero}{couples}(x, zero = TRUE, index = NULL, value)
}
\arguments{
\item{x}{object of class \code{couples}}

\item{zero}{logical, if \code{TRUE} (the default) all negative temporal
lags are replaced with zero. If \code{x@typetest} is equal to
\code{"sym"} (symmetry test) the argument \code{setzero} is ignored because
both positive and negative temporal lags are required for computing the test}

\item{index}{two column matrix. Each row of the matrix \code{index} contains
the specific row and column, of the slot \code{@couples.st}, for which the
spatio-temporal covariance is not required}

\item{value}{numeric, the value to be replaced. Note that this method is reasonable
to be used only to replace a value equal to zero}
}
\description{
Through the function {\link{couples}}, \code{m} spatial couples and \code{n}
temporal lags are provided, hence a set of \code{m x n} spatio-temporal lags
are defined. If some of these lags are not required for the specific test, they
can be set equal to zero by using the \code{setzero} method for object of class
\code{couples}
}
\examples{
sel.staz.sym <- c("DERP016", "DENW065", "DEHE051", "DETH026", "DENW063", "DENI019",
"DENW068", "DEHE046", "DEUB029", "DEBY047", "DETH061", "DESN049")

sp.couples.in.sym <- matrix(data = c("DERP016", "DENW065", "DEHE051", "DETH026",
"DENW063", "DENI019", "DENW068", "DEHE046", "DEUB029", "DEBY047", "DETH061", "DESN049"),
ncol = 2, byrow = TRUE)

t.couples.in.sym <- c(1, 2)

couples.sym <- couples(sel.staz = sel.staz.sym, sp.couples.in = sp.couples.in.sym,
t.couples.in = t.couples.in.sym, typetest = "sym", typecode = character())

zero.index <- matrix(data=c(1,3,1,4,2,5,2,6), ncol=2, byrow = TRUE)

setzero(couples.sym, zero = FALSE, index = zero.index, value = 0)

}
\seealso{
\code{\link{couples}}
}
