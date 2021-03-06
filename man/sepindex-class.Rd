% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sepindex.R
\docType{class}
\name{sepindex-class}
\alias{sepindex-class}
\alias{sepindex}
\alias{boxplot,sepindex-method}
\alias{sepindex-method}
\alias{show,sepindex-method}
\alias{[,sepindex-method}
\alias{summary,sepindex-method}
\title{Class "sepindex"}
\usage{
sepindex(vario_st, nt, ns, globalSill)

\S4method{boxplot}{sepindex}(x, ...)

\S4method{show}{sepindex}(object)

\S4method{[}{sepindex}(x, i, j)

\S4method{summary}{sepindex}(object)
}
\arguments{
\item{vario_st}{spatio-temporal sample variogram, output from
\code{\link[gstat]{variogramST}}}

\item{nt}{integer, the number of temporal lags in \code{vario_st}}

\item{ns}{integer, the number of spatial lags in \code{vario_st}}

\item{globalSill}{numeric, the value of the sample variance}

\item{x}{object of class \code{sepindex} for methods \code{boxplot} and \code{extract}}

\item{...}{any arguments that will be passed to the panel plotting functions}

\item{object}{object of class \code{sepindex} for methods \code{show} and \code{summary}}

\item{i}{index specifing elements to extract. Each row includes data for specific
spatio-temporal lags}

\item{j}{index specifing elements to extract. Set \code{1} for spatial lags (hs),
\code{2} for temporal lags (ht) and \code{3} for the non-separability index
(SepIndex)}
}
\description{
A class for the non-separability index (r) for different spatial
and temporal lags:
\deqn{r(h, u, \Theta)= \rho(h, u;\Theta)/ [\rho(h,0;\Theta)\rho(0,u;\Theta)]}
with \eqn{\rho(h, u;\Theta)>0}; \eqn{\rho(h,0;\Theta)>0} and \eqn{\rho(0,u;\Theta)>0}.
On the basis of this index, the type of non-separability of the covariance
function can be analyzed.
}
\section{Slots}{

\describe{
\item{\code{sep.index.ratio}}{the empirical non-separability index ratio and the
corresponding spatio-temporal lags}

\item{\code{cov.st}}{the spatio-temporal sample covariance function and the
corresponding spatio-temporal lags}

\item{\code{cov.tm}}{the purely temporal sample covariance function and the
corresponding temporal lags}

\item{\code{cov.sp}}{the purely spatial sample covariance function and the
corresponding spatial lags}
}}

\examples{
# --start define the STFDF rr_13-- #
library(sp)
library(spacetime)
library(gstat)
data(air)
ls()
if (!exists("rural")) rural = STFDF(stations, dates, data.frame(PM10 =
as.vector(air)))
rr = rural[,"2005::2010"]
unsel = which(apply(as(rr, "xts"), 2, function(x) all(is.na(x))))
r5to10 = rr[-unsel,]
rr_13 <- r5to10[c("DEHE046","DESN049","DETH026","DENW063","DETH061","DEBY047",
"DENW065","DEUB029","DENW068","DENI019","DEHE051","DERP016","DENI051"),
"2005::2006"]
# --end define the STFDF rr_13-- #

#compute the Global Sill
C00_13<-var(rr_13[,,"PM10"]@data[[1]], na.rm = TRUE)

#estimate the spatio-temporal variogram
data(vv_13)
nonsep.index<-sepindex(vario_st = vv_13, nt = 16, ns = 4, globalSill = C00_13)

##methods for sepindex

#1. show
nonsep.index

#2. summary
summary(nonsep.index)

#3. boxplot
boxplot(nonsep.index, ylab="Non-separability ratio")

#4. [ extract
nonsep.index[1:8, ] #selection of the first 8 rows
nonsep.index[1:8, 1:2] #selection of the first 2 columns
}
\references{
Cappello, C., De Iaco, S., Posa, D., 2020, {covatest}: An {R} Package for
Selecting a Class of Space-Time Covariance Functions.
Journal of Statistical Software, \bold{94(1)} 1--42.

De Iaco, S., Posa, D., 2013, Positive and negative non-separability
for space-time covariance models. Journal of Statistical Planning and
Inference, \bold{143} 378--391.

{Gräler}, B., Pebesma, E.J., Heuvelink G., 2016, Spatio-Temporal
Interpolation Using gstat. The {R} Journal, \bold{8(1)} 204--218.

Pebesma, E.J., 2004, Multivariable geostatistics in S: the gstat
package. Computers & Geosciences, \bold{30} 683--691.

Rodriguez, A., Diggle, P.J., 2010, A class of convolution-based
models for spatio-temporal processes with non-separable covariance structure.
Scandinavian Journal of Statistics, \bold{37(4)} 553--567.
}
\seealso{
\code{\link[gstat]{variogramST}}
}
