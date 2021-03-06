% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMarginalZ.R
\name{getMarginalZ}
\alias{getMarginalZ}
\title{Construct a (smooth) marginal z density approximation from a model
information list}
\usage{
getMarginalZ(
  info,
  method = c("linear", "spline", "logspline", "normalspline", "normal"),
  verbose = FALSE,
  plot = FALSE
)
}
\arguments{
\item{info}{the model information list}

\item{method}{method for approximating the marginal density:
\describe{
\item{linear}{Linearly interpolate the points.}
\item{spline}{The saved points of the unnormalized density approximation are joined
by a \dQuote{monotonic} spline. The density is smoothed out to zero at the
tails. Since the spline might be slightly negative for extreme values, the positive part
is returned.}
\item{logspline}{The saved points of the log unnormalized density approximation are joined 
by a \dQuote{monotonic} spline, which is then exponentiated.}
\item{normalspline}{A \dQuote{monotonic} spline is fitted to the differences of the saved
log density values and the log normal approximation. The resulting spline function is
exponentiated and then multiplied with the normal density.}
\item{normal}{Just take the normal approximation.}
This may also be a vector with more than one method names, to select the modify the preference
sequence: If the first method does not work, the second is tried and so on. The normal
approximation always \dQuote{works} (but may give bad results).
}}

\item{verbose}{Echo the chosen method? (not default)}

\item{plot}{produce plots of the different approximation steps? (not default)}
}
\value{
a list with the log of the normalized density approximation (\dQuote{logDens}) and 
the random number generator (\dQuote{gen}).
}
\description{
Construct a (smooth) marginal z density approximation from a model
information list
}
\keyword{internal}
