% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stopper.R
\name{stopper}
\alias{stopper}
\title{Check bins against stop criteria}
\usage{
stopper(binList, criteria)
}
\arguments{
\item{binList}{a list of bins, each a list which can be cast as an
environment for evaluation}

\item{criteria}{string of logical expressions separated by pipes
to be evaluated within each bin of `binList`}
}
\value{
A logical vector of the same length as `binList`.
}
\description{
Evaluate the stop `criteria` for each bin in
`binList`
}
\details{
This function makes use of R's lexical scoping to
evaluate `criteria` (a string), within each bin of `binList`.
}
\examples{
crits <- makeCriteria(depth >= 5, n < 1)
binList1 <- list(makeBin(x = c(1,2), y = c(3,1), depth = 1, n = 2),
                makeBin(x = c(3,4), y = c(2,4), depth = 1, n = 2))
binList2 <- list(makeBin(x = c(1,2), y = c(3,1), depth = 6, n = 2),
                makeBin(x = c(), y = c(), depth = 1, n = 0))
stopper(binList1, crits)
stopper(binList2, crits)
}
\author{
Chris Salahub
}
