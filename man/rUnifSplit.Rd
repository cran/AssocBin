% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitters.R
\name{rUnifSplit}
\alias{rUnifSplit}
\title{Random uniform splitting}
\usage{
rUnifSplit(bin, minExp = 0, squarify = FALSE)
}
\arguments{
\item{bin}{a bin to be split with elements `x`, `y`, `depth`,
`bnds` (list with elements `x` and `y`), `expn`, `n`}

\item{minExp}{numeric giving the minimum expected count allowed
in a bin}

\item{squarify}{logical value, should we force splitting on
the longer side?}
}
\value{
A list of two bins resulting from the split of `bin`
at a random location on a random margin
}
\description{
Split bins randomly and uniformly
}
\details{
This function samples a coordinate uniformly along a
random margin and splits a bin at that coordinate. In contrast to
maxScoreSplit with randScores, this can introduce splits at
locations other than the points.
}
\examples{
bin <- makeBin(x = 1:10, y = sample(1:10))
rUnifSplit(bin, minExp = 2)
}
\author{
Chris Salahub
}
