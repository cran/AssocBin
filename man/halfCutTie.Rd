% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitters.R
\name{halfCutTie}
\alias{halfCutTie}
\title{Halve continuously to break ties}
\usage{
halfCutTie(bin, xscore, yscore)
}
\arguments{
\item{bin}{a bin to be split with elements `x`, `y`, `depth`,
`bnds` (list with elements `x` and `y`), `expn`, `n`}

\item{xscore}{numeric value giving the score for all splits along
x}

\item{yscore}{numeric value giving the score for all splits along
y}
}
\value{
A list of two bins resulting from the split of `bin` in
half along the margin corresponding to the larger score.
}
\description{
This function halves a bin based on the midpoint of
the bounds along whichever margin produces the larger score.
}
\details{
The goal of this function is to break ties within bin
splitting in a way which prevents very small or lopsided bins from
forming, a common problem with the `halfSplit` function
}
\examples{
bin <- list(x = 1:10, y = sample(1:10),
            bnds = list(x = c(0, 10), y = c(0, 10)),
            expn = 10, n = 10, depth = 0)
halfCutTie(bin, 1, 2) # splits on y
halfCutTie(bin, 2, 1) # splits on x
halfCutTie(bin, 1, 1) # ties are random
}
\author{
Chris Salahub
}
