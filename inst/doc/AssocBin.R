## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(AssocBin)

## -----------------------------------------------------------------------------
set.seed(9023831)
n <- 100
randx <- rnorm(n)
randy <- rnorm(n)
plot(randx, randy)

## -----------------------------------------------------------------------------
rankx <- rank(randx, ties.method = "random")
ranky <- rank(randy, ties.method = "random")
plot(rankx, ranky)

## -----------------------------------------------------------------------------
criteria <- makeCriteria(expn <= 10, n == 0, depth >= d)
str(criteria)

## -----------------------------------------------------------------------------
stopper
stopFn <- function(bns) stopper(bns, criteria)

## -----------------------------------------------------------------------------
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
miSplit <- function(bn) maxScoreSplit(bn, miScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)

## -----------------------------------------------------------------------------
d <- 5
randBin <- binner(x = rankx, y = ranky, stopper = stopFn, splitter = chiSplit)

## -----------------------------------------------------------------------------
plotBinning(randBin, pch = 19, cex = 1,
            col = adjustcolor("gray50", 0.5))

## -----------------------------------------------------------------------------
## first fill by depth
plotBinning(randBin, pch = 19, cex = 1,
            fill = depthFill(randBin),
            col = adjustcolor("gray50", 0.5))

## -----------------------------------------------------------------------------
## next fill by residual
plotBinning(randBin, pch = 19, cex = 1,
            fill = residualFill(randBin, colrng = c("steelblue", "white", "firebrick")),
            col = adjustcolor("gray50", 0.5))

## -----------------------------------------------------------------------------
## next fill by residual
plotBinning(randBin, pch = 19, cex = 1,
            fill = residualFill(randBin, colrng = c("steelblue", "white", "firebrick"),
                                nbr = 50),
            col = adjustcolor("gray50", 0.5))

## -----------------------------------------------------------------------------
depx <- rnorm(10*n)
depy <- depx + rnorm(10*n, sd = 0.4)
plot(depx, depy)

## -----------------------------------------------------------------------------
d <- 10 # change maximum depth due to larger sample size
depx.rank <- rank(depx)
depy.rank <- rank(depy)
depBins <- binner(depx.rank, depy.rank, stopper = stopFn, splitter = chiSplit)
plotBinning(depBins, pch = ".", cex = 1)

## -----------------------------------------------------------------------------
plotBinning(depBins, pch = ".", cex = 1,
            fill = residualFill(depBins))

## -----------------------------------------------------------------------------
set.seed(591241)
depBins.rand <- binner(depx.rank, depy.rank, stopper = stopFn, splitter = rndSplit)
plotBinning(depBins.rand, pch = ".", cex = 1, 
            fill = residualFill(depBins.rand))

## -----------------------------------------------------------------------------
binChi(randBin)$stat
binChi(depBins)$stat
binChi(depBins.rand)$stat

