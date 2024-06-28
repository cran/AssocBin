## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(72604204)

## -----------------------------------------------------------------------------
library(AssocBin)

## -----------------------------------------------------------------------------
data(heart)
summary(heart)

## -----------------------------------------------------------------------------
heartClean <- heart
heartClean$thal <- NULL
heartClean$ca <- NULL
heartClean$slope <- NULL
heartClean <- na.omit(heartClean)
str(heartClean)

## -----------------------------------------------------------------------------
stopCrits <- makeCriteria(depth >= 6, n < 1, expn <= 10)
assocs <- inDep(heartClean, stopCriteria = stopCrits)

## -----------------------------------------------------------------------------
summary(assocs)

## ----fig.width = 5, fig.height = 8--------------------------------------------
plot(assocs) # by default this displays the 5 strongest relationships

## ----fig.width = 5, fig.height = 8--------------------------------------------
## by specifying which indices should be displayed, others can be plotted
## the binnings are returned in increasing order of p-value, so the indices
## chosen give the rank of the strength a particular pair's relationship
plot(assocs, which = 6:10) # like the next 5 strongest
plot(assocs, which = 62:66) # or the 5 weakest relationships

## -----------------------------------------------------------------------------
maxCatCon <- function(bn) uniMaxScoreSplit(bn, chiScores)
maxConCon <- function(bn) maxScoreSplit(bn, chiScores)
maxPairs <- inDep(data = heartClean, stopCriteria = stopCrits, 
                  catCon = maxCatCon, conCon = maxCatCon)
summary(maxPairs)

## ----fig.width = 5, fig.height = 8--------------------------------------------
plot(maxPairs)
plot(maxPairs, which = 6:10)

## ----fig.width = 5, fig.height = 5--------------------------------------------
randOrd <- match(assocs$pairs, assocs$pairs)
maxOrd <- match(assocs$pairs, maxPairs$pairs)
plot(randOrd, maxOrd, xlab = "Random rank", ylab = "Max chi rank", 
     main = "Rankings of pair significance between methods")

## -----------------------------------------------------------------------------
heartCleve <- heartClean[heartClean$study == "cleveland", ]
heartCleve$study <- NULL
cleveAssoc <- inDep(heartCleve, stopCriteria = stopCrits, catCon = maxCatCon,
                    conCon = maxConCon)
summary(cleveAssoc)

## ----fig.width = 5, fig.height = 8--------------------------------------------
plot(cleveAssoc)
plot(cleveAssoc, which = 6:10)
plot(cleveAssoc, which = 11:15)

