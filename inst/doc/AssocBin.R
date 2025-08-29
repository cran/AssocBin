## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(AssocBin)

## -----------------------------------------------------------------------------
data(heart)

## -----------------------------------------------------------------------------
str(heart)

## -----------------------------------------------------------------------------
heartClean <- heart
heartClean$thal <- NULL
heartClean$ca <- NULL
heartClean$slope <- NULL
heartClean <- na.omit(heartClean)
str(heartClean)

## ----fig.width=6, fig.height=6, fig.align = 'center'--------------------------
depDisplay(heartClean$sex, heartClean$num)

## ----fig.width=6, fig.height=6, fig.align = 'center'--------------------------
SexVsNum <- depDisplay(heartClean$sex, heartClean$num, xlab = "Sex", 
                       ylab = "Number of arteries >50% obstructed", 
                       pch = 20)

## -----------------------------------------------------------------------------
rbind(cbind(table(num = heartClean$num, sex = heartClean$sex), total = table(heartClean$num)),
            total = c(table(heartClean$sex), nrow(heartClean)))

## -----------------------------------------------------------------------------
str(SexVsNum, 1)

## -----------------------------------------------------------------------------
str(SexVsNum[[1]])

## -----------------------------------------------------------------------------
binChi(SexVsNum)

## ----fig.width=6, fig.height=6, fig.align = 'center'--------------------------
set.seed(1235) # more on this later
# the depDisplay function also has a method for data.frames
AgeVsNum <- depDisplay(x = heartClean, pair="age:num", xlab = "Age", 
                       ylab = "Number of arteries >50% obstructed", 
                       pch = 20, col = adjustcolor('gray50', alpha.f=0.5))

## ----fig.width=6, fig.height=6, fig.align = 'center'--------------------------
set.seed(812)
AgeVsChol <- depDisplay(heartClean$thalach, heartClean$oldpeak, 
                        xlab = "Maximum heart rate during exercise",
                        ylab = "ST wave depression during exercise",
                        pch = 20, col = adjustcolor('gray50', alpha.f=0.5))

## ----fig.width=6, fig.height=6, fig.align = 'center'--------------------------
set.seed(812)
AgeVsNum <- depDisplay(heartClean$thalach, heartClean$oldpeak, 
                       xlab = "Maximum heart rate during exercise",
                       ylab = "ST wave depression during exercise",
                       pch = 20, col = adjustcolor('gray50', alpha.f=0.5),
                       border = "black")

## -----------------------------------------------------------------------------
heartAssociations <- DepSearch(heartClean)

## -----------------------------------------------------------------------------
summary(heartAssociations)

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
plot(heartAssociations)

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
plot(heartAssociations, which = 62:66)

## -----------------------------------------------------------------------------
stopCrits <- makeCriteria(depth >= 10, # maximum depth of 10
                          expn <= 10, # smallest possible bin size of 5
                          n < 1 # don't split empty bins
                          )
stopCrits

## -----------------------------------------------------------------------------
greedyCrits <- makeCriteria(abs(expn - n)/sqrt(expn) > 4,
                            expn <= 10,
                            n < 1)
greedyCrits

## -----------------------------------------------------------------------------
conConChi <- function(bn) maxScoreSplit(bin = bn, scorer = chiScores)
# the univariate splitter requires an additional argument specifying which
# margin should be split
catConChi <- function(bn, on) uniMaxScoreSplit(bin = bn, scorer = chiScores,
                                               on = on)

## -----------------------------------------------------------------------------
heartAssociations_greedy <- DepSearch(heartClean,
                                      stopCriteria=greedyCrits,
                                      catCon=catConChi,
                                      conCon=conConChi)

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
plot(heartAssociations_greedy)

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
# a final way to use depDisplay is on a depSearch object
depDisplay(heartAssociations, pair="thalach:oldpeak",
           xlab = "Maximum heart rate during exercise",
           ylab = "ST wave depression during exercise",
           pch = "+", col = adjustcolor('purple', alpha.f=0.5),
           border = "black")

## -----------------------------------------------------------------------------
thalachOldpeak <- heartAssociations$binnings[["thalach:oldpeak"]]

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
# note that plotBinning does not have access to the marginal information to plot
# quantiles and so the marginal labels give the ranks
plotBinning(thalachOldpeak, pch = 20, 
            xlab = "Maximum heart rate during exercise",
            ylab = "ST wave depression during exercise",
            showXax = TRUE, showYax = TRUE,
            fill=depthFill(thalachOldpeak))

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
x <- rnorm(1000)
y <- 2*x + rnorm(1000, sd = 0.3)
rankx <- rank(x, ties.method = "random")
ranky <- rank(y, ties.method = "random")

# set up splitting criteria: depth stop limits run time (not necessary here)
criteria <- makeCriteria(expn <= 10, n == 0, depth >= 10)
# define the stop function using these criteria
stopFn <- function(bns) stopper(bns, criteria)
# use binner to run the algorithm
xyBins <- binner(x = rankx, y = ranky, stopper = stopFn, splitter = rIntSplit)

# plot with depthfill
set.seed(2119)
plotBinning(xyBins, fill=depthFill(xyBins), pch = 20)

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
plotBinning(thalachOldpeak, pch = 20, 
            xlab = "Maximum heart rate during exercise",
            ylab = "ST wave depression during exercise",
            showXax = TRUE, showYax = TRUE,
            fill=residualFill(thalachOldpeak, nbr = 10))

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
plotBinning(thalachOldpeak, pch = 20, 
            xlab = "Maximum heart rate during exercise",
            ylab = "ST wave depression during exercise",
            showXax = TRUE, showYax = TRUE,
            fill=residualFill(thalachOldpeak, nbr = 50,
                              resFun=binMI,
                              colrng = c("orange", "pink", "blue")))

