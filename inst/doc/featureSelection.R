## ----makeHappy,echo=FALSE------------------------------------------------
if (!require(Umpire)) {
  knitr::opts_chunk$set(eval = FALSE)
}

## ----lib-----------------------------------------------------------------
library(GenAlgo)
library(Umpire)
library(oompaBase)

## ----survmodel-----------------------------------------------------------
set.seed(391629)
sm <- SurvivalModel(baseHazard=1/5, accrual=5, followUp=1)

## ----cancModel-----------------------------------------------------------
nBlocks <- 20    # number of possible hits
cm <- CancerModel(name="cansim",
                  nPossible=nBlocks,
                  nPattern=6,
                  OUT = function(n) rnorm(n, 0, 1), 
                  SURV= function(n) rnorm(n, 0, 1),
                  survivalModel=sm)
### Include 100 blocks/pathways that are not hit by cancer
nTotalBlocks <- nBlocks + 100

## ----hyper---------------------------------------------------------------
### block size
blockSize <- round(rnorm(nTotalBlocks, 100, 30))
### log normal mean hypers
mu0    <- 6
sigma0 <- 1.5
### log normal sigma hypers
rate   <- 28.11
shape  <- 44.25
### block corr
p <- 0.6
w <- 5

## ----engine--------------------------------------------------------------
rho <- rbeta(nTotalBlocks, p*w, (1-p)*w)
base <- lapply(1:nTotalBlocks,
               function(i) {
                 bs <- blockSize[i]
                 co <- matrix(rho[i], nrow=bs, ncol=bs)
                 diag(co) <- 1
                 mu <- rnorm(bs, mu0, sigma0)
                 sigma <- matrix(1/rgamma(bs, rate=rate, shape=shape), nrow=1)
                 covo <- co *(t(sigma) %*% sigma)
                 MVN(mu, covo)
               })
eng <- Engine(base)

## ----alter---------------------------------------------------------------
altered <- alterMean(eng, normalOffset, delta=0, sigma=1)
object <- CancerEngine(cm, eng, altered)
summary(object)

## ----clean1--------------------------------------------------------------
rm(altered, base, blockSize, cm, eng, mu0, nBlocks, nTotalBlocks,
   p, rate, rho, shape, sigma0, sm, w)

## ----traind--------------------------------------------------------------
train <- rand(object, 198)
tdata <- train$data
pid <- paste("PID", sample(1001:9999, 198+93), sep='')
rownames(train$clinical) <- colnames(tdata) <- pid[1:198]

## ----noise---------------------------------------------------------------
noise <- NoiseModel(3, 1, 1e-16)
train$data <- log2(blur(noise, 2^(tdata)))
sum(is.na(train$data))
rm(tdata)
summary(train$clinical)
summary(train$data[, 1:3])

## ----validd--------------------------------------------------------------
valid <- rand(object, 93)
vdata <- valid$data
vdata <- log2(blur(noise, 2^(vdata))) # add noise
sum(is.na(vdata))
vdata[is.na(vdata)] <- 0.26347
valid$data <- vdata
colnames(valid$data) <- rownames(valid$clinical) <- pid[199:291]
rm(vdata, noise, object, pid)
summary(valid$clinical)
summary(valid$data[, 1:3])

## ----measureFitness------------------------------------------------------
measureFitness <- function(arow, context) {
  predictors <- t(context$data[arow, ]) # space defined by features
  groups <- context$clinical$Outcome    # good or bad outcome
  maha(predictors, groups, method='var')
}

## ----mutator-------------------------------------------------------------
mutator <- function(allele, context) {
   sample(1:nrow(context$data),1)
}

## ----initialize----------------------------------------------------------
set.seed(821831)
n.individuals <- 200
n.features <- 10
y <- matrix(0, n.individuals, n.features)
for (i in 1:n.individuals) {
  y[i,] <- sample(1:nrow(train$data), n.features)
}

## ----round1--------------------------------------------------------------
my.ga <- GenAlg(y, measureFitness, mutator, context=train) # initialize
summary(my.ga)

## ----save0---------------------------------------------------------------
recurse <- my.ga
pop0 <- sort(table(as.vector(my.ga@data)))

## ----recurse-------------------------------------------------------------
NGEN <- 20
diversity <- meanfit <- fitter <- rep(NA, NGEN)
for (i in 1:NGEN) {
  recurse <- newGeneration(recurse)
  fitter[i] <- recurse@best.fit
  meanfit[i] <- mean(recurse@fitness)
  diversity[i] <- popDiversity(recurse)
}

## ----fig.cap="Fitness by generation."------------------------------------
plot(fitter, type='l', ylim=c(0, 1.5), xlab="Generation", ylab="Fitness")
abline(h=max(fitter), col='gray', lty=2)
lines(fitter)
lines(meanfit, col='gray')
points(meanfit, pch=16, col=jetColors(NGEN))
legend("bottomleft", c("Maximum", "Mean"), col=c("black", "blue"), lwd=2)

## ----fig.cap="Diversity."------------------------------------------------
plot(diversity, col='gray', type='l', ylim=c(0,10), xlab='', ylab='', yaxt='n')
points(diversity, pch=16, col=jetColors(NGEN))

## ----div-----------------------------------------------------------------
sort(table(as.vector(recurse@data)))

