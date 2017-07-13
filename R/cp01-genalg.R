# Copyright (C) Kevin R. Coombes, 2007-2012

###
### GENALG.R
###


########################################################################
# Here are the functions to implement the genetic algorithm generically

# The first basic "mating" function, which implements crossover as 
# a single cut in the list of features. The default is not to
# shuffle the feature order, and simply make a single cut.
crossover <- function(a, b, shuffle=FALSE) {
  cutpoint <- sample(1:(length(a)-1), 1)
  if (shuffle) {
    which <- rep(FALSE, length(a))
    which[sample(1:length(a), cutpoint)] <- TRUE
    alpha <- c(a[which], b[!which])
    beta  <- c(b[which], a[!which])
  } else {
    alpha <- c(a[1:cutpoint], b[(cutpoint+1):length(a)])
    beta <- c(b[1:cutpoint], a[(cutpoint+1):length(a)])
  }
  list(a=alpha, b=beta)
}

# The second basic "mating" function. A discrete cumulative probability
# vector is supplied as input, and we convert from a uniform distribution
# to select a pair of individuals from the population.
selectPair <- function(probvec) {
  i <- 1 + sum(probvec < runif(1))
  if (i > length(probvec)) {
    i <- 1
  }
  j <- 1 + sum(probvec < runif(1))
  if (i == j) {
    j <- j + 1
  }
  if (j > length(probvec)) {
    j <- 1
  }
  list(i=i, j=j)
}

# The "culling" routine. The use of "solve" in the computation of the
# mahalanobis distance means that the covariance matrix for the selected
# features has to be nonsingular. Although we should probably test that
# directly, we have so far only run into problems when we manage to select
# the same feature twice in the same set. We work around this by "mutating"
# any individuals that suffer from this 'fatal' gene combination.
removeDuplicates <- function(arow, mf, context) {
  needs <- length(arow) - length(unique(arow))
  while (needs > 0) {
    arow <- c(unique(arow), unlist(lapply(rep(0,needs), mf, context)))
    needs <- length(arow) - length(unique(arow))
  }
  arow
}

# This routine serves two roles. First, it is a constructor function for
# GenAlg objects. You have to supply a data matrix, where we think of each
# row as an individual in the population, and the columns mean whatever you
# want them to mean. In our application, the entries are the indexes into the
# list of features. You also have to supply
#	fitfun = a fuction that will be applied to each row in order to compute
#		the fitness of that individual
#	mutfun = a function that will be applied to each individual produced
#		by mating to cause random mutations in the inherited alleles
#       context = a list containing auxiilary data to be passed throught
#               to the the mutation and fitness functions
#	pm = the per-allele probability of a mutation
#       pc = the per-mating probability of a crossover event
#	gen = an integer represnting which generation this is in the running of
#		the algorithm.
#
# When the GenAlg is initialized, we compute the fitness of all the
# individuals in the population. Thus, the GenAlg function also serves as a
# critical routine when you iterate the generate algorithm to process several
# generations.
setClass("GenAlg",
         representation(data="matrix",
                        fitfun="function",
                        mutfun="function",
                        p.mutation="numeric",
                        p.crossover="numeric",
                        generation="numeric",
                        fitness="numeric",
                        best.fit="numeric",
                        best.individual="matrix",
                        context="list"))

GenAlg <- function(data, fitfun, mutfun, context, pm=0.001, pc=0.50, gen=1) {
  fitness <- apply(data, 1, fitfun, context=context)
  best.fit <- max(fitness)
  cheat <- (1:length(fitness))[fitness==best.fit][1]
  b1 <- as.matrix(data[cheat,])
  new("GenAlg",
      data=data, fitfun=fitfun, mutfun=mutfun, p.mutation=pm,
      p.crossover=pc, generation=gen, fitness=fitness,
      best.fit=best.fit, best.individual=b1, context=context)
}

# This is the iterative routine that controls the overall genetic algorithm.
# Having first initialized the genetic algorithm to produce a GenAlg object,
# you feed it to this routine and get back anotheer GenAlg object that
# represents the next generation.
newGeneration <- function(ga) {
  x <- ga@data
  fit <- ga@fitness
  probs <- cumsum(fit)/sum(fit) # prob of selection given by relative fitness
  temp <- matrix(unlist(lapply(1:(nrow(x)/2), function(i, p, data) {
    ij <- selectPair(p)
    if (runif(1) <=  ga@p.crossover) {
      ab <- crossover(data[ij$i,], data[ij$j,])
    } else {
      ab <- list(a=data[ij$i,], b=data[ij$j,])
    }
  }, probs, x)), nrow=nrow(x), ncol=ncol(x), byrow=T)
  targets <- sample(1:(nrow(x)*ncol(x)), ceiling(ga@p.mutation*nrow(x)*ncol(x)))
  start.size <- dim(temp)[2]
  temp[targets] <- ga@mutfun(temp[targets], ga@context)
  temp <- as.matrix(t(apply(temp, 1, removeDuplicates, ga@mutfun, ga@context)))
  GenAlg(temp, ga@fitfun, ga@mutfun, ga@context,
         ga@p.mutation, ga@p.crossover, ga@generation+1)
}

setMethod("summary", signature(object="GenAlg"),
          function(object, ...) {
  cat(paste("An object representing generation", object@generation,
            "in a genetic algorithm.\n"))
  cat(paste("Population size:", nrow(object@data), "\n"))
  cat(paste("Mutation probability:", object@p.mutation, "\n"))
  cat(paste("Crossover probability:", object@p.crossover, "\n"))
  cat("Fitness distribution:\n")
  print(summary(object@fitness))
})

setMethod('as.data.frame', signature(x='GenAlg'),
          function(x, row.names=NULL, optional=FALSE, ...) {
  val <- data.frame(x@fitness, x@data)
  size <- dim(x@data)[2]
  colnames(val) <- c('Fitness', paste('Feature', 1:size, sep=''))
  val
})

setMethod("as.matrix", signature(x="GenAlg"),
          function(x, ...) {
	as.matrix(as.data.frame(x))
})

popDiv <- function(x) {
  N <- nrow(x)
  ndiff <- 0
  for (i in 1:(N-1)) {
    more <- sapply((i+1):N, function(j) {
      length(unique(c(x[i,], x[j,]))) - ncol(x)
    })
    ndiff <- ndiff+sum(more)
  }
  ndiff/(N*(N-1)/2)  
}

popDiversity <- function(ga) {
  popDiv(ga@data)
}

########################################################################
# this is a residue of a previous demonstration version of this code.
# it only used alleles taking values 0 or 1 in a binary representation
# of an approximate real number solution to a problem. The fitness function
# that acompanied it has long since vanished.
simpleMutate <- function(allele, context) { 1-allele }

########################################################################
# Now we add additional functions for a specific genetic algorithm
# that we are often interested in. Note that the fitness and mutation
# functions assume that they have access to the following variables inside
# the 'context' list:
#	dataset = the matrix of feature values for individual samples. The
#		entries in the genetic algorithm data matrix (ie, the alleles)
#		are indices into the rows of the dataset matrix.
#	gps = a logical vector classifying the columns of the dataset into
#		two distinct groups, such as cancer or normal.
#
# In order to create an instance of the class of genetic algorithms that
# uses the fitness and mutation functions defined in this section, you must
# create the "dataset" and "gps" objects, and must then also initialize
# a starting population for the gentic algorithm.


# This function uses Mahalanobis distance to compute the fitness of
# an individual.
selectionFitness <- function(arow, context) {
  maha(t(context$dataset[arow,]), context$gps, method='var')
}

# This is the mutation routine. It assumes that every feature is equally
# likely to mutate into every other feature.
selectionMutate <- function(allele, context) {
	sample(1:nrow(context$dataset),1)
}

