\name{GenAlg-tools}
\alias{GenAlg-tools}
\alias{simpleMutate}
\alias{selectionFitness}
\alias{selectionMutate}
\title{Utility functions for selection and mutation in genetic algorithms}
\description{
  These functions implement specific forms of mutation and fitness
  that can be used in genetic algorithms for feature selection.
}
\usage{
simpleMutate(allele, context)
selectionMutate(allele, context)
selectionFitness(arow, context)
}
\arguments{
  \item{allele}{
    In the \code{simpleMutate} function, \code{allele} is a binary
    vector filled with 0's and 1's.  In the \code{selectionMutate}
    function, \code{allele} is an integer (which is silently ignored;
    see Details). 
  }
  \item{arow}{
    A vector of integer indices identifying the rows (features) to be
    selected from the \code{context$dataset} matrix.
  }
  \item{context}{
    A list or data frame containing auxiliary information that is needed
    to resolve references from the mutation or fitness code.  In both 
    \code{selectionMutate} and \code{selectionFitness}, \code{context}
    must contain a \code{dataset} component that is either a matrix or a
    data frame.  In \code{selectionFitness}, the \code{context} must
    also include a grouping factor (with two levels) called \code{gps}.
  }
}
\details{
  These functions represent 'callbacks'. They can be used in the
  function \code{\link{GenAlg}}, which creates objects. They will then
  be called repeatedly (for each individual in the population) each time
  the genetic algorithm is updated to the next generation.

  The \code{simpleMutate} function assumes that chromosomes are binary
  vectors, so alleles simply take on the value 0 or 1. A mutation of an
  allele, therefore, flips its state between those two possibilities.

  The \code{selectionMutate} and \code{selectionFitness} functions, by
  contrast, are specialized to perform feature selection assuming a
  fixed number K of features, with a goal of learning how to
  distinguish between two different groups of samples. We assume that
  the underlying data consists of a data frame (or matrix), with the
  rows representing features (such as genes) and the columns
  representing samples. In addition, there must be a grouping vector
  (or factor) that assigns all of the sample columns to one of two
  possible groups. These data are collected into a list,
  \code{context}, containing a \code{dataset} matrix and a \code{gps}
  factor. An individual member of the population of potential
  solutions is encoded as a length K vector of indices into the rows
  of the \code{dataset}. An individual \code{allele}, therefore, is a
  single index identifying a row of the \code{dataset}. When mutating
  it, we assume that it can be changed into any other possible allele;
  i.e., any other row number. To compute the fitness, we use the
  Mahalanobis distance between the centers of the two groups defined by
  the \code{gps} factor.
}
\value{
  Both \code{selectionMutate} and \code{simpleMutate} return an integer
  value; in the simpler case, the value is guaranteed to be a 0 or 1.
  The \code{selectionFitness} function returns a real number.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{GenAlg}},
  \code{\link{GenAlg-class}},
  \code{\link{maha}}.
}
\examples{
# generate some fake data
nFeatures <- 1000
nSamples <- 50
fakeData <- matrix(rnorm(nFeatures*nSamples), nrow=nFeatures, ncol=nSamples)
fakeGroups <- sample(c(0,1), nSamples, replace=TRUE)
myContext <- list(dataset=fakeData, gps=fakeGroups)

# initialize population
n.individuals <- 200
n.features <- 9
y <- matrix(0, n.individuals, n.features)
for (i in 1:n.individuals) {
  y[i,] <- sample(1:nrow(fakeData), n.features)
}

# set up the genetic algorithm
my.ga <- GenAlg(y, selectionFitness, selectionMutate, myContext, 0.001, 0.75)

# advance one generation
my.ga <- newGeneration(my.ga)

}
\keyword{optimize}

