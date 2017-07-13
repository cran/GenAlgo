\name{GenAlg-class}
\alias{GenAlg-class}
\alias{as.data.frame,GenAlg-method}
\alias{as.matrix,GenAlg-method}
\alias{summary,GenAlg-method}
\docType{class}
\title{Class "GenAlg"}
\description{
  Objects of the \code{GenAlg} class represent one step (population) in the
  evolution of a genetic algorithm. This algorithm has been customized to
  perform feature selection for the class prediction problem.
}
\usage{
\S4method{as.data.frame}{GenAlg}(x, row.names=NULL, optional=FALSE, \dots)
\S4method{as.matrix}{GenAlg}(x, \dots)
\S4method{summary}{GenAlg}(object, \dots)
}
\arguments{
  \item{object}{object of class \code{GenAlg}}
  \item{x}{object of class \code{GenAlg}}
  \item{row.names}{character vector giving the row names for the data frame,
    or \code{NULL}}
  \item{optional}{logical scalar. If \code{TRUE}, setting row names and
    converting column names to syntactic names is optional.}
  \item{\dots}{extra arguments for generic routines}
} 
\section{Objects from the Class}{
  Objects should be created by calls to the \code{\link{GenAlg}} generator;
  they will also be created automatically as a result of applying the function
  \code{\link{newGeneration}} to an existing \code{GenAlg} object.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{The initial population of potential solutions,
      in the form of a data matrix with one individual per row.}
    \item{\code{fitfun}:}{A function to compute the fitness of an
      individual solution. Must take two input arguments: a vector of
      indices into the rows of the population matrix, and a \code{context}
      list within which any other items required by the function can be
      resolved. Must return a real number; higher values indicate better
      fitness, with the maximum fitness occurring at the optimal solution
      to the underlying numerical problem.} 
    \item{\code{mutfun}:}{A function to mutate individual alleles in the
      population. Must take two arguments: the starting allele and a
      \code{context} list as in the fitness function.}
    \item{\code{p.mutation}:}{numeric scalar between \code{0} and \code{1},
      representing the probability that an individual allele will be mutated.}
    \item{\code{p.crossover}:}{numeric scalar between \code{0} and \code{1},
      representing the probability that crossover will occur during
      reproduction.}
    \item{\code{generation}:}{integer scalar identifying the current
      generation.}
    \item{\code{fitness}:}{numeric vector containing the fitness of
      all individuals in the population.}
    \item{\code{best.fit}:}{A numeric value; the maximum fitness.}
     \item{\code{best.individual}:}{A matrix (often with one row)
       containing the individual(s) achieving the maximum fitness. }
    \item{\code{context}:}{A list of additional data required to perform
      mutation or to compute fitness. This list is passed along as the
      second argument when \code{fitfun} and \code{mutfun} are called.}
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "GenAlg")}: Converts the
      \code{GenAlg} object into a data frame.  The first column contains
      the fitness ; remaining columns contain three selected features,
      given as integer indices into the rows of the original data matrix. }
    \item{as.matrix}{\code{signature(x = "GenAlg")}: Converts the GenAlg
      object into a matrix, following the conventions of \code{as.data.frame}.}
    \item{summary}{\code{signature(object = "GenAlg")}: Print a summary
      of the GenAlg object.}
  }
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\references{
  David Goldberg.\cr
  "Genetic Algorithms in Search, Optimization and Machine Learning."\cr
  Addison-Wesley, 1989.
}
\seealso{
  \code{\link{GenAlg}},
  \code{\link{GenAlg-tools}},
  \code{\link{maha}}.
}
\examples{
showClass("GenAlg")
}
\keyword{classes}
\keyword{optimize}
\keyword{classif}

