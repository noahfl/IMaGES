\name{IMaGES}
\alias{IMaGES}
\encoding{UTF-8}
\concept{Independent Multi-sample Greedy Equivalence Search}
\concept{essential graph}
\concept{graphNEL}
\title{Estimate the DAG and Markov equivalence class of multiple datasets using
the IMaGES algorithm}
\description{
  Estimate the DAG using the IMaGES
  algorithm of Hanson et al (2009).
}
\usage{
new("IMaGES", matrices = list(matrix1, matrix2....), penalty=3, num.markovs=5)
}
\arguments{
  \item{matrices}{A list of (at least one) matrix object(s) to run in the algorithm}
  \item{penalty}{Node labels; by default, they are determined from the scoring
    object.}
  \item{num.markovs}{Number of graphs to store in the Markov equivalence class.}
}
\details{

  IMaGES (Independent Multi-sample Greedy Equivalence Search) is a score-based 
  algorithm that greedily maximizes a score function similar to the one used
  in the pcalg implementation of GES. It modifies the scoring by creating a
  global score across all datasets and uses this score to determine which 
  step from the individual datasets best represents all of the datasets. It 
  accomplishes this by using forward, backward, and turning steps as described
  below:
  \describe{
    \item{Forward phase}{In the forward phase, IMaGES moves through the space of
    essential graphs in steps that correspond to the addition of a single edge
    in the space of DAGs; the phase is aborted as soon as the score cannot be
    augmented any more.}
    \item{Backward phase}{In the backward phase, the algorithm performs moves
    that correspond to the removal of a single edge in the space of DAGs until
    the score cannot be augmented any more.}
    \item{Turning phase}{In the turning phase, the algorithm performs moves
    that correspond to the reversal of a single arrow in the space of DAGs
    until the score cannot be augmented any more.}
    
    
  }
    These operations are carried out and result in a global representative graph
    and a Markov Equivalence Class.
    
    TODO: finish

}
\value{
  \code{IMaGES} returns a list with the following two components:
  \item{essgraph}{An object of class \code{\linkS4class{EssGraph}} containing an
    estimate of the equivalence class of the underlying DAG.}
  \item{repr}{An object of a class derived from \code{\linkS4class{ParDAG}}
    containing a (random) representative of the estimated equivalence class.}
}
\references{
  D.M. Chickering (2002).  Optimal structure identification with greedy search.
  \emph{Journal of Machine Learning Research} \bold{3}, 507--554

  A. Hauser and P. Bühlmann (2012).  Characterization and greedy learning of
  interventional Markov equivalence classes of directed acyclic graphs.
  \emph{Journal of Machine Learning Research} \bold{13}, 2409--2464.

  P. Nandy, A. Hauser and M. Maathuis (2015).  Understanding consistency in 
  hybrid causal structure learning.  \emph{arXiv preprint} 1507.02608
  
  P. Spirtes, C.N. Glymour, and R. Scheines (2000).  \emph{Causation,
  Prediction, and Search}, MIT Press, Cambridge (MA).
}
\author{
  Noah Frazier-Logue (\email{n.frazier.logue@nyu.edu})
}
\seealso{
  \code{\link{plotIMGraph}}, \code{\link{plotAll}}, \code{\link{plotMarkovs}}
}
\examples{
## Load predefined data
data(gmG)

## Define the score (BIC)
score <- new("GaussL0penObsScore", gmG8$x)

## Estimate the essential graph
ges.fit <- ges(score)

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(ges.fit$essgraph, main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
} else { ## alternative:
  str(ges.fit, max=2)
  as(as(ges.fit$essgraph,"graphNEL"),"Matrix")
}
}
\keyword{models}
\keyword{graphs}
