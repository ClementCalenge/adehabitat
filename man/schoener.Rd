\name{schoener}
\alias{schoener}
\alias{schoener.rtest}
\title{Compute Schoener's ratio }
\description{
  \code{schoener} computes the Schoener's ratio on radio-tracking data.\cr
  \code{schoener.rtest} performs a randomization test of the equality of
  the Schoener's ratio to 2\cr 
}
\usage{
schoener(tr, keep, byburst = TRUE)
schoener.rtest(tr, keep, byburst = TRUE, nrep = 500)
}
\arguments{
  \item{tr}{an object of class \code{traj} }
  \item{keep}{a vector of length 2, giving the lower and the upper bound
    of the time interval for which a pair of relocations is considered
    in the computation of t\^2 (see details).  These values are given in
    seconds.}
  \item{byburst}{logical.  If \code{TRUE}, the Schoener's ratio is
    computed by burst.  If \code{FALSE}, the ratio is computed by
    animal.}
  \item{nrep}{the number of randomisations of the test. }
  
}
\details{
  The Schoener's ratio is a measure of time-autocorrelation in the data.
  This ratio is computed as the squared mean distance between "neighbour"
  relocations (t\^2) divided by the squared mean distance between the
  relocations and their barycenter (r\^2).  The theoretical value of this
  ratio under the hypothesis of independance of the relocations is 2.\cr

  Swihart and Slade (1985) consider as neighbour two successive
  relocations.  However, the Schoener's ratio computed in this way makes
  sense biologically only if the relocations are equally spaced in
  time.  However, as indicated by these authors, "such a data set
  probably is the exception rather than the rule because many problems
  may arise in taking a locational reading at a specified time".\cr

  In this function, we define as "neighbour" two relocations (not
  necessarily successive relocations) separated by a time interval
  comprised within the bounds specified in the vector \code{keep} (in
  seconds).  For example, if \code{keep = c(60, 300)}, all relocations
  separated by a time interval comprised between 1 and 5 minutes are
  considered in the computation.  Thus, the total number of pairs of
  relocations m taken into account in the computation may be larger than
  the number of relocations n (m can be at most equal to n*(n-1)/2).
   
}
\value{
  returns an object of class \code{schoener}.
}
\references{
  Schoener, T.W. (1981) An empirically based estimate of home
  range. \emph{Theoretical Population Biology}, \bold{20}, 281--325.

  Swihart, R.K. and Slade, N.A. (1985) Testing for independence of
  observations in animal movements. \emph{Ecology}, \bold{66},
  1176--1184.
  
  Solow, A.R. (1989) A randomization test for independence of animal
  locations. \emph{Ecology}, \bold{70}, 1546--1549. 
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}

\examples{

data(puechcirc)
puechcirc <- ltraj2traj(puechcirc)
puechcirc$date[1:10]

## Relocations are taken every 10 minutes
## For example we consider relocations as
## neighbour when they are separated by a time
## interval comprised between 5 and 15 minutes
schoener(puechcirc, keep = c(5*60, 15*60))
\dontrun{
schoener.rtest(puechcirc, keep = c(5*60, 15*60))
}



}
\keyword{math}

