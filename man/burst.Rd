\name{burst}
\alias{burst}
\alias{id}
\alias{burst<-}
\alias{id<-}
\title{ID and Bursts of an Object of Class ltraj}
\description{
  Functions to get or set the attribute \code{"id"} or \code{"burst"} of
  the components of an object of class \code{ltraj}.
}
\usage{
burst(ltraj)
burst(ltraj) <- value
id(ltraj)
id(ltraj) <- value
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{value}{a character vector of up to the same length as
    \code{ltraj}}
}
\details{
     The functions \code{id} and \code{burst} are accessor
     functions, and \code{id<-} and \code{burst<-} are replacement
     function.
}
\value{
  For \code{id} and \code{burst}, a character vector of the same length
  as \code{ltraj}.\cr
  
  For \code{id<-} and \code{burst<-}, the updated object.  (Note that
  the value of \code{burst(x) <- value} is that of the assignment,
  \code{value}, not the return value from the left-hand side.)
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link{ltraj}}, \code{\link{names}}}
\examples{
data(puechcirc)
puechcirc

## To see the ID and the burst
id(puechcirc)
burst(puechcirc)

## Change the burst
burst(puechcirc) <- c("glou", "toto", "titi")
puechcirc

burst(puechcirc)[2] <- "new name"
puechcirc

## Change the ID
id(puechcirc)[id(puechcirc)=="CH93"] <- "WILD BOAR"
puechcirc

}
\keyword{programming}
