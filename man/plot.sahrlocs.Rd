\name{plot.sahrlocs}
\alias{plot.sahrlocs}
\alias{print.plotsahr}
\title{Exploratory Analysis of Habitat Selection}
\description{
  \code{plot.sahrlocs} applies the function \code{widesII} or
  \code{widesIII} for each  variable in the object of class
  \code{sahrlocs}, and the results are stored in a list.  Then, the
  function \code{plot.wi} is applied to each component of
  the list.  This allows to investigate habitat selection by animals at
  several scales for design II and III data.
}
\usage{
\method{plot}{sahrlocs}(x, ani = names(x$hr), var = names(x$sa),
              type = c("hr.in.sa", "locs.in.hr", "locs.in.sa"),
              ncla = 4, ylog = FALSE, caxis = 0.7, clab = 0.7,
              errbar = c("SE", "CI"), alpha = 0.05, draw = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{sahrlocs}}
  \item{ani}{a character vector.  This vector contains the names
    of the animals in \code{x} for which habitat selection should be
    displayed.  At least two animals are required}
  \item{var}{a character vector.  This vector contains the names
    of the variables in \code{x} for which habitat selection should be
    displayed}
  \item{type}{a character string. Type of habitat selection that should
    be investigated. 
    If \code{"hr.in.sa"}, the selection of the home-range within the study
    area is displayed.  If \code{"locs.in.sa"}, the selection of the
    relocations within the study area is displayed. If
    \code{"locs.in.hr"}, the selection of the relocations within the home
    range is displayed}
  \item{ncla}{numeric variables are converted to factors.  This
    parameter controls the number of classes of these factors}
  \item{ylog}{logical.  If \code{TRUE}, the selection ratios are plotted
    on a log scale}
  \item{caxis}{character size on axes to be passed to
    \code{par("cex.axis")}}
  \item{clab}{character size of axes labels to be passed to
          \code{par("cex.lab")}}
  \item{errbar}{a character string.  Type of error bars: either
    \code{"CI"} for confidence intervals or \code{"SE"} for standard
    errors}
  \item{alpha}{the alpha-level for the tests}
  \item{draw}{logical. If \code{FALSE}, no plot is drawn}
  \item{\dots}{further arguments to be passed to the function
    \code{plot.wi}}
}
\value{
  The function returns a list of objects of class \code{wiII} or
  \code{wiIII} (one component per animal).
}
\references{
  Manly B.F.J., McDonald L.L., Thomas, D.L., McDonald, T.L. & Erickson,
  W.P. (2003) \emph{Resource selection by animals - Statistical design and
    Analysis for field studies. Second edition}.
  London: Kluwer academic publishers.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{widesII}} and \code{\link{widesIII}} for further
  information on objects of class \code{wiII} and \code{wiIII},
  \code{\link{as.sahrlocs}} for further information on objects of class
  \code{sahrlocs}.}

\examples{
data(puechabon)
sahr <- puechabon$sahr
toto <- plot(sahr)






### Note that the wild boars named Brock and Calou
### have only one herbaceous cover class available
### (the second one), and they use it exclusively.
### So they have identical
### selection ratios (that's why the curve of Brock
### does not appear: it is hidden behind the curve
### of the boar named Calou).

toto
toto$Elevation

}
\keyword{hplot}
\keyword{math}
