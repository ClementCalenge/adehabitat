\name{biv.test}
\alias{biv.plot}
\alias{biv.test}
\title{Bivariate Test}
\description{
  \code{biv.plot} displays a bivariate plot.
  \code{biv.test} displays the results of a bivariate randomisation test.
}
\usage{
biv.plot(dfxy, br = 10, points = TRUE, density = TRUE, 
         kernel = TRUE, o.include = FALSE, pch, cex, col, h, sub, 
         side = c("top", "bottom", "none"), \dots)
biv.test(dfxy, point, br = 10, points = TRUE, density = TRUE, 
         kernel = TRUE, o.include = FALSE, pch, cex, col, Pcol, h, sub, 
         side = c("top", "bottom", "none"), \dots)
}
\arguments{
  \item{dfxy}{a data frame with N lines (couples of values) and two
    columns}
  \item{br}{a parameter used to define the numbers of breaks of the
    histograms.  A larger value leads to a larger number of breaks}
  \item{points}{logical.  Whether the points should be displayed}
  \item{density}{logical.  Whether the kernel density estimation should be 
    displayed for the marginal histograms}
  \item{kernel}{logical.  Whether the kernel density estimation should be 
    displayed for the bivariate plot}
  \item{o.include}{logical.  If \code{TRUE}, the origin is included in
    the plot}
  \item{pch}{plotting "character", i.e., symbol to use for the points.  (see 
    \code{?points})}
  \item{cex}{character expansion for the points}
  \item{col}{color code or name for the points, see \code{?par}}
  \item{h}{vector of bandwidths for x and y directions, used in the
    function \code{kde2d} of the package \code{MASS}.  Defaults to normal
    reference bandwidth (see \code{?kde2d})}
  \item{sub}{a character string to be inserted in the plot as a title}
  \item{side}{if \code{"top"}, the x and y scales of the grid are upside, 
    if \code{"bottom"} they are downside, if \code{"none"} no legend}
  \item{point}{a vector of length 2, representing the observation to be 
    compared with the simulated values of the randomisation test}
  \item{Pcol}{color code or name for the observed point}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  \code{biv.test} is used to display the results of a bivariate
  randomisation test.  An example of use of the function is provided in
  the function \code{niche.test}.
  
  The x-axis of the main window corresponds to the first column of
  \code{dfxy}; the y-axis corresponds to the second column.  Kernel
  density is estimated to indicate the contours of the distribution of
  randomised values.  The two marginal histograms correspond to the
  univariate tests on each axis, for which the p-values are computed with
  \code{as.randtest} (package ade4, one-sided tests).
}
\author{Mathieu Basille \email{basille@ase-research.org}}
\section{Warning}{
  \code{biv.plot} and \code{biv.test} uses the function \code{kde2d} of the 
    package \code{MASS}.
}
\seealso{
  \code{\link[ade4]{as.randtest}} (package ade4), \code{\link{niche.test}}
}
\examples{

x = rnorm(1000,2)
y = 2*x+rnorm(1000,2)
dfxy = data.frame(x, y)

biv.plot(dfxy)
biv.plot(dfxy, points=FALSE, col="lightblue", br=20)

p = c(3, 4)
biv.test(dfxy, p)
biv.test(dfxy, p, points=FALSE, Pcol="darkred", col="lightblue", br=20)

}
\keyword{multivariate}
