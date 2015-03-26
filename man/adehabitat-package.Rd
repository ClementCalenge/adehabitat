\name{adehabitat-package}
\alias{adehabitat-package}
\alias{adehabitat}
\docType{package}
\title{
  adehabitat: a Package for the Analysis of the Space Use by Animals
}
\description{
  The package adehabitat has been designed to allow the analysis of the
  space use by animals.  This package is organised in four main parts:
  (i) management of raster maps, (ii) habitat selection/ecological niche
  analysis,  (iii) home range estimation and (iv) analysis of animals
  trajectories.  These four parts interact with each other to facilitate the
  analysis.  Note that adehabitat strongly relies on the package ade4,
  which provides numerous functions for the analysis of multivariate
  data.
}
\details{
  \tabular{ll}{
    Package: \tab adehabitat\cr
    Type: \tab Package\cr
    Version: \tab 1.8.3\cr
    Date: \tab 2009-05-04\cr
    License: \tab GPL version 2 or newer\cr
  }

  The four parts of the packages are described more in detail below:
  \itemize{
    \item{Management of raster maps}{However, it provides numerous
      functions allowing to import and export raster maps from/to Geographic
      information system, to compute buffers around points or lines, to
      identify the value of environmental variables at given spatial
      locations, to count the number of points of a pattern in each pixel of
      a map, etc.  Two basic object classes allow to manage ratser maps with
      adehabitat: the class \code{asc} is intended to store basic raster
      map (see \code{help(import.asc)}), whereas the class \code{kasc} is
      intended to store multi-layer maps (all covering the same area with
      the same resolution, see \code{help(as.kasc)}).  For additional
      information this part of the package, see the tutorial available in
      the package.  Type \code{demo(rastermaps)} for a demonstration of the
      package capabilities.  Note that the package sp also provides many
      interesting functions to manage raster maps, and adehabitat provide
      functions of conversion to the classes of the package (see
      \code{help(kasc2spixdf)}).}
    
    \item{Habitat selection/ecological niche analysis}{many methods have
      been included in the package to render statistical methods allowing
      the analysis of habitat selection available to wildlife
      ecologists. Many of them are factor analyses of the niche or of the
      habitat selection, but other methods are also available. They 
      include the selection ratios (see \code{help(wi)}), the Ecological
      niche factor analysis (see \code{help(enfa)}), the Mahalanobis
      distances (see \code{help(mahasuhab)}) and their factorial
      decomposition (the MADIFA, see \code{help(madifa)}) or the algorithm
      DOMAIN (see \code{help(domain)}). Other common methods, such as the
      resource selection functions can also be used with the rest of the R
      environment.  Note that the package also include functions allowing
      the analysis of habitat selection using radio-tracking data, such as
      the compositional analysis (see \code{help(compana)}), the
      eigenanalysis of selection ratios (see \code{help(eisera)}) or the
      K-select analysis (see \code{help(kselect)}). An overview of these
      methods is available by typing \code{demo(nichehs)}.}
    
    \item{home range estimation}{many methods are available to estimate
      the home range of animals relocated using radio-tracking data. Common
      methods of estimation are available, such as the Minimum convex
      polygon (see \code{help(mcp)}), the kernel estimation of the
      utilization distribution (see \code{help(kernelUD)}), the cluster home
      range (see \code{help(clusthr)}) or the nearest neighbour convex hull
      (see \code{help(NNCH)}). Note that Paolo Cavallini has designed a
      website dedicated to the analysis of space use by animals, which
      contain a wiki page, a tutorial for the home range estimation using R
      and adehabitat and a forum (URL:
      http://www.faunalia.it/animov/index.php). Several methods of these
      part of the package have been included following discussions that
      arose on this forum (especially, the nearest neighbour convex hull and
      the brownian bridge kernel). Type \code{demo(homerange)} for examples
      of use of these functions.}
    
    \item{The analysis of animals' trajectories}{This part is the most recent
      one in the package. A new class designed to store animals'
      trajectories has been included in the package, the class
      \code{ltraj} (see \code{help(as.ltraj)}). Two types of
      trajectories can be handled with adehabitat: for trajectories of
      type I, the time is not recorded for the relocations (e.g. the
      sampling of the tracks of an animal in the snow). For trajectories
      of type II, the time has been recorded during sampling
      (e.g. radio-tracking, GPS, Argos monitoring). Many descriptive
      parameters are automatically computed (relative or turning angles,
      distance between successive relocations, mean squared
      displacement). Many functions allow the management and analysis of 
      trajectories, through the analysis of these parameters (e.g. tests of
      independence, see \code{help(wawotest,indmove)}, first passage time,
      see \code{fpt}). The rediscretisation of trajectories of type I is also
      possible (\code{help(redisltraj)}). Many graphical functions are
      available for the  exploration of the trajectory properties
      (\code{plot,plotltr,sliwinltr}), etc.  A new partitioning algorithm
      has been added (but it is still under research) to partition animals
      trajectories into segment with homogeneous properties (see
      \code{modpartltraj}).  Further details can be found on
      the help page of the function \code{as.ltraj}. For a demonstration,
      type \code{demo(managltraj)} or \code{demo(analysisltraj)}.}
  } 
}
\author{
Clement Calenge, with contributions from Mathieu Basille, Stephane Dray,
Manuela Royer and Scott Fortmann-Roe

Maintainer: Clement Calenge <clement.calenge@oncfs.gouv.fr>}

\references{
  Calenge, C. (2006) The package adehabitat for the R software: a tool
  for the analysis of space and habitat use by animals. Ecological
  Modelling, 197, 516-519
  
}
\keyword{ package }
\seealso{
\code{\link[ade4:ade4-package]{ade4}}
}
\examples{
\dontrun{

## For examples of use of mapping capabilities
demo(rastermaps)

## For examples of use of functions for
## habitat selection and niche analysis
demo(nichehs)

## For example of home range estimation
demo(homerange)

## For example of trajectory management and analysis
demo(managltraj)
demo(analysisltraj)


}
}
