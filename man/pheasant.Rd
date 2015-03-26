\name{pheasant}
\docType{data}
\alias{pheasant}
\title{Radio-Tracking of Pheasants}
\description{
This data set describes the use and availability of 5 habitat types for
13 pheasants monitored using radio-tracking. 
}
\usage{data(pheasant)}
\format{
  This list has three components:
  \describe{
    \item{\code{studyarea}}{a data frame giving the proportion of each
      habitat type (columns) on the study area.  These habitat types are
      \code{Scrub}, \code{Broadleaf}, \code{Coniferous},
      \code{Grassland} and \code{Crop}.  These proportions are
      repeated by rows, for all animals (rows)}
    \item{\code{mcp}}{a data frame giving the proportion of
      each habitat type (columns) in the home range of each animal
      (rows)}
    \item{\code{locs}}{a data frame giving the proportion of relocations
      of each animal (rows) reported in 3 of the 5 habitat types
      (columns). Indeed, \code{Coniferous} and \code{Crops} were not
      used by most of the animals.}
    }
}
\source{
  Aebischer, N. J., Robertson, P. A. and Kenward, R. E. (1993)
  Compositional analysis of habitat use from animal radiotracking
  data. \emph{Ecology}, \bold{74}, 1313--1325.
}
\keyword{datasets}
