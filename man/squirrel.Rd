\name{squirrel}
\docType{data}
\alias{squirrel}
\title{Radio-Tracking Data of Squirrels}
\description{
This data set describes the use and availability of 5 habitat types for
17 squirrels monitored using radio-tracking.  See also the dataset
\code{squirreloc}.
}
\usage{data(squirrel)}
\format{
  This list has three components: \describe{
    \item{\code{studyarea}}{a data frame giving the proportion of each
      habitat type (columns) on the study area. These proportions are
      repeated by rows, for all animals}
    \item{\code{mcp}}{a data frame giving the proportion of
      each habitat type (columns) in the home range of each animal (rows)}
    \item{\code{locs}}{a data frame giving the proportion of relocations
      of each animal (rows) reported in each habitat type (columns).}
    }
}
\source{
  Aebischer, N. J., Robertson, P. A. and Kenward, R. E. (1993)
  Compositional analysis of habitat use from animal radiotracking
  data. \emph{Ecology}, \bold{74}, 1313--1325.
}
\keyword{datasets}
