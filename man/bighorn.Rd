\name{bighorn}
\docType{data}
\alias{bighorn}
\title{ Radio-Tracking of Bighorn Sheeps}
\description{
This data set describes the habitat use and availability for 6 bighorn
sheeps monitored by radio-tracking (Arnett et al. 1989, in Manly et al.,
2003, p. 67-74). 10 habitat types are considered.
}
\usage{
data(bighorn)
}
\details{
  The object \code{bighorn} is a list, with the following components:
  \describe{
  \item{\code{used}}{the number of resource units used by each animal (in rows)
    in each habitat category (in columns).}
  \item{\code{availTrue}}{the availability of each habitat category.}
  \item{\code{availEstimated}}{a sample of available resource units in each
    habitat category.}
  }
}
\references{
  Manly, B.F.J., McDonald, L.L., Thomas, D.L., McDonald, T.L. & Erickson,
  W.P. (2003) \emph{Resource selection by animals - Statistical design and
    Analysis for field studies. Second edition.}
  London: Kluwer academic publishers. 
}
\keyword{datasets}
