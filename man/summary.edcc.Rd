\name{summary.edcc}
\alias{summary.edcc}
\title{summary a "edcc" class object}
\usage{
  \method{summary}{edcc} (object, ...)
}
\arguments{
  \item{object}{an object of class "edcc"}

  \item{...}{not used in this function}
}
\value{
  A list containing the optimum reslut and a data frame is
  returned.
}
\description{
  summary a "edcc" class object
}
\details{
  Both the optimum result and a data frame(cost.frame) are
  returned. Each row of the cost.frame instead for an
  optimum result corresponding to the specified n value.
}
\examples{
x <- ecoXbar(P0=100,P1=0,nlevels=50)
summary(x)
}

