\name{print.edcc}
\alias{print.edcc}
\title{Print a "edcc" class object}
\usage{
  \method{print}{edcc} (x, ...)
}
\arguments{
  \item{x}{an object of class "edcc"}

  \item{...}{not used in this function}
}
\value{
  Invisible, the object itself.
}
\description{
  Print a "edcc" class object
}
\details{
  Print the optimum parameters of a "edcc" object
}
\examples{
x <- ecoXbar(P0=100,P1=0,nlevels=50)
print(x)
}

