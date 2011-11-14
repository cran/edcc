\name{contour.edcc}
\alias{contour.edcc}
\title{contour plot of "edcc" class}
\usage{
  \method{contour}{edcc} (x, ...)
}
\arguments{
  \item{x}{an object of "edcc" class}

  \item{...}{arguments to be passed to contour plot, not
  supported now}
}
\value{
  a contour plot
}
\description{
  contour plot of "edcc" class
}
\details{
  contour plot for "edcc" class object
}
\examples{
z=ecoXbar(P0=100,P1=0,nlevels=50)
contour(z)
}

