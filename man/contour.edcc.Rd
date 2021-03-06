\name{contour.edcc}
\alias{contour.edcc}
\title{Contour plot of an "edcc" class object}
\usage{
  \method{contour}{edcc} (x, call.print = TRUE, ...)
}
\arguments{
  \item{x}{an object of "edcc" class}

  \item{call.print}{a logical value indicating whether the
  the R command should be printed on the contour plot.
  Default is TRUE}

  \item{...}{arguments to be passed to contour plot, see
  \code{\link[graphics]{contour}} for details}
}
\value{
  a contour plot
}
\description{
  contour plot of an "edcc" class object
}
\details{
  S3 method of contour plot for "edcc" class object
}
\examples{
x <- ecoXbar(h=seq(0.7,0.9,by=.01),L=seq(2.8,3.3,by=.01),n=4:6,P0=110,
P1=10,nlevels=50,contour.plot=TRUE)
contour(x,nlevels=20,lty=2,col=2,call.print=FALSE)
}
\seealso{
  \code{\link{ecoXbar}}, \code{\link{ecoCusum}},
  \code{\link{ecoEwma}}, \code{\link{update.edcc}},
  \code{\link[graphics]{contour}}
}

