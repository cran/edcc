\name{ecoXbar}
\alias{ecoXbar}
\title{Economic design for X-bar chart}
\usage{
  ecoXbar(h = seq(0.1, 1, by = 0.01), L = seq(2, 4.5, by =
  0.01), n = 1:15, lambda = 0.05, delta = 2, P0 = NULL, P1
  = NULL, C0 = NULL, C1 = NULL, Cr = 25, Cf = 50, T0 =
  0.0167, Tc = 1, Tf = 0, Tr = 0, a = 1, b = 0.1, d1 = 1,
  d2 = 1, contour.plot = FALSE, call.print = TRUE, ...)
}
\arguments{
  \item{h}{sampling interval}

  \item{L}{number of standard deviations from control
  limits to center line.}

  \item{n}{sample size.}

  \item{lambda}{we assume the in-control time follows a
  exponential distribution with mean 1/lambda. Default
  value is 0.05.}

  \item{delta}{critical value: the extent of the shift when
  assignable cause occurs. delta = |(mu1 - mu0)/sigma|.
  Default value is 2.}

  \item{P0}{profit per hour earned by the process operating
  in control. See 'Details'.}

  \item{P1}{profit per hour earned by the process operating
  out of control}

  \item{C0}{cost per hour due to nonconformities produced
  while the process is in control.}

  \item{C1}{cost per hour due to nonconformities produced
  while the process is out of control.(C1 > C0)}

  \item{Cr}{cost for searching and repairing the assignable
  cause, including any downtime.}

  \item{Cf}{cost per false alarm, including the cost of
  searching for the cause and the cost of downtime if
  production ceases during search.}

  \item{T0}{time to sample and chart one item.}

  \item{Tc}{expected time to discover the assignable
  cause.}

  \item{Tf}{expected search time when false alarm occures.}

  \item{Tr}{expected time to repair the process.}

  \item{a}{fixed cost per sample.}

  \item{b}{cost per unit sampled.}

  \item{d1}{flag for whether production continues during
  searches (1-yes, 0-no). Default value is 1.}

  \item{d2}{flag for whether production continues during
  repairs (1-yes, 0-no). Default value is 1.}

  \item{contour.plot}{a logical value indicating whether a
  contour plot should be drawn. Default is FALSE.}

  \item{call.print}{a logical value indicating whether the
  "call" should be drawn on the contour plot. Default is
  TRUE}

  \item{...}{other arguments to be passed to contour
  function.}
}
\value{
  return the optimum parameters and the corresponding cost
  value
}
\description{
  Calculate the optimum parameters, n(sample size),
  h(sampling interval) and L(number of s.d. from control
  limits to center line) for Econimic Design of X-bar
  control chart .
}
\details{
  For cost parameters P0, P1 and C0, C1, only one pair is
  needed. If P0 and P1 are given, they will be used first,
  else C0 and C1 will be used. For economic design of x-bar
  chart, only if the difference between P0 and P1 keeps the
  same, the results are identical. If the difference
  between C0 and C1 keeps the same, the optimum parameters
  are identical but the cost values will change.
}
\examples{
# Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
# In control profit per hour is 110, out of control profit per hour is 10
x=ecoXbar(P0=110,P1=10,contour.plot=TRUE,nlevels=50)
summary(x)
contour(x)
# In control profit per hour is 150, out of control profit per hour
# is 50, the result is identical with the previous one, because the
#difference between P0 and P1 are the same
ecoXbar(P0=150,P1=50,contour.plot=TRUE,nlevels=50)
# In control cost per hour is 0, out of control cost per hour is 100.
#The result is the same with the previous one
ecoXbar(C0=0,C1=100,contour.plot=TRUE,nlevels=50)
# The optimum parameters are the same with the previous one,
# but Cost values are different. See 'details'
ecoXbar(C0=10,C1=110,contour.plot=TRUE,nlevels=50)
}
\references{
  Douglas (2009). Statistical quality control: a modern
  introduction, sixth edition, 463-471.

  Lorenzen and Vance (1986). The economic design of control
  charts: a unified approach, Technometrics, 28. 3-10.
}
\seealso{
  \code{\link{ecoCusum}}, \code{\link{ecoEwma}}
}

