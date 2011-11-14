\name{ecoEwma}
\alias{ecoEwma}
\title{Economic design for EWMA control chart}
\usage{
  ecoEwma(h = seq(0.7, 1, by = 0.1), w = seq(0.7, 1, by =
  0.1), k = seq(2, 4, by = 0.1), n = 4:8, delta = 2, lambda
  = 0.05, P0 = NULL, P1 = NULL, C0 = NULL, C1 = NULL, Cr =
  25, Cf = 10, T0 = 0.0167, Tc = 1, Tf = 0, Tr = 0, a = 1,
  b = 0.1, d1 = 1, d2 = 1, sided = "two", contour.plot =
  FALSE, call.print = TRUE, ...)
}
\arguments{
  \item{h}{sampling interval}

  \item{w}{the weight value given to the latest sample}

  \item{k}{control limit coefficient}

  \item{n}{sample size}

  \item{delta}{critical value: the extent of the shift when
  assignable cause occurs. delta = |(mu1 - mu0)/sigma|.
  Default value is 2.}

  \item{lambda}{we assume the in-control time follows a
  exponential distribution with mean 1/lambda. Default
  value is 0.05.}

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

  \item{sided}{distinguish between one- and two-sided EWMA
  control chart by choosing "one" and "two", respectively.
  See details in \code{\link[spc]{xewma.arl}}}

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
  h(sampling interval), w(weight to the present sample) and
  k(number of s.d. from control limits to center line) for
  econimic Design of EWMA control chart .
}
\details{
  For cost parameters P0, P1 and C0, C1, only one pair is
  needed. If P0 and P1 are given, they will be used first,
  else C0 and C1 will be used.
}
\examples{
#Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
x = ecoEwma(n=4:8,P0=110,P1=10,Cf=50,contour.plot=TRUE)
summary(x)
contour(x)
#ecoEwma(P0=150,P1=50,Cr=30,delta=0.5,h=seq(2.5,3,by=0.01),n=25:35)
}
\seealso{
  \code{\link{ecoXbar}}, \code{\link{ecoCusum}}
  \code{\link[spc]{xewma.arl}}
}

