\name{ecoCusum}
\alias{ecoCusum}
\title{Economic design for Cusum control chart}
\usage{
  ecoCusum(h = seq(0.1, 2, by = 0.1), H = seq(0.4, 0.7, by
  = 0.01), n = 3:7, delta = 2, lambda = 0.01, P0 = NULL, P1
  = NULL, C0 = NULL, C1 = NULL, Cr = 20, Cf = 10, T0 = 0,
  Tc = 0.1, Tf = 0.1, Tr = 0.2, a = 0.5, b = 0.1, d1 = 1,
  d2 = 1, sided = "one", contour.plot = FALSE, call.print =
  TRUE, ...)
}
\arguments{
  \item{h}{sampling interval}

  \item{H}{decision interval}

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

  \item{sided}{distinguish between one-, two-sided and
  Crosier's modified two-sided CUSUM scheme by choosing
  "one", "two", and "Crosier", respectively. See details in
  \code{\link[spc]{xcusum.arl}}}

  \item{contour.plot}{a logical value indicating wether a
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
  Calculate the optimum parameters of n(sample size),
  h(sampling interval), k(reference value) and H(decision
  interval) for Econimic Design of Cusum control chart. For
  more information about the reference value see 'Details'.
}
\details{
  There is strong numerical and theoretical evidence that
  for given L1, the value of L0 approaches its maximum when
  k(reference value) is chosen mid-way the between AQL and
  the RQL: $k = mu0 + 0.5*delta*sigma (Appl. Statist.(1974)
  23, No. 3, p. 420). For this reason we treat k as a
  constant value and optimize n, h and H. For cost
  parameters P0, P1 and C0, C1, only one pair is needed. If
  P0 and P1 are given, they will be used first, else C0 and
  C1 will be used.
}
\examples{
#Taylor (1968). Technoinetrics, 10, p427 Table3, row 1-4,14
y=ecoCusum(P0=150,P1=50,n=3:7,Cr=30,h=seq(1.3,1.5,by=0.01),d1=0,d2=0)
summary(y)
contour(y)
#ecoCusum(P0=150,P1=50,n=3:7,lambda=0.05,Cr=30,d1=0,d2=0)
#ecoCusum(P0=150,P1=140,n=3:7,h=seq(4,5,by=0.01),Cr=30,d1=0,d2=0)
#ecoCusum(P0=2000,P1=1000,n=4:8,Cr=30,h=seq(0.3,0.6,by=0.01),H=seq(0.8,1,by=0.01),d1=0,d2=0)
#ecoCusum(P0=150,P1=50,Cr=30,delta=0.5,h=seq(2.5,3,by=0.01),n=25:35,H=seq(0.2,0.4,by=0.01),d1=0,d2=0)

#Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
ecoCusum(h=seq(0.6,1,by=.01),H=seq(.1,1,by=.1),n=3:8,lambda=.05,
P0=110,P1=10,Cr=25,Cf=50,Tr=0,Tf=0,Tc=1,T0=.0167,a=1)
}
\references{
  Lorenzen and Vance (1986). The economic design of control
  charts: a unified approach, Technometrics, 28. 3-10.

  Taylor (1968). The economic design of cumulative sum
  control charts, Technoinetrics, 10 479-488.
}
\seealso{
  \code{\link{ecoXbar}}, \code{\link{ecoEwma}},
  \code{\link[spc]{xcusum.arl}}
}

