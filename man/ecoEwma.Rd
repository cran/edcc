\name{ecoEwma}
\alias{echEwma}
\alias{ecoEwma}
\title{Economic design for EWMA control chart}
\usage{
  ecoEwma(h = seq(0.7, 1, by = 0.1),
    w = seq(0.7, 1, by = 0.1), k = seq(2, 4, by = 0.1),
    n = 4:8, delta = 2, lambda = 0.05, P0 = NULL,
    P1 = NULL, C0 = NULL, C1 = NULL, Cr = 25, Cf = 10,
    T0 = 0.0167, Tc = 1, Tf = 0, Tr = 0, a = 1, b = 0.1,
    d1 = 1, d2 = 1, nlevels = 30, sided = "two",
    par = NULL, contour.plot = FALSE, call.print = TRUE,
    ...)

  echEwma(h, w, k, n, delta = 2, lambda = 0.05, P0 = NULL,
    P1 = NULL, C0 = NULL, C1 = NULL, Cr = 25, Cf = 10,
    T0 = 0.0167, Tc = 1, Tf = 0, Tr = 0, a = 1, b = 0.1,
    d1 = 1, d2 = 1, sided = "two")
}
\arguments{
  \item{h}{sampling interval. It can be a numeric vector or
  left undefined. See 'Details'}

  \item{w}{the weight value between 0 and 1 given to the
  latest sample. It must be specified.}

  \item{k}{control limit coefficient. It can be a numeric
  vector or left undefined. See 'Details'}

  \item{n}{sample size. It can be an integer vector or left
  undefined. See 'Details'}

  \item{delta}{shift in process mean in standard deviation
  units when assignable cause occurs (delta = |mu1 -
  mu0|/sigma), where sigma is the standard deviation of
  observations; mu0 is the in-control process mean; mu1 is
  the out-of-control process mean. Default value is 2.}

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

  \item{Tf}{expected search time when false alarm occurs.}

  \item{Tr}{expected time to repair the process.}

  \item{a}{fixed cost per sample.}

  \item{b}{cost per unit sampled.}

  \item{d1}{flag for whether production continues during
  searches (1-yes, 0-no). Default value is 1.}

  \item{d2}{flag for whether production continues during
  repairs (1-yes, 0-no). Default value is 1.}

  \item{nlevels}{number of contour levels desired. Default
  value is 30. It works only when \code{contour.plot} is
  TRUE.}

  \item{sided}{distinguish between one- and two-sided EWMA
  control chart by choosing "one" and "two", respectively.
  See details in \code{\link[spc]{xewma.arl}}}

  \item{par}{initial values for the parameters to be
  optimized over. It can be a vector of length 2 or 3. See
  'Details'}

  \item{contour.plot}{a logical value indicating whether a
  contour plot should be drawn. Default is FALSE.}

  \item{call.print}{a logical value indicating whether the
  "call" should be printed on the contour plot. Default is
  TRUE}

  \item{...}{other arguments to be passed to contour
  function.}
}
\value{
  The function returns a list of elements \code{optimum},
  \code{cost.frame}, \code{FAR} and \code{ATS}.
  \code{optimum} is a vector with the optimum parameters
  and the corresponding ECH value; \code{cost.frame} is a
  dataframe with the optimum parameters and the
  corresponding ECH values for all given \code{n}(if
  \code{n} is not specified, \code{cost.frame} won't be
  returned); \code{FAR} indicates the false alarm rate
  during the in-control time, which is calculated as
  lambda*(average number of false alarm); \code{ATS}
  indicates the average time to signal after the occurrence
  of an assignable cause, calculated as h*ARL2 - tau, where
  tau is the expected time of occurrence of the assignable
  cause given it occurs between the i-th and (i+1)st
  samples
}
\description{
  Calculate the optimum parameters, n(sample size),
  h(sampling interval), w(weight to the present sample) and
  k(number of s.d. from control limits to center line) for
  economic Design of EWMA control chart .
}
\details{
  Parameter \code{w} should always be given, because the
  range of \code{w} is so restricted that optimization
  algorithms usually don't converge.

  When parameter \code{par} is specified, optimization
  algorithms would be used as default. \code{par} can be
  specified as: \code{par = c(h, k)} where \code{h} and
  \code{k} are the initial values of smapling interval and
  control limit when \code{n} is specified; or \code{par =
  c(h, k, n)}. Good inital values may lead to good optimum
  results.

  When parameters \code{h}, \code{k}, \code{n} are all
  undefined, \code{ecoEwma} function will try to find the
  global optimum point to minimize the ECH (Expected Cost
  per Hour) using optimization algorithms (\code{optim}
  function), but in this case \code{n} would not be
  integer. It is usually helpful for the experimenter to
  find the region where the optimum point may exist
  quickly. When \code{h} and \code{k} are undefined but
  \code{n} is given as an integer vector, \code{ecoEwma}
  function will try to find the optimum point for each
  \code{n} value using optimization algorithms. When
  \code{h}, \code{k} and \code{n} are all given,
  \code{ecoEwma} function will use a "grid method" way to
  calculate the optimum point, that is ECH for all the
  combinations of the parameters will be calculated. The
  "grid method" way is much slower than using optimization
  algorithms, but it would be a good choice when
  optimization algorithms fail to work well.

  For cost parameters either {P0, P1} or {C0, C1} is
  needed. If P0 and P1 are given, they will be used first,
  else C0 and C1 will be used.

  \code{echEwma} is used to calculate the ECH (Expected
  Cost per Hour) for a given design point
}
\examples{
#Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
## Set w from 0.1 to 1 by 0.1 to catch the trend.
ecoEwma(w=seq(0.1,1,by=0.1),P0=110,P1=10,Cf=50)
#yy = ecoEwma( h=seq(.7,1,by=.01), w=seq(0.8,1,by=.01),k=seq(2.9,3.3,by=0.01),n=4:5,P0=110,P1=10,Cf=50,contour.plot=TRUE)
##$optimum
##Optimum h Optimum w Optimum k     n       ECH
##0.81000   0.95000   2.99000   5.00000  10.36482
#contour(yy)
}
\references{
  Lorenzen and Vance (1986). The economic design of control
  charts: a unified approach, Technometrics, 28. 3-10.
}
\seealso{
  \code{\link{ecoXbar}}, \code{\link{ecoCusum}},
  \code{\link[spc]{xewma.arl}}, \code{\link{update.edcc}},
  \code{\link[stats]{optim}},\code{\link{contour}}
}

