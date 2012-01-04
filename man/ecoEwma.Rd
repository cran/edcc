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
  Parameter \code{w} should always be given, because the
  range of \code{w} is so restricted that optimization
  algorithms usually don't converge.  When parameters
  \code{h}, \code{k}, \code{n} are all undefined,
  \code{ecoEwma} function will try to find the global
  optimum point to minimize the ECH (Expected Cost per
  Hour) using optimization algorithms (\code{optim}
  function), but in this case \code{n} would not be
  integer. It is usually helpful for the experimenter to
  find the region where the optimum point may exist
  quickly. When \code{h} and \code{k} are undefined but
  \code{n} is given as an integer vector, \code{ecoEwma}
  function will try to find the optimum point for each
  \code{n} value using optimization algorithms. When
  \code{h}, \code{k} and \code{n} are all given,
  \code{ecoEwma} function will use a "robust" way to
  calculate the optimum point, that is ECH for all the
  combinations of the parameters will be calculated. The
  "robust" way is much slower than using optimization
  algorithms, but it would be a good choice when
  optimization algorithms fail to work well.

  For cost parameters P0, P1 and C0, C1, only one pair is
  needed. If P0 and P1 are given, they will be used first,
  else C0 and C1 will be used.
}
\examples{
#Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
## Set w from 0.1 to 1 by 0.1 to catch the trend.
y <- list()
w <- seq(0.1,0.9,by=0.1)
for(i in 1:length(w))
y[[i]] <- ecoEwma(w=w[i],P0=110,P1=10,Cf=50,contour.plot=TRUE)$optimum
y <- do.call("rbind",y)
#yy = ecoEwma( h=seq(.7,1,by=.01), w=seq(0.8,1,by=.01),k=seq(2.9,3.3,by=0.01),n=4:5,P0=110,P1=10,Cf=50,contour.plot=TRUE)
##$optimum
##Optimum h Optimum w Optimum k     n       ECH
##0.81000   0.95000   2.99000   5.00000  10.36482
#contour(yy)
}
\seealso{
  \code{\link{ecoXbar}}, \code{\link{ecoCusum}},
  \code{\link[spc]{xewma.arl}}, \code{\link{update.edcc}},
  \code{\link[stats]{optim}},\code{\link{contour}}
}

