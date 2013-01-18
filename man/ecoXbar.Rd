\name{ecoXbar}
\alias{echXbar}
\alias{ecoXbar}
\title{Economic design for the X-bar control chart}
\usage{
  ecoXbar(h, L, n, delta = 2, lambda = 0.05, P0 = NULL,
    P1 = NULL, C0 = NULL, C1 = NULL, Cr = 25, Cf = 50,
    T0 = 0.0167, Tc = 1, Tf = 0, Tr = 0, a = 1, b = 0.1,
    d1 = 1, d2 = 1, nlevels = 30, sided = "two",
    par = NULL, contour.plot = FALSE, call.print = TRUE,
    ...)

  echXbar(h, L, n, delta = 2, lambda = 0.05, P0 = NULL,
    P1 = NULL, C0 = NULL, C1 = NULL, Cr = 25, Cf = 50,
    T0 = 0.0167, Tc = 1, Tf = 0, Tr = 0, a = 1, b = 0.1,
    d1 = 1, d2 = 1, sided = "two")
}
\arguments{
  \item{h}{sampling interval. It can be a numeric vector or
  left undefined. See 'Details'}

  \item{L}{number of standard deviations from control
  limits to center line. It can be a numeric vector or left
  undefined. See 'Details'}

  \item{n}{sample size. It can be an integer vector or left
  undefined. See 'Details'}

  \item{lambda}{we assume the in-control time follows an
  exponential distribution with mean 1/lambda. Default
  value is 0.05.}

  \item{delta}{shift in process mean in standard deviation
  units when assignable cause occurs (delta = (mu1 -
  mu0)/sigma), where sigma is the standard deviation of
  observations; mu0 is the in-control process mean; mu1 is
  the out-of-control process mean. Default value is 2.}

  \item{P0}{profit per hour earned by the process operating
  in control. See 'Details'.}

  \item{P1}{profit per hour earned by the process operating
  out of control(P0 > P1).}

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

  \item{sided}{distinguish between one- and two-sided X-bar
  chart by choosing ``one'' or ``two'' respectively.  When
  \code{sided = ``one''}, \code{delta > 0} means the
  control chart for detecting a positive shift, and vice
  versa. Default is ``two''.}

  \item{par}{initial values for the parameters to be
  optimized over. It can be a vector of length 2 or 3. See
  'Details'}

  \item{contour.plot}{a logical value indicating whether a
  contour plot should be drawn. Default is FALSE. Only
  works when the parameters h, L and n are all specified.}

  \item{call.print}{a logical value indicating whether the
  "call" should be printed on the contour plot. Default is
  TRUE}

  \item{...}{other arguments to be passed to \code{optim}
  function.}
}
\value{
  The \code{ecoXbar} function returns an object of class
  "edcc", which is a list of elements \code{optimum},
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
  samples. The \code{echXbar} function returns the
  calculated ECH value only.
}
\description{
  Calculate the optimum parameters, n(sample size),
  h(sampling interval) and L(number of s.d. from control
  limits to center line) for Economic Design of the X-bar
  control chart .
}
\details{
  When parameter \code{par} is specified, optimization
  algorithms would be used as default. \code{par} can be
  specified as: \code{par = c(h, L)} where \code{h} and
  \code{L} are the initial values of smapling interval and
  control limit when \code{n} is specified; or \code{par =
  c(h, L, n)}. Good inital values may lead to good optimum
  results.

  When parameters \code{h}, \code{L}, \code{n} are all
  undefined, \code{ecoXbar} function will try to find the
  global optimum point to minimize the ECH (Expected Cost
  per Hour) using optimization algorithms (\code{optim}
  function), but in this case \code{n} would not be
  integer. It is usually helpful for the experimenter to
  find the region where the optimum point may exist
  quickly. When \code{h} and \code{L} are undefined but n
  is given as an integer vector, \code{ecoXbar} function
  will try to find the optimum point for each \code{n}
  value using optimization algorithms. When \code{h},
  \code{L} and \code{n} are all given, \code{ecoXbar}
  function will use the "grid method" to calculate the
  optimum point, that is ECH for all the combinations of
  the parameters will be calculated. The "grid method" way
  is much slower than using optimization algorithms, but it
  would be a good choice when optimization algorithms fail
  to work well.

  For cost parameters either {P0, P1} or {C0, C1} is
  needed.  If P0 and P1 are given, they will be used first,
  else C0 and C1 will be used.  For economic design of the
  X-bar chart, when \code{d1} and \code{d2} are both 1,
  only if the difference between P0 and P1 keeps the same,
  the results are identical. If the difference between C0
  and C1 keeps the same, the optimum parameters are almost
  the same but the ECH(Expected Cost per Hour) values will
  change.

  \code{echXbar} is used to calculate the ECH (Expected
  Cost per Hour) for one given design point.
}
\examples{
# Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
# In control profit per hour is 110, out of control profit per hour is 10
ecoXbar(P0=110,P1=10)
# In control profit per hour is 150, out of control profit per hour
# is 50, the result is identical with the previous one, because the
#difference between P0 and P1 are the same
ecoXbar(P0=150,P1=50)
# In control cost per hour is 0, out of control cost per hour is 100.
# The result is the same with the previous one
ecoXbar(C0=0,C1=100)
# The optimum parameters are the same with the previous one,
# but Cost values are different. See 'details'
ecoXbar(C0=10,C1=110)
# Based on the global optimum above, we specify the range of the
# parameters like this
x <- ecoXbar(h=seq(0.7,0.9,by=.01),L=seq(2.8,3.3,by=.01),n=4:6,P0=110,
P1=10,nlevels=50,contour.plot=TRUE)
x
# Modify the contour plot
contour(x,nlevels=20,lty=2,col=2,call.print=FALSE)
# update the parameters
update(x,P0=NULL,P1=NULL,C0=10,C1=110)
}
\references{
  Weicheng Zhu, Changsoon Park (2013), {edcc: An R Package
  for the Economic Design of the Control Chart}.
  \emph{Journal of Statistical Software}, 52(9), 1-24.
  \url{http://www.jstatsoft.org/v52/i09/}

  Douglas (2009). \emph{Statistical quality control: a
  modern introduction}, sixth edition, 463-471.

  Lorenzen and Vance (1986). The economic design of control
  charts: a unified approach, \emph{Technometrics}, 28.
  3-10.
}
\seealso{
  \code{\link{ecoCusum}}, \code{\link{ecoEwma}},
  \code{\link{contour}}, \code{\link[stats]{optim}},
  \code{\link{update.edcc}}
}

