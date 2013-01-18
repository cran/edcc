\name{ecoCusum}
\alias{echCusum}
\alias{ecoCusum}
\title{Economic design for the CUSUM control chart}
\usage{
  ecoCusum(h, H, n, delta = 2, lambda = 0.01, P0 = NULL,
    P1 = NULL, C0 = NULL, C1 = NULL, Cr = 20, Cf = 10,
    T0 = 0, Tc = 0.1, Tf = 0.1, Tr = 0.2, a = 0.5, b = 0.1,
    d1 = 1, d2 = 1, nlevels = 30, sided = "one",
    par = NULL, contour.plot = FALSE, call.print = TRUE,
    ...)

  echCusum(h, H, n, delta = 2, lambda = 0.01, P0 = NULL,
    P1 = NULL, C0 = NULL, C1 = NULL, Cr = 20, Cf = 10,
    T0 = 0, Tc = 0.1, Tf = 0.1, Tr = 0.2, a = 0.5, b = 0.1,
    d1 = 1, d2 = 1, sided = "one")
}
\arguments{
  \item{h}{sampling interval. It can be a numeric vector or
  left undefined. See 'Details'}

  \item{H}{decision interval. It can be a numeric vector or
  left undefined. See 'Details'}

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

  \item{nlevels}{30. It works only when \code{contour.plot}
  is TRUE.}

  \item{sided}{distinguish between one-, two-sided and
  Crosier's modified two-sided CUSUM scheme by choosing
  "one", "two", and "Crosier", respectively. See details in
  \code{\link[spc]{xcusum.arl}}}

  \item{par}{initial values for the parameters to be
  optimized over. It can be a vector of length 2 or 3. See
  'Details'}

  \item{contour.plot}{a logical value indicating wether a
  contour plot should be drawn. Default is FALSE.}

  \item{call.print}{a logical value indicating whether the
  "call" should be printed on the contour plot. Default is
  TRUE}

  \item{...}{other arguments to be passed to \code{optim}
  function.}
}
\value{
  The \code{ecoCusum} function returns an object of class
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
  samples. The \code{echCusum} function returns the
  calculated ECH value only.
}
\description{
  Calculate the optimum parameters of n(sample size),
  h(sampling interval), k(reference value) and H(decision
  interval) for Economic Design of the CUSUM control chart.
  For more information about the reference value see
  'Details'.
}
\details{
  When parameter \code{par} is specified, optimization
  algorithms would be used as default. \code{par} can be
  specified as: \code{par = c(h, H)} where \code{h} and
  \code{H} are the initial values of smapling interval and
  decision interval when \code{n} is specified; or
  \code{par = c(h, H, n)}. Good inital values may lead to
  good optimum results.

  When parameters \code{h}, \code{H}, \code{n} are all
  undefined, \code{ecoCusum} function will try to find the
  global optimum point to minimize the ECH (Expected Cost
  per Hour) using optimization algorithms (\code{optim}
  function), but in this case \code{n} would not be
  integer. It is usually helpful for the experimenter to
  find the region where the optimum point may exist
  quickly. When \code{h} and \code{H} are undefined but
  \code{n} is given as an integer vector, \code{ecoCusum}
  function will try to find the optimum point for each
  \code{n} value using optimization algorithms. When
  \code{h}, \code{H} and \code{n} are all given,
  \code{ecoCusum} function will use a "grid method" way to
  calculate the optimum point, that is ECH for all the
  combinations of the parameters will be calculated. The
  "grid method" way is much slower than using optimization
  algorithms, but it would be a good choice when
  optimization algorithms fail to work well.

  There is strong numerical and theoretical evidence that
  for given L1, the value of L0 approaches its maximum when
  k(reference value) is chosen mid-way the between AQL and
  the RQL: $k = mu0 + 0.5*delta*sigma (Appl. Statist.(1974)
  23, No. 3, p. 420). For this reason we treat k as a
  constant value and optimize n, h and H. For cost
  parameters either {P0, P1} or {C0, C1} is needed.  If P0
  and P1 are given, they will be used first, else C0 and C1
  will be used. For economic design of the CUSUM chart,
  when \code{d1} and \code{d2} are both 1, only if the
  difference between P0 and P1 keeps the same, the results
  are identical. If the difference between C0 and C1 keeps
  the same, the optimum parameters are almost the same but
  the ECH(Expected Cost per Hour) values will change.

  \code{echCusum} is used to calculate the ECH (Expected
  Cost per Hour) for one given design point.
}
\examples{
#Chiu, W.K. (1974). Applied Statistics, 23, p427 Table3, row 1-2,14
## LINE 1
## global optimization to h, H and n, when lambda = 0.01, "Nelder-Mead" optimization algorithm doesn't work
#(y <- ecoCusum( P0=150,P1=50,Cr=30,d1=0,d2=0))
## we can try other algorithms:
(y1 <- ecoCusum( P0=150,P1=50,Cr=30,d1=0,d2=0,method="BFGS"))
# Based on the global optimum above, we specify the range of the
# parameters like this
(yy1 <- ecoCusum( h=seq(1.3,1.45,by=.01), H=seq(.5,0.6,by=.01),n=4:6,
P0=150,P1=50,Cr=30,d1=0,d2=0))
## LINE 2
(y2 <- ecoCusum( P0=150,P1=50,Cr=30,d1=0,d2=0,lambda=0.05))
(yy2 <- ecoCusum( h=seq(.6,0.7,by=.01), H=seq(.5,0.6,by=.01),n=3:6,
P0=150,P1=50,Cr=30,d1=0,d2=0,lambda=0.05))
contour(yy2)
## LINE 14
(y14 <- ecoCusum(n=30,P0=150,P1=50,Cr=30,delta=0.5,d1=0,d2=0,method="L-BFGS-B"))
(yy14 <- ecoCusum(h=seq(2.55,2.65,by=0.01),H=seq(0.3,0.4,by=0.01),
n=28:30,P0=150,P1=50,Cr=30,delta=0.5,d1=0,d2=0))
#Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
ecoCusum(lambda=.05,P0=110,P1=10,Cr=25,Cf=50,Tr=0,Tf=0,Tc=1,T0=.0167,a=1)
ecoCusum(h=seq(0.75,0.85,by=.01),H=seq(.55,0.65,by=.01),n=4:6,lambda=.05,
P0=110,P1=10,Cr=25,Cf=50,Tr=0,Tf=0,Tc=1,T0=.0167,a=1)
}
\references{
  Weicheng Zhu, Changsoon Park (2013), {edcc: An R Package
  for the Economic Design of the Control Chart}.
  \emph{Journal of Statistical Software}, 52(9), 1-24.
  \url{http://www.jstatsoft.org/v52/i09/}

  Lorenzen and Vance (1986). The economic design of control
  charts: a unified approach, \emph{Technometrics}, 28.
  3-10.

  Chiu, W.K. (1974). The economic design of CUSUM charts
  for controlling normal means, \emph{Journal of the Royal
  Statistical Society. Series C (Applied Statistics)},
  23(3), 420-433.
}
\seealso{
  \code{\link{ecoXbar}}, \code{\link{ecoEwma}},
  \code{\link[spc]{xcusum.arl}},
  \code{\link[stats]{optim}}, \code{\link{update.edcc}},
  \code{\link{contour}}
}

