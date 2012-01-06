\name{update.edcc}
\alias{update.edcc}
\title{Update for a "edcc" class object}
\usage{
  \method{update}{edcc} (object, ..., evaluate = TRUE)
}
\arguments{
  \item{object}{an object of "edcc" class}

  \item{...}{additional arguments to the call, or arguments
  with anged values.}

  \item{evaluate}{If true evaluate the new call else return
  the call.}
}
\value{
  the fitted object
}
\description{
  'update' will update and (by default) re-fit a model.  It
  does this by extracting the call stored in the object,
  updating the call and (by default) evaluating that call.
}
\details{
  S3 method for update.
}
\examples{
x <- ecoXbar(P0=110,P1=10)
update(x,P0=NULL,P1=NULL,C0=10,C1=110)
}

