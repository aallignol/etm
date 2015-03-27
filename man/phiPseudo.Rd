\name{phiPseudo}
\alias{phiPseudo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Pseudo Value Regression for the Expected Excess Length of Stay
}
\description{
  Pseudo value regression for the expected excess length of stay for
  each landmark time
}
\usage{
phiPseudo(data, state.names, tra, cens.name, s = 0, formula, timepoints, ncores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data.frame in a format suitable for \code{\link{etm}}.}
  \item{state.names}{A vector of characters giving the states names.}
  \item{tra}{A quadratic matrix of logical values describing the possible
    transitions within the multistate model.}
  \item{cens.name}{ A character giving the code for censored
    observations in the column 'to' of \code{data}. If there is no
    censored observations in your data, put \code{NULL}}
  \item{s}{Starting time for computing the transition probabilities.}
  \item{formula}{A formula with the covariates at the right of a
    \code{~} operator. The left hand side can be left empty.}
  \item{timepoints}{Landmark time points at which the pseudo values are
    computed.}
  \item{ncores}{Number of cores used if doing parallel computation using
  the \pkg{parallel} package}
}

\details{
  The function calculates the pseudo-observations for the extra
  length-of-stay at several landmark time points for each
  individual. These pseudo-observations can then
  be used to fit a direct regression model using generalized estimating
  equation (e.g., package \pkg{geepack}).

  Computation of the pseudo-observations can be parallelised using the
  \code{mclapply} function of the \pkg{parallel} package. See argument
  \code{ncores}.
}
\value{
  An object of class \code{phiPseudo} with the following components:
  \item{pseudoData}{a data.frame containing \code{id}, computed pseudo
    values (see details) and the covariates as specified in the formula}
  \item{phi}{Estimates of excess LoS in the whole sample}
  \item{ps}{}
}
\references{
  Andersen, P.K, Klein, J.P, Rosthoj, S. (2003). Generalised
  linear  models for correlated pseudo-observations,
  with applications to multi-state models. \emph{Biometrika},
  90(1):15--27.
}

\author{
  Arthur Allignol \email{arthur.allignol@uni-ulm.de}
}

\seealso{
  \code{\link[parallel]{mclapply}}, \code{\link[etm]{clos}}
}

\examples{
## TODO
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
