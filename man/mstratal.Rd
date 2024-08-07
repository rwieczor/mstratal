% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mstratal.R
\name{mstratal}
\alias{mstratal}
\title{Multivariate optimal stratification and allocation}
\usage{
mstratal(
  data,
  L,
  cc,
  takeall = TRUE,
  method = "rnabox",
  opt_alg = "subplex",
  p_min = 0.1,
  p_max = 0.9,
  maxit1 = 10,
  maxit2 = 100,
  rel_tol = 0.01,
  min_size = 2,
  verbose = TRUE,
  history = FALSE,
  ratios = NULL
)
}
\arguments{
\item{data}{- matrix or data.frame with stratified numerical variables (in columns)}

\item{L}{- number of strata}

\item{cc}{- vector with given coefficients of variation for stratified variables}

\item{takeall}{- if TRUE take-all stratum is used (default)}

\item{method}{- string parameter, choice of algorithm for optimal allocation with given
strata and box constraints (generalization of Neyman allocation),
default method="rnabox", which uses algorithm "RNABOX"
from the paper Wesolowski, Wieczorkowski, Wojciak (2024),
other parameter could be method="capacity" i.e. algoritm "Capacity Scaling"
from the paper Friedrich, Münnich, Vries, Wagner (2015).}

\item{opt_alg}{- string parameter, choice of algorithm for numerical optimization
can be: simplex for Nelder-Mead simplex, or 'subplpex' for subplex algorithm}

\item{p_min}{- minimal value for search of optimal value of
power in cumulative power density rule (default 0.1)}

\item{p_max}{- maximal value for search of optimal value of
power in cumulative power density rule (default 0.9)}

\item{maxit1}{- maximal number of iterations at first step of stratification,
where optimal value of power in cumulative power density rule is found
(default 10)}

\item{maxit2}{- maximal number of iterations in second step of stratification, where
Nelder-Mead or subplex algorithm for minimization of objective function
(default 100)}

\item{rel_tol}{- relative tolerance, used as stopping rule in using sequentially
selected optimization alorithm (default 0.01)}

\item{min_size}{- minimal sample size in strata (default 2), if min_size < 1 then
minimal sample fraction in strata}

\item{verbose}{- if TRUE then diagnostic output is printed}

\item{history}{- if TRUE then output contains list of sample sizes
from consecutive generations of algorithm}

\item{ratios}{- estimates for ratios, used when stratified variables
are given after linearization}
}
\value{
list with elements: bh - data frame with columns of strata boundaries for stratification variables (bh1,bh2,...),
  and nh - corresponding sample allocation in obtained strata;
  if parameter 'history' is set to TRUE then output have additional
  element - vector 'n_history' with sample sizes obtained in the process
  of optimization.
}
\description{
Function for optimal stratification and allocation of multivariate population.
New numerical algorithm is based on ideas from the papers:
Lednicki, Wieczorkowski (2003),  Friedrich et al. (2015), and
Wesolowski, Wieczorkowski, Wojciak (2024)
}
\examples{
# Example of solving bi-variate stratification and allocation problem

# Generation of correlated lognormal variables (x,y)
set.seed(3456)
ro <- 0.5 # correlation coefficient for lognormal vaiables
#
# correlation coefficient for normal vaiables (theoretical formula)
(ro_norm <- log(.5 * (exp(1) - 1) + 1))
#
x <- rnorm(10000)
z <- rnorm(10000)
y <- x * ro_norm + z * sqrt(1 - ro_norm^2)
x <- exp(x)
y <- exp(y)
cor(x, y)
#
L <- 5
c <- 0.01

ex <- mstratal(cbind(x, y), L, c(c, c),
  opt_alg = "simplex",
  maxit1 = 20, maxit2 = 100, rel_tol = 0.01,
  verbose = TRUE,
  history = TRUE
)
ex
sum(ex$nh) # total sample size
# Plot for optimization history
n_history <- ex$n_history
plot(n_history,
  cex = 0.5,
  ylim = c(min(n_history) - 10, max(n_history))
)
lines(n_history)
abline(h = min(n_history), col = 3)

}
\references{
- Lednicki B., Wieczorkowski R., Optimal stratification and sample
allocation between subpopulations and strata, Statistics in Transition (2003),
Vol. 6, No. 2,  287-305.

- Ulf Friedrich, Ralf Münnich, Sven de Vries, Matthias Wagner,
Fast integer-valued algorithms for optimal allocations under constraints in stratified sampling,
Computational Statistics and Data Analysis 92 (2015), 1-12

- Wesolowski J., Wieczorkowski R., Wojciak W. (2024),
Recursive Neyman Algorithm for Optimum Sample Allocation under Box Constraints on
Sample Sizes in Strata, (to be published in *Survey Methodology*),
arxiv version available at: \url{https://arxiv.org/abs/2304.07034}
}
