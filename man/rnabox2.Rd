% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mstratal.R
\name{rnabox2}
\alias{rnabox2}
\title{Optimal univariate allocation under lower and upper constraints
  for stratified sampling, version with constraint imposed on the variance}
\usage{
rnabox2(V0, Nh, Sh, mh = NULL, Mh = NULL)
}
\arguments{
\item{V0}{- upper limit for value of variance which must be attained for computed optimal allocation}

\item{Nh}{- population sizes in strata}

\item{Sh}{- standard deviations for given variable in strata}

\item{mh}{- lower constraints for sample sizes in strata}

\item{Mh}{- upper constraints for sample sizes in strata}
}
\value{
vector of optimal allocation sizes, which can be rounded to
    integers using functions: stratallo::ran_round or stratallo::round_oric
}
\description{
Modified recursive Neyman algorithm for optimal allocation in stratified sampling
  with lower and upper constraints
  (used RNABOX algorithm from 'stratallo' package)
}
\examples{

N <- c(454, 10, 116, 2500, 2240, 260, 39, 3000, 2500, 400)
S <- c(0.9, 5000, 32, 0.1, 3, 5, 300, 13, 20, 7)
A <- N * S
m <- c(322, 3, 57, 207, 715, 121, 9, 1246, 1095, 294) # lower bounds
M <- N # upper bounds

n <- 6000
# Allocation using RNABOX algorithm with required sample size
nh <- stratallo::rnabox(n, A, M, m)
var0 <- stratallo::var_st_tsi(nh, N, S) # variance of stratified estimator

nh2 <- rnabox2(var0, N, S, m, M)
print(nh)
print(nh2)

#' @export

}
\references{
- Wesolowski J., Wieczorkowski R., Wojciak W. (2024),
Recursive Neyman Algorithm for Optimum Sample Allocation under Box Constraints on
Sample Sizes in Strata,(to be published in *Survey Methodology*),
arxiv version available at: \url{https://arxiv.org/abs/2304.07034}
}
