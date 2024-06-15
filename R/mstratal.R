#' @title Optimal univariate allocation under lower and upper constraints
#'   for stratified sampling, version with constraint imposed on the variance
#' @description Modified recursive Neyman algorithm for optimal allocation in stratified sampling
#'   with lower and upper constraints
#'   (used RNABOX algorithm from 'stratallo' package)
#'
#' @param V0 - upper limit for value of variance which must be attained for computed optimal allocation
#' @param Nh - population sizes in strata
#' @param Sh - standard deviations for given variable in strata
#' @param mh - lower constraints for sample sizes in strata
#' @param Mh - upper constraints for sample sizes in strata
#'
#' @return  vector of optimal allocation sizes, which can be rounded to
#'     integers using functions: stratallo:ran_round or stratallo::round_oric
#' @export
#'
#' @references
#' - Wesolowski J., Wieczorkowski R., Wojciak W. (2024),
#' Recursive Neyman Algorithm for Optimum Sample Allocation under Box Constraints on
#' Sample Sizes in Strata,(to be published in *Survey Methodology*),
#' arxiv version available at: \url{https://arxiv.org/abs/2304.07034}
#'
#' @examples
#'
#' N <- c(454, 10, 116, 2500, 2240, 260, 39, 3000, 2500, 400)
#' S <- c(0.9, 5000, 32, 0.1, 3, 5, 300, 13, 20, 7)
#' A <- N * S
#' m <- c(322, 3, 57, 207, 715, 121, 9, 1246, 1095, 294) # lower bounds
#' M <- N # upper bounds
#'
#' n <- 6000
#' # Allocation using RNABOX algorithm with required sample size
#' nh <- stratallo::rnabox(n, A, M, m)
#' var0 <- stratallo::var_st_tsi(nh, N, S) # variance of stratified estimator
#'
# Allocation using version of algorithm with constraint on variance
#' nh2 <- rnabox2(var0, N, S, m, M)
#' print(nh)
#' print(nh2)
#'
#' #' @export
#'
rnabox2 <- function(V0, Nh, Sh, mh = NULL, Mh = NULL) {
  Sh[Sh == 0] <- 1e-8
  Ah <- Sh * Nh
  Ah2 <- Ah^2
  nh <- stratallo::rnabox(
    n = V0 + sum(Ah * Sh), A = Ah,
    bounds1 = Ah2 / mh,
    bounds2 = Ah2 / Mh,
    check_violations1 = .Primitive(">="), # RRNA variant
    check_violations2 = .Primitive("<=")
  )
  nh <- Ah2 / nh

  return(nh)
}




#' @title Cumulative power density rule for initial stratification
#' @description Function for computation initial stratification by generalized
#' cumroot density rule i.e. cumulative power density rule
#'
#' @param x - stratification numerical variable
#' @param L - number of strata
#' @param p - value of power,  from interval (0,1)
#'
#' @return initial strata boundaries
#'
#' @examples
#' # Examples for cumfp
#' # generation of artificial data
#' x1 <- rchisq(10000, 1)
#'
#' cumfp(x = x1, L = 5, p = 0.5) # classical cumroot rule
#'
#' cumfp(x = x1, L = 5, p = 0.3)
#'
#' @export
#'
cumfp <- function(x, L, p) {
  hst <- graphics::hist(x,
    plot = FALSE,
    nclass = grDevices::nclass.FD
  )

  mids <- hst$mids
  fy <- hst$dens

  cf <- cumsum((fy)^(p))

  min <- min(cf)
  max <- max(cf)
  delta <- (max - min) / L
  # apr<-approx(mids,cf)

  gr <- double(L - 1)
  xi <- double(L - 1)
  for (i in 1:(L - 1))
  {
    gr[i] <- min + i * delta
    xi <- min(which(gr[i] <= cf))
    # xi<-min(which(gr[i]<=apr$y))
    gr[i] <- mids[xi]
    # gr[i]<-(apr$x)[xi]
    # cat("xi gri",xi," ",gr[i],"\n")
  }

  return(gr)
}




#' @title Optimal sample allocation for given strata
#' @description Function for optimal sample allocation in L strata defined by boundaries
#'
#' @param gr - strata boundaries, as vector concatenated with (L-1) values for each variable
#' @param xx - matrix with considered variables (in columns)
#' @param L - number of strata
#' @param cc - vector of given coefficients of variation for stratification variables
#' @param cv - if TRUE then coefficients of variations for variables for obtained allocations are printed
#' @param method - algorithm used for optimal allocation: 'rnabox' (default) or
#'    'capacity'; method 'rnabox' from *stratallo* package is based on generalization of
#'    Neyman optimal allocation formula taking into account lower an upper bounds for
#'    sample sizes in strata; method 'capacity' uses integer optimal allocation
#'    'CapacityScaling' algorithm from Friedrich et al. (2015) paper
#' @param min_size - minimal sample size in strata (default 2); if min_size < 1 then minimal
#'     sample fraction in strata
#'
#' @return numerical vector with optimal sample allocation in strata
#' defined by input boundaries;
#'
#' @references
#' - Ulf Friedrich, Ralf Münnich, Sven de Vries, Matthias Wagner,
#' Fast integer-valued algorithms for optimal allocations under constraints in stratified sampling,
#' Computational Statistics and Data Analysis 92 (2015), 1-12
#'
#' - Wesolowski J., Wieczorkowski R., Wojciak W. (2024),
#' Recursive Neyman Algorithm for Optimum Sample Allocation under Box Constraints on
#' Sample Sizes in Strata, (to be published in *Survey Methodology*),
#' arxiv version available at: \url{https://arxiv.org/abs/2304.07034}
#'
#' @examples
#'
#' # Generation of correlated lognormal variables (x,y)
#' set.seed(3456)
#' ro <- 0.5 # correlation coefficient for lognormal variables
#' #
#' # correlation coefficient for normal variables (theoretical formula)
#' (ro_norm <- log(.5 * (exp(1) - 1) + 1))
#' #
#' x <- rnorm(10000)
#' z <- rnorm(10000)
#' y <- x * ro_norm + z * sqrt(1 - ro_norm^2)
#' x <- exp(x)
#' y <- exp(y)
#' cor(x, y)
#' #
#' dataxy <- cbind(x, y) # data matrix
#' L <- 3 # number of strata
#' cv <- 0.05 # coefficient of variation for estimators
#' boundaries <- c(cumfp(dataxy[, 1], L, 0.5), cumfp(dataxy[, 2], L, 0.5))
#' nh <- al_nh(boundaries, # strata boundaries
#'   xx = dataxy,
#'   L = L,
#'   cc = c(cv, cv),
#'   cv = TRUE,
#'   method = "rnabox"
#' )
#' cat("Optimal sample allocation for given strata = ", nh, "\n")
#' cat("Optimal sample size = ", sum(nh), "\n")
#' #
#' nh <- al_nh(boundaries,
#'   xx = dataxy,
#'   L = L,
#'   cc = c(cv, cv),
#'   cv = TRUE,
#'   method = "rnabox",
#'   min_size = 20 # change for minimal sample size
#' )
#' cat("Optimal sample allocation for given strata = ", nh, "\n")
#' cat("Optimal sample size = ", sum(nh), "\n")
#'
#' @export
#'
al_nh <- function(gr, xx, L, cc, cv = FALSE, method = "rnabox", min_size = 2) {
  # count_alok <<- count_alok + 1 # global counter of function calls

  big_nh <- 1e10

  jitter4gr <- function(x) {
    # jittering for unique vector 'x'
    y <- c(NA, diff(x))
    iy <- which(y == 0)
    iy <- union(iy - 1, iy)
    jx <- jitter(x)
    jx[-iy] <- x[-iy]
    return(sort(jx))
  }


  ndim <- ncol(xx)
  N <- nrow(xx)

  h <- rep(1, N)
  for (i in 1:ndim)
  {
    gri <- gr[(1 + (i - 1) * (L - 1)):(i * (L - 1))]
    grix <- c(min(xx[, i]), gri, max(xx[, i]))
    # while (!identical(grix, unique(grix)) && (max(gri) < max(xx[, i]))) {
    #   # gri <- runif(length(gri),min = min(xx[,i])+0.01 , max = max(xx[,i]*1.1))
    #   gri <- jitter4gr(gri)
    #   grix <- c(min(xx[, i]), gri, max(xx[, i]))
    # }

    if (length(unique(grix)) < L + 1) {
      gri <- seq(min(xx[, i]), max(xx[, i]), length.out = L + 1)
      gri <- sort(gri)[2:(L - 1)]
    }
    h <- pmax(h, cut(xx[, i], c(min(xx[, i]), sort(gri), max(xx[, i])), include.lowest = TRUE, labels = FALSE))
  }

  # Xbar<-apply(xx,2,mean)
  Xbar <- matrixStats::colMeans2(xx)

  Nh <- table(h)
  ## if (length(Nh) != L) cat("Too many strata L !","\n")
  if (any(is.na(Nh)) | any(is.nan(Nh)) | length(Nh) != L) {
    return(rep(1e15, L))
  }

  # Wh<-Nh/N
  S2h <- stats::aggregate(xx, by = list(h = h), FUN = function(x) {
    matrixStats::colVars(matrix(x))
  })
  names(S2h) <- c("h", paste0("S2", 1:ndim))

  Xh <- stats::aggregate(xx, by = list(h = h), FUN = function(x) {
    matrixStats::colMeans2(matrix(x))
  })
  names(Xh) <- c("h", paste0("X", 1:ndim))


  nh <- matrix(0, L, ndim)

  lo.str <- rep(min_size, L)
  if (min_size < 1) {
    # condition for minimal sample fraction in strata
    lo.str <- round(min_size * Nh)
  }
  lo.str <- pmin(pmax(2, lo.str), Nh)
  lo.str[L] <- Nh[L]

  if (any(is.na(S2h))) {
    return(rep(big_nh, L))
  }

  for (i in 1:ndim) {
    # V0<-(cc[i]^2)*(Xbar[i]*N*Xbar[i]*N)
    V0 <- (cc[i] * Xbar[i] * N)^2

    # cat("i: V0= ",i,":",V0," Nh= ",Nh,"  S2 = ",S2h[,i+1]," lo.str= ",lo.str,"\n")
    if (method == "capacity") {
      nh[, i] <- stratallo:::CapacityScaling2(V0, Nh, sqrt(S2h[, i + 1]), lo.str, Nh)
    } else if (method == "rnabox") nh[, i] <- stratallo::round_oric(rnabox2(V0, Nh, sqrt(S2h[, i + 1]), lo.str, Nh))
  }
  nh <- apply(nh, 1, max)
  n <- sum(nh)

  if (any(is.na(nh)) | any(is.nan(nh))) {
    return(rep(big_nh, L))
  }

  if (any(nh > Nh)) cat("nh > Nh !!! ", "\n")

  if (nh[L] < Nh[L]) cat("nh[L] < Nh[L] !!! ", "\n")

  if (cv == TRUE) {
    cvout <- double(ndim)
    for (i in 1:ndim)
    {
      cvout[i] <- (sqrt(sum((Nh / nh) * (Nh - nh) * S2h[, i + 1])) / (N * Xbar[i]))
    }

    cat("\nCVs for stratification variables = ", round(cvout, 4), "\n")
  }

  return(nh)
}




#' @title Multivariate optimal stratification and allocation
#' @description Function for optimal stratification and allocation of multivariate population.
#' New numerical algorithm is based on ideas from the papers:
#' Lednicki, Wieczorkowski (2003),  Friedrich et al. (2015), and
#' Wesolowski, Wieczorkowski, Wojciak (2024)
#'
#' @param xx - matrix or data.frame with stratified numerical variables (in columns)
#' @param L - number of strata
#' @param cc - vector with given coefficients of variation for stratified variables
#'
#' @param method - string parameter, choice of algorithm for optimal allocation with given
#'    strata and box constraints (generalization of Neyman allocation),
#'    default metod="rnabox", which uses algorithm "RNABOX"
#'    from the paper Wesolowski, Wieczorkowski, Wojciak (2024),
#'    other parameter could be metod="capacity" i.e. algoritm "Capacity Scaling"
#'    from the paper Friedrich, Münnich, Vries, Wagner (2015).
#'
#' @param opt_alg - string parameter, choice of algorithm for numerical optimization
#'   can be: simplex for Nelder-Mead simplex, or 'sublpex' for subplex algorithm
#' @param p_min - minimal value for search of optimal value of
#'         power in cumulative power density rule (default 0.1)
#' @param p_max - maximal value for search of optimal value of
#'         power in cumulative power density rule (default 0.9)
#' @param maxit1 - maximal number of iterations at first step of stratification,
#'    where optimal value of power in cumulative power density rule is found
#'    (default 10)
#' @param maxit2 - maximal number of iterations in second step of stratification, where
#'    Nelder-Mead or subplex algorithm for minimization of objective function
#'    (default 100)
#'
#' @param rel_tol - relative tolerance, used as stopping rule in using sequentially
#'       selected optimization alorithm (default 0.01)
#' @param min_size - minimal sample size in strata (default 2), if min_size < 1 then
#'    minimal sample fraction in strata
#' @param verbose - if TRUE then diagnostic output is printed
#' @param history - if TRUE then output contains list of sample sizes
#' from consecutive generations of algorithm
#'
#' @return list with elements: bh - data frame with columns of strata boundaries for stratification variables (bh1,bh2,...),
#'   and nh - corresponding sample allocation in obtained strata;
#'   if parameter 'history' is set to TRUE then output have additional
#'   element - vector 'n_history' with sample sizes obtained in the process
#'   of optimization.
#'
#' @references
#' - Lednicki B., Wieczorkowski R., Optimal stratification and sample
#' allocation between subpopulations and strata, Statistics in Transition (2003),
#' Vol. 6, No. 2,  287-305.
#'
#' - Ulf Friedrich, Ralf Münnich, Sven de Vries, Matthias Wagner,
#' Fast integer-valued algorithms for optimal allocations under constraints in stratified sampling,
#' Computational Statistics and Data Analysis 92 (2015), 1-12
#'
#' - Wesolowski J., Wieczorkowski R., Wojciak W. (2024),
#' Recursive Neyman Algorithm for Optimum Sample Allocation under Box Constraints on
#' Sample Sizes in Strata, (to be published in *Survey Methodology*),
#' arxiv version available at: \url{https://arxiv.org/abs/2304.07034}
#'
#' @examples
#' # Example of solving bi-variate stratification and allocation problem
#'
#' # Generation of correlated lognormal variables (x,y)
#' set.seed(3456)
#' ro <- 0.5 # correlation coefficient for lognormal vaiables
#' #
#' # correlation coefficient for normal vaiables (theoretical formula)
#' (ro_norm <- log(.5 * (exp(1) - 1) + 1))
#' #
#' x <- rnorm(10000)
#' z <- rnorm(10000)
#' y <- x * ro_norm + z * sqrt(1 - ro_norm^2)
#' x <- exp(x)
#' y <- exp(y)
#' cor(x, y)
#' #
#' L <- 5
#' c <- 0.01
#'
#' ex <- mstratal(cbind(x, y), L, c(c, c),
#'   opt_alg = "simplex",
#'   maxit1 = 20, maxit2 = 100, rel_tol = 0.01,
#'   verbose = TRUE,
#'   history = TRUE
#' )
#' ex
#' sum(ex$nh) # total sample size
#' # Plot for optimization history
#' n_history <- ex$n_history
#' plot(n_history,
#'   cex = 0.5,
#'   ylim = c(min(n_history) - 10, max(n_history))
#' )
#' lines(n_history)
#' abline(h = min(n_history), col = 3)
#'
#' @export
#'

mstratal <- function(xx, L, cc,
                     method = "rnabox",
                     opt_alg = "subplex",
                     p_min = 0.1, p_max = 0.9,
                     maxit1 = 10,
                     maxit2 = 100,
                     rel_tol = 0.01,
                     min_size = 2,
                     verbose = TRUE,
                     history = FALSE) {

  xx <- apply(xx,2,as.double)

  ndim <- ncol(xx)

  lower <- NULL
  for (i in 1:ndim) lower[(1 + (i - 1) * (L - 1)):(i * (L - 1))] <- min(xx[, i])
  upper <- NULL
  for (i in 1:ndim) upper[(1 + (i - 1) * (L - 1)):(i * (L - 1))] <- max(xx[, i])

  # method used in optimization for cumulative power rule
  # if (ndim == 1) method0 <- "Brent" else method0 <- "Nelder-Mead"

  p0 <- rep(0.5, ndim)
  if (verbose) cat("Initial powers in cumulative rule = ", p0, "\n")

  # gr0<-NULL
  # for (i in 1:ndim) gr0<-c(gr0,cumfp(xx[,i],L,p0[i]))
  gr0 <- double(ndim * (L - 1))
  for (i in 1:ndim) gr0[(1 + (i - 1) * (L - 1)):(i * (L - 1))] <- cumfp(xx[, i], L, p0[i])

  nh0 <- al_nh(gr0, xx, L, cc, method = method, min_size = min_size)
  if (verbose) cat("Initial sample size   = ", sum(nh0), "\n")

  if (history) n_history <- sum(nh0)

  if (opt_alg == "simplex") {
    # p0 <- optim(p0,
    p0 <- nloptr::neldermead(p0,
      function(r) {
        # gr<-NULL ;for (i in 1:ndim) gr<-c(gr,cumfp(xx[,i],L,r[i]))
        gr <- double(ndim * (L - 1))
        for (i in 1:ndim) gr[(1 + (i - 1) * (L - 1)):(i * (L - 1))] <- cumfp(xx[, i], L, r[i])
        sum(al_nh(gr, xx, L, cc, method = method, min_size = min_size))
      },
      # method = method0,
      lower = rep(p_min, ndim), upper = rep(p_max, ndim),
      # control = list(maxit = maxit1)
      control = list(maxeval = maxit1)
    )$par
  } else if (opt_alg == "subplex") {
    p0 <- nloptr::sbplx(p0,
      function(r) {
        # gr<-NULL ;for (i in 1:ndim) gr<-c(gr,cumfp(xx[,i],L,r[i]))
        gr <- double(ndim * (L - 1))
        for (i in 1:ndim) gr[(1 + (i - 1) * (L - 1)):(i * (L - 1))] <- cumfp(xx[, i], L, r[i])
        sum(al_nh(gr, xx, L, cc, method = method, min_size = min_size))
      },
      lower = rep(p_min, ndim), upper = rep(p_max, ndim),
      control = list(maxeval = maxit1)
    )$par
  } else {
    stop("Bad parameter 'opt_alg': should be 'simplex' or 'subplex'")
  }

  if (verbose) cat("Optimal powers in cumulative rule = ", p0, "\n")

  # gr0<-NULL ; for (i in 1:ndim) gr0<-c(gr0,cumfp(xx[,i],L,p0[i]))
  for (i in 1:ndim) gr0[(1 + (i - 1) * (L - 1)):(i * (L - 1))] <- cumfp(xx[, i], L, p0[i])

  nh0 <- al_nh(gr0, xx, L, cc, method = method, min_size = min_size)
  if (verbose) cat("Sample size for optimal cumulative power rule   = ", sum(nh0), "\n")

  sumpop <- sum(nh0)
  if (history) n_history <- c(n_history, sumpop)

  # gropt<-gr0; nhopt<-nh.nd1(gropt,xx,L,cc,cv=TRUE)

  if (opt_alg == "simplex") {
    while (1) {
      # gropt <- optim(gr0, function(z) {
      gropt <- nloptr::neldermead(gr0, function(z) {
        sum(al_nh(z, xx, L, cc, method = method, min_size = min_size))
      },
      lower = lower, upper = upper,
      # gr = function(z) {pracma::grad(function(z) {sum( al_nh(z,xx,L,cc, method=method, min_size=min_size) )},z)},
      # method = "Nelder-Mead", control = list(maxit = maxit2),
      control = list(maxeval = maxit2)
      )$par

      nhopt <- al_nh(gropt, xx, L, cc, cv = TRUE, method = method, min_size = min_size)
      # break
      if (verbose) cat("Sample size from sequential simplex optimization = ", sum(nhopt), "\n")
      if (history) n_history <- c(n_history, sum(nhopt))

      if (abs(sum(nhopt) - sumpop) / sumpop <= rel_tol) {
        break
      }
      gr0 <- gropt
      sumpop <- sum(nhopt)
    }
    nhopt <- al_nh(gropt, xx, L, cc, cv = TRUE, method = method, min_size = min_size)
    if (verbose) cat("Final sample size from sequential simplex optimization = ", sum(nhopt), "\n")
    if (history) n_history <- c(n_history, sum(nhopt))
  } else if (opt_alg == "subplex") {
    while (1) {
      gropt <- nloptr::sbplx(gr0, function(z) {
        sum(al_nh(z, xx, L, cc, method = method, min_size = min_size))
      },
      lower = lower, upper = upper,
      control = list(maxeval = maxit2)
      # gropt <- subplex::subplex(gr0, function(z) {
      #    sum(al_nh(z, xx, L, cc, method = method, min_size = min_size))
      #  },
      # hessian = TRUE,
      # control = list(maxit = maxit2)
      )$par

      # gropt<-sort(gropt)

      nhopt <- al_nh(gropt, xx, L, cc,
        cv = TRUE, method = method,
        min_size = min_size
      )
      # break
      if (verbose) cat("Sample size from sequential subplex optimization = ", sum(nhopt), "\n")
      if (history) n_history <- c(n_history, sum(nhopt))
      if (abs(sum(nhopt) - sumpop) / sumpop <= rel_tol) {
        break
      }
      gr0 <- gropt
      sumpop <- sum(nhopt)
    }
    nhopt <- al_nh(gropt, xx, L, cc, cv = TRUE, method = method, min_size = min_size)
    if (verbose) cat("Final sample size from sequential subplex optimization = ", sum(nhopt), "\n")
    if (history) n_history <- c(n_history, sum(nhopt))
  } else {
    stop("Bad parameter 'opt_alg': should be 'simplex' or 'subplex'")
  }

  # info for output: strata boundaries and sample sizes in strata
  bh <- matrix(0, L, ndim)

  for (i in 1:ndim) bh[1:(L - 1), i] <- sort(gropt[(1 + (i - 1) * (L - 1)):(i * (L - 1))])
  bh[L, ] <- apply(xx, 2, max)

  bh <- data.frame(bh)
  names(bh) <- paste0("bh", 1:ndim)
  # bh <- cbind(bh, nh = nhopt)

  if (history) {
    return(list(bh = bh, nh = nhopt, n_history = n_history))
  } else {
    return(list(bh = bh, nh = nhopt))
  }
}
