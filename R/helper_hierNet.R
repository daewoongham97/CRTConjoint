# slight computational modifications of functions from hierNet package to quicken computation
# All functions here are originally from author Jacob Bien source: https://github.com/cran/hierNet
# All functions here are used with explicit permission from Jacob Bien
predict_new <- function(object, newx, newzz=NULL, ...) {
  n <- nrow(newx)

  newx <- scale(newx, center=object$mx, scale=object$sx)

  p <- ncol(newx)
  cp2 <- p * (p - 1) / 2

  out <- .C("ComputeInteractionsWithIndices",
            as.double(newx),
            as.integer(n),
            as.integer(p),
            z=rep(0, n * cp2),
            i1=as.integer(rep(0, cp2)),
            i2=as.integer(rep(0, cp2)), PACKAGE="CRTConjoint")

  unscaled_z = out$z
  # modifications to speed things more faster
  newzz = demean_center(unscaled_z, n, cp2, object$mzz)$vec
  newx <- as.numeric(newx)
  stopifnot(is.finite(newzz), is.finite(newx))

  yhatt <- Compute.yhat.c(newx, newzz, object) + object$my

  b0 <- object$b0
  yhatt <- b0 + yhatt
  pr <- 1 / (1 + exp(-yhatt))
  return(pr)
}


hierNet_logistic = function(x, y, lam, delta=1e-8, diagonal=FALSE, strong=FALSE, aa=NULL, zz=NULL, center=TRUE,
                            stand.main=TRUE, stand.int=FALSE,
                            rho=nrow(x), niter=100, sym.eps=1e-3,# ADMM params
                            step=2, maxiter=2000, backtrack=0.1, tol=1e-3, # GG descent params
                            trace=0) {
  this.call <- match.call()
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(y %in% c(0,1))
  stopifnot(length(y) == n, lam >= 0, delta >= 0, delta <= 1)
  stopifnot(!is.null(step) && !is.null(maxiter))
  stopifnot(class(lam) == "numeric")
  stopifnot(class(delta) == "numeric")
  stopifnot(class(step) == "numeric", step > 0, maxiter > 0)
  stopifnot(is.finite(x), is.finite(y), is.finite(lam), is.finite(delta))
  lam.l1 <- lam * (1 - delta)
  lam.l2 <- lam * delta
  if (!center)
    cat("WARNING: center=FALSE should almost never be used.  This option is available for special uses only.", fill = TRUE)
  x <- scale(x, center = center, scale = stand.main)
  mx <- attr(x, "scaled:center")
  sx <- attr(x, "scaled:scale")
  if (is.null(aa)) aa <- list(b0=0, bp=rep(0, p), bn=rep(0, p), th=matrix(0, p, p), diagonal=diagonal)
  cp2 <- p * (p - 1) / 2
  out <- .C("ComputeInteractionsWithIndices",
            as.double(x),
            as.integer(n),
            as.integer(p),
            z=rep(0, n * cp2),
            i1=as.integer(rep(0, cp2)),
            i2=as.integer(rep(0, cp2)), PACKAGE="CRTConjoint")
  out_zz = fast_demean(out$z, n, cp2)
  szz = NULL
  mzz = out_zz$means
  zz = out_zz$vec

  xnum <- as.numeric(x)
  out <- ggdescent.logistic(xnum=xnum, zz=zz, y=y, lam.l1=lam.l1, lam.l2=lam.l2, diagonal=diagonal, rho=0, V=matrix(0,p,p),
                            stepsize=step, backtrack=backtrack, maxiter=maxiter,
                            tol=tol, aa=aa, trace=trace)
  out$call <- this.call
  out$lam <- lam
  out$delta <- delta
  out$type <- "logistic"
  out$diagonal <- diagonal
  out$strong <- strong
  out$step <- step
  out$maxiter <- maxiter
  out$backtrack <- backtrack
  out$tol <- tol
  out$mx <- mx
  out$my <- 0
  out$sx <- sx
  out$mzz = mzz
  out$szz = szz
  class(out) <- "hierNet"
  return(out)
}

hierNet_logistic_CV = function(lambda, nfolds = 3, X, y_var, seed = sample(1:1000, size = 1), constraint = TRUE, tol= 1e-3, fold_idx = NULL) {
  if (is.null(fold_idx)) {
    if (!constraint) {
      half = nrow(X)
      random_idx = suppressWarnings(split(sample(half, half, replace = FALSE), as.factor(1:nfolds)))
    } else {
      half = (nrow(X)/2)
      random_idx = suppressWarnings(split(sample(half, half, replace = FALSE), as.factor(1:nfolds)))
      for (i in 1:nfolds) {
        random_idx[[i]] = c(random_idx[[i]], random_idx[[i]] + half)
      }
    }
  } else {
    random_idx = fold_idx
  }


  error_list_prob = list()
  for (i in 1:nfolds) {
    errors_prob = vector()
    test_idx = random_idx[[i]]
    train_idx = (1:nrow(X))[-test_idx]

    for (j in 1:length(lambda)) {
      if( j == 1) {
        invisible(capture.output(cv_hiernets <- hierNet_logistic(X[train_idx, ], y_var[train_idx], lam = lambda[j], tol = tol)))
      } else {
        invisible(capture.output(cv_hiernets <- hierNet_logistic(X[train_idx, ], y_var[train_idx], lam = lambda[j], tol = tol, aa = cv_hiernets)))
      }
      predicted_y_prob = predict_new(cv_hiernets, X[test_idx, ])
      errors_prob[j] = mean((y_var[test_idx] - predicted_y_prob)^2)
    }
    error_list_prob[[i]] = errors_prob
  }


  cv_errors_prob = vector()
  for (i in 1:length(lambda)) {
    cv_errors_prob[i] = mean(sapply(error_list_prob, "[[", i))
  }

  gotten_lam = lambda[which.min(cv_errors_prob)]
  return(gotten_lam)
}

## these functions are taken directly from the hierNet package
ggdescent.logistic <- function(xnum, zz, y, lam.l1, lam.l2, diagonal, rho, V, stepsize, backtrack=0.2, maxiter=100,
                               tol = 1e-3, aa=NULL, trace=1) {
  # See ADMM4 pdf and logistic.pdf for the problem this solves.
  #
  # xnum, zz, y: data (note: zz is a length n*cp2 vector, not a matrix) xnum is x as a (n*p)-vector
  # lam.l1: l1-penalty parameter
  # lam.l2: l2-penalty parameter
  # rho: admm parameter
  # V: see ADMM4 pdf
  # stepsize: step size to start backtracking with
  # backtrack: factor by which step is reduced on each backtrack.
  # maxiter: number of generalized gradient steps to take.
  # tol: stop gg descent if change in objective is below tol.
  # aa: initial estimate of (b0, th, bp, bn)
  # trace: how verbose to be
  #
  #void ggdescent_logistic_R(double *x, int *n, int *p, double *zz, int * diagonal, double *y,
  #			     double *lamL1, double *lamL2, double *rho, double *V, int *maxiter,
  #			     double *curb0, double *curth, double *curbp, double *curbn,
  #			     double *t, int *stepwindow, double *backtrack, double *tol, int *trace,
  #			     double *b0, double *th, double *bp, double *bn) {
  n <- length(y)
  p <- length(xnum) / n
  stopifnot(p == round(p))
  if (diagonal) stopifnot(length(zz) == n * (choose(p,2)+p))
  else stopifnot(length(zz) == n * choose(p,2))
  stepwindow <- 10
  if (is.null(aa)) aa <- list(b0=0, th=matrix(0,p,p), bp=rep(0,p), bn=rep(0,p))
  out <- .C("ggdescent_logistic_R",
            xnum,
            as.integer(n),
            as.integer(p),
            zz,
            as.integer(diagonal),
            as.double(y), # convert from integer to double
            as.double(lam.l1),
            as.double(lam.l2),
            as.double(rho),
            as.double(V),
            as.integer(maxiter),
            as.double(aa$b0),
            as.double(aa$th),
            aa$bp,
            aa$bn,
            as.double(stepsize),
            as.integer(stepwindow),
            as.double(backtrack),
            as.double(tol),
            as.integer(trace),
            b0=as.double(0),
            th=rep(0, p*p),
            bp=rep(0, p),
            bn=rep(0, p),
            PACKAGE="CRTConjoint")
  list(b0=out$b0, bp=out$bp, bn=out$bn, th=matrix(out$th, p, p))
}

Compute.yhat.c <- function(xnum, zz, aa) {
  # aa: list containing bp, bn, th, diagonal
  # note: zz is the n by cp2 matrix, whereas z is the n by p^2 one.
  p <- length(aa$bp)
  n <- length(xnum) / p
  stopifnot(n==round(n))
  stopifnot("diagonal" %in% names(aa))
  if (aa$diagonal) stopifnot(length(zz) == n * (choose(p,2) + p))
  else stopifnot(length(zz) == n * choose(p,2))

  out <- .C("compute_yhat_zz_R",
            xnum,
            as.integer(n),
            as.integer(p),
            zz,
            as.integer(aa$diagonal),
            as.double(aa$th),
            aa$bp,
            aa$bn,
            yhat=rep(0, n),
            PACKAGE="CRTConjoint")
  out$yhat
}

