# Metropolis-Hasting sample (like MATLAB)
# Input:
#   start: start value of the Markov chain.
#   nsample: number of samples to be generated.
#   pdf: the target distribution density.
#   logpdf: the log of "pdf".
#   proppdf: the proposal distribution density.
#   logproppdf: the log of "proppdf".
#   proprnd: the proposal distribution sampler.
#   symmetric: a logical value indicating whether the proposal distribution is symmetric. The default value is FALSE. If it is TRUE (i.e., proppdf is symmetric), proppdf is optional.
#   burnin: burnin number. The default value is 0.
#   thin: thin number. The default value is 1.
#   show: a logical value indicating whether to show the progress bar. The progress bar will greatly increase the running time of this function, which can be alleviated by using a large "showinterval" argument.
#   showinterval: an integer indicating the sample interval to update the progress bar.
# Output:
#   smpl: generated samples.
mhsample <- function(start, nsample, ...) {
  # Parse input
  varargin <- list(...)
  pdf <- varargin$pdf
  logpdf <- varargin$logpdf
  proppdf <- varargin$proppdf
  logproppdf <- varargin$logproppdf
  proprnd <- varargin$proprnd
  sym <- if (is.null(varargin$symmetric)) FALSE else varargin$symmetric
  burnin <- if (is.null(varargin$burnin)) 0 else varargin$burnin
  thin <- if (is.null(varargin$thin)) 1 else varargin$thin
  show <- if (is.null(varargin$show)) FALSE else varargin$show
  showinterval <- if (is.null(varargin$showinterval)) 1 else varargin$showinterval

  if (!is.array(start)) {
    start <- array(start, c(1,length(start)))
  }
  nDim <- ncol(start)
  if (is.null(logpdf))
    logpdf <- function(x) log(pdf(x))
  if (is.null(logproppdf))
    logproppdf <- function(x,y) log(proppdf(x,y))
  
  # Perform Metropolis-Hasting algorithm
  smpl <- array(0, dim = c(nsample, nDim))
  x0 <- start
  accept <- 0
  U <- log(runif(nsample * thin + burnin))
  if (show){
    library(tcltk)
    pb <- tkProgressBar(title = sprintf('%d / %d', 0, nsample * thin + burnin),
                        label = sprintf('Sampling: %d / %d', 0, nsample * thin + burnin),
                        max = nsample * thin + burnin, initial = 0)
  }
  for (i in (1 - burnin):(nsample * thin)) {
    y <- array(proprnd(x0), c(1,nDim))
    if (sym) {
      rho <- logpdf(y) - logpdf(x0)
    } else {
      rho <- logpdf(y) - logpdf(x0) + logproppdf(x0, y) - logproppdf(y, x0)
    }
    Ui <- U[i + burnin]
    acc <- Ui <= min(0, rho)
    x0[acc,] <- y[acc,]
    accept <- accept + acc

    if (i>0 && i%%thin==0) { # burnin and thin
      smpl[i/thin,] <- x0
    }
    if (show && (i+burnin)%%showinterval==0) {
      if (!environment(pb[["getVal"]])[[".killed"]])
        setTkProgressBar(pb, value = i + burnin, title = sprintf('%d / %d', i + burnin, nsample * thin + burnin), label = sprintf('Sampling: %d / %d', i + burnin, nsample * thin + burnin))
    }
  }
  
  # Clean and Return
  if (show) {
    close(pb)
  }
  
  return(smpl)
}
