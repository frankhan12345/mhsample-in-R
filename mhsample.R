# Metropolis-Hasting sample (like MATLAB)
# Input:
#   start: start value of the Markov chain.
#   nsample: number of samples to be generated.
#   pdf: the target distribution density.
#   proppdf: the proposal distribution density.
#   proprnd: the proposal distribution sampler.
#   symmetric: a logical value indicating whether the proposal distribution is symmetric. The default value is FALSE. If it is TRUE (i.e., proppdf is symmetric), proppdf is optional.
#   burnin: burnin number. The default value is 0.
#   thin: thin number. The default value is 1.
#   show: a logical value indicating whether to show the real-time progress
# Output:
#   smpl: generated samples
mhsample = function(start, nsample, ...) {
  # Parse input
  varargin = list(...)
  pdf = varargin$pdf
  proppdf = varargin$proppdf
  proprnd = varargin$proprnd
  sym = if (is.null(varargin$symmetric)) FALSE else varargin$symmetric
  burnin = if (is.null(varargin$burnin)) 0 else varargin$burnin
  thin = if (is.null(varargin$thin)) 1 else varargin$thin
  show = if (is.null(varargin$show)) FALSE else varargin$show
  
  nDim = length(start)
  
  # Perform Metropolis-Hasting algorithm
  smpl = array(0, dim = c(burnin + nsample * thin, nDim))
  x0 = start
  i = 1
  if (show){
    library(tcltk)
    pb = tkProgressBar(max = nrow(smpl))
  }
  while (i <= nrow(smpl)) {
    y = proprnd(x0)
    if (sym) {
      rho = pdf(y) / pdf(x0)
    } else {
      rho = pdf(y) / pdf(x0) * proppdf(x0,y) / proppdf(y,x0)
    }
    acc = min(1,rho)
    u = runif(1)
    if (u<acc) {
      smpl[i,] = y
      x0 = y
      if (show) {
        if (!environment(pb[["getVal"]])[[".killed"]])
          setTkProgressBar(pb, value = i, title = sprintf('%d / %d', i, nrow(smpl)), label = sprintf('Sampling: %d / %d', i, nrow(smpl)))
      }
      i = i + 1
    }
  }
  
  # Burnin
  smpl = smpl[(1+burnin):nrow(smpl),]
  smpl = array(smpl, dim = c(NROW(smpl), NCOL(smpl)))
  
  # Thin
  smpl = smpl[seq(thin,nsample*thin,thin),]
  
  # Clean and Return
  if (show) {
    close(pb)
  }
  
  return(smpl)
}
