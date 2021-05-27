# Metropolis-Hasting sample (like MATLAB)
# Input:
#   start: start value of the Markov Chain
#   nsample: number of samples to be generated
#   pdf: the target distribution density
#   proppdf: the proposal distribution density
#   proprnd: the proposal distribution sampler
#   burnin: burnin number
#   thin: thin number
# Output:
#   smpl: generated samples
mhsample = function(start, nsample, ...) {
  # Parse input
  varargin = list(...)
  pdf = varargin$pdf
  proppdf = varargin$proppdf
  proprnd = varargin$proprnd
  burnin = if (is.null(varargin$burnin)) 0 else varargin$burnin
  thin = if (is.null(varargin$thin)) 1 else varargin$thin
  
  
  nDim = length(start)
  
  # Perform Metropolis-Hasting algorithm
  smpl = array(0, dim = c(burnin + nsample * thin, nDim))
  x0 = start
  i = 1
  while (i <= nrow(smpl)) {
    y = proprnd(x0)
    rho = pdf(y) / pdf(x0) * proppdf(x0,y) / proppdf(y,x0)
    acc = min(1,rho)
    u = runif(1)
    if (u<acc) {
      smpl[i,] = y
      x0 = y
      i = i + 1
    }
  }
  
  # Burnin and Thin
  smpl = smpl[(1+burnin):nrow(smpl),]
  smpl = smpl[seq(thin,nsample*thin,thin),]
  
  return(smpl)
}
