library(rgl)
library(ggplot2)

# Define the target PDF (2-dimensional Rastrigin function within range [-2,2]*[-2,2])
pdf = function(x) {
  (20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])))/362.6667
}

# View the shape of the target PDF
curv_mesh = expand.grid(x=seq(-2,2,length.out=1e2),y=seq(-2,2,length.out=1e2),z=NA)
for (i in 1:nrow(curv_mesh)) {
  curv_mesh$z[i] = pdf(c(curv_mesh$x[i],curv_mesh$y[i]))
}
plot3d(x=curv_mesh$x,y=curv_mesh$y,z=curv_mesh$z)

# Sample
s = mhsample(c(0,0), 1e4, pdf=pdf, proppdf=function(x,y) 1/16, proprnd=function(x) c(runif(1,-2,2),runif(1,-2,2)))

# View the sampling result
ggplot() +
  geom_point(data = data.frame(x=s[,1],y=s[,2]), mapping = aes(x=x,y=y))
