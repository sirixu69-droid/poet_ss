library(matrixStats) 
library(parallel)
library(MASS)
library("mvtnorm")
library("mnormt")
library(ICSNP)
library(ICS)
library("SpatialNP")
library("MNM")
library(flare)
#library(glasso)
library(glassoFast)
library(stats)
norm2=function(x){sqrt(sum(x^2))}

tt1<-proc.time()[1]

result<-matrix(0,4,10)
for (ip in 1:4)
{
n<-250
p<-100*(ip+1)
m<-3
mi=4
SIM<-1000
emER<-rep(0,SIM)
emGR<-rep(0,SIM)

for (ii in 1:SIM)
{

  indices <- matrix(1:p, nrow = p, ncol = p)
  Sigma_u <- 0.9^abs((indices - t(indices)))
  Sigma_joint <- matrix(0, nrow = m + p, ncol = m + p)
  Sigma_joint[1:m, 1:m] <- diag(m)  
  Sigma_joint[(m+1):(m+p), (m+1):(m+p)] <- Sigma_u  
  s = c(1, 0.75^2, 0.5^2)
  B <- matrix(0, nrow = p, ncol = m)
  for (k in 1:m) {
    B[, k] <- rnorm(p, mean = 0, sd = sqrt(s[k]))
  }

  Sigma <- B %*% t(B) + Sigma_u
  c <- p / sum(diag(Sigma))
  Sigma0 <- c * Sigma
  Sigma_u0 <- c * Sigma_u
if(mi==1){
    joint_sample <- mvtnorm::rmvnorm(n, mean = rep(0, m + p), sigma = Sigma_joint)
  }     
if (mi==2) {
    df=4
    scale_factor <- sqrt(df / (df - 2))
    joint_sample <- mvtnorm::rmvt(n, sigma = Sigma_joint, df = 4)
    joint_sample <- joint_sample/scale_factor
   }
if (mi==3) {
    df=2.2
    scale_factor <- sqrt(df / (df - 2))
    joint_sample <- mvtnorm::rmvt(n, sigma = Sigma_joint, df = 2.2)
    joint_sample <- joint_sample/scale_factor
   }

if (mi==4)
{ 
   df<-10
   joint_sample<-rmvnorm(n,rep(0,m+p),Sigma_joint)
   jn<-rbinom(n,1,0.8)
   joint_sample<-jn*joint_sample+(1-jn)*joint_sample*df
   scale_factor<-0.8+df^2*0.2
   joint_sample <- joint_sample/scale_factor
}


  f_t <- joint_sample[, 1:m, drop = FALSE]
  u_t <- joint_sample[, (m + 1):(m + p), drop = FALSE]

  Y <- f_t %*% t(B) + u_t
  Y <- sqrt(c) * Y


sy<-SCov(Y)

syeig<-eigen(sy)$values

M<-10

emER[ii]<-which.max(syeig[1:(M-1)]/syeig[2:M])

mnp<-min(n,p)

seig<-sum(syeig[1:mnp])

vj<-rep(0,M)
for (i in 1:(M))
{
vj[i]<-sum(syeig[i:mnp])
}


emGR[ii]<-which.max(log(1+syeig[1:(M-1)]/vj[1:(M-1)])/log(1+syeig[2:M]/vj[2:M]))

}

hist(emER)

hist(emGR)

result[ip,1:5]<-tabulate(emER,nbins=5)/SIM
result[ip,6:10]<-tabulate(emGR,nbins=5)/SIM

}

result

write.table(result,sep="&","num4.txt")


tt2<-proc.time()[1]
tt2-tt1


