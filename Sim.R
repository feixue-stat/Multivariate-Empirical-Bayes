library(MASS)
p=30   #Number of SNPs
N=100   #Number of samples, N should be larger than p under low-dimensional setting
ntrain=50
ntest=50
nval=10
nt=50    #Number of tissues
q1=10
beta1=c(0.01, 0.1, 0.5, 1)
q2=10
beta2=c(0.005, 0.05, 0.25, 0.5)

a=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
prob=0.5    #Prior probability for I=1
sigmasq=1    #variance for error term

source('membayes.R')

set.seed(1) 

s1=1
b=1

sigma1=matrix(0,p,p)          # Covariance matrix for x
sigma1[1:p,1:p]=diag(1-a[s1],p,p)+matrix(a[s1],p,p)
sigma2=matrix(0,p,p)          # Covariance matrix for beta_t
sigma2[1:p,1:p]=(diag(1-a[s1],p,p)+matrix(a[s1],p,p))
beta=c(beta1[b]*rep(1,q1), beta2[b]*rep(1,q2), rep(0, p-q1-q2))

x=mvrnorm(n=N,rep(0,p),sigma1)
y=matrix(0,N,nt)

train.inds <- sample(1:N, ntrain)
X.train <- x[train.inds,]
X.test <- x[-train.inds,]

scov=t(X.train)%*%X.train
invscov=solve(scov)
beta_t=matrix(0,nt,p)
I_t=sample(c(0,1), nt, prob = c(1-prob,prob), replace = TRUE)
for (j in 1:nt) {
  if (I_t[j]==1) {
    beta_t[j,]=mvrnorm(n=1,beta,sigma2)
    noise=rnorm(n=N,0,sigmasq)
    y[,j]=x%*%beta_t[j,]+noise
  } else {
    beta_t[j,]=rep(0,p)
    noise=rnorm(n=N,0,sigmasq)
    y[,j]=x%*%beta_t[j,]+noise
  }
}

Y.train <- y[train.inds,]
Y.test <- y[-train.inds,]

# -----------------------------------------------------------
# OLS method
# -----------------------------------------------------------
olsbeta=matrix(0,nt,p)
for (j in 1:nt) {
  olsbeta[j,]=invscov%*%t(X.train)%*%Y.train[,j]
}


# -----------------------------------------------------------
# Proposed method
# -----------------------------------------------------------
model=membayesCM(X.train, Y.train, olsbeta)
beta_em=model$postbeta





