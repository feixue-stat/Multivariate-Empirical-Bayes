library(glmnet)
membayesCM <- function (x0, y0, olsbeta=NULL, ols_intercept=NULL, scov=NULL, eta_pre=NULL, itmax=1000, prob0=rep(0.5,2), eps=1e-4, eps2=1e-8, standardize=FALSE, standardize.x=FALSE, sigma_error_pre_lower=0) {
  # olsbeta: OLS beta estimates for all tissues
  # scov: X^TX sample covariance matrix
  # itmax: Number of iterations
  # prob0: initial values of prior probabilities
  
  TT=dim(y0)[2]      # Total number of tissues
  p=dim(x0)[2]       # Total number of predictors
  N=dim(y0)[1]       # Number of total samples
  
  centerx=rep(0,p)
  scalex=rep(1,p)
  centery=rep(0,TT)
  scaley=rep(1,TT)
  x=x0
  y=y0
  if (standardize) {
    centerx=apply(x0,2,mean)
    for (j in 1:TT) {
      centery[j]=mean(y0[,j], na.rm=TRUE)
      y[,j]=(y0[,j]-centery[j])/scaley[j]
    }
    for (j in 1:p) {
      x[,j]=(x0[,j]-centerx[j])/scalex[j]
    }
  } 
  if (standardize.x) {
    centerx=apply(x0,2,mean)
    for (j in 1:p) {
      x[,j]=(x0[,j]-centerx[j])/scalex[j]
    }
  }
  intercept=rep(0,TT)
  
  obs=!is.na(y)           # Observaed patterns
  obs_num=apply(obs,2,sum)
  y2=y
  y2[is.na(y)]=0
  vary=apply(y2^2,2,sum)
  
  if (is.null(scov)) {
    scov=t(x)%*%x
  }
  invscov=solve(scov)     # (X^TX)^{-1}
  H=x%*%invscov%*%t(x)
  I=diag(1,N,N)
  temp_invers=array(0, dim = c(N,N,TT))
  
  posttau=matrix(0,TT,2)  # Posterior probabilities
  posttau[,1]=prob0[1]
  posttau[,2]=prob0[2]
  prob=prob0
  Gpdf=matrix(0,TT,2)
  logG=matrix(0,TT,2)
  Q_pre=-Inf
  
  
  if (is.null(olsbeta) | is.null(ols_intercept)) {
    olsbeta=matrix(NA,TT,p)
    ols_intercept=rep(NA,TT)
    for (t in 1:TT) {
      model1=lm(y[obs[,t],t] ~ x[obs[,t],])
      olsbeta[t,]=model1$coefficients[-1]
      ols_intercept[t]=model1$coefficients[1]
      
      if (sum(apply(as.matrix(x[obs[,t],]), 2, sd)==0)>0) {
        out <- tryCatch(cv.glmnet(y = y[obs[,t],t], x = x[obs[,t],,drop=FALSE], family = "gaussian") , error = function(e) e)
        if (!any(class(out) == "error")) {
          fit <- out$glmnet.fit
          olsbeta[t,]=t(as.numeric(coef(fit, s=out$lambda.min))[-1])
          ols_intercept[t]=as.numeric(coef(fit, s=out$lambda.min))[1]
        } else {
          olsbeta[t,is.na(olsbeta[t,])]=0
        }
      } else if (min(eigen(cor(x[obs[,t],,drop=FALSE]))$values)<0.01 | is.na(sum(model1$coefficients))) {
        out <- tryCatch(cv.glmnet(y = y[obs[,t],t], x = x[obs[,t],,drop=FALSE], family = "gaussian") , error = function(e) e)
        if (!any(class(out) == "error")) {
          fit <- out$glmnet.fit
          olsbeta[t,]=t(as.numeric(coef(fit, s=out$lambda.min))[-1])
          ols_intercept[t]=as.numeric(coef(fit, s=out$lambda.min))[1]
        } else {
          olsbeta[t,is.na(olsbeta[t,])]=0
        }
      }
    }
  }
  
  
  if (!is.null(olsbeta) & !is.null(ols_intercept)) {
    weight1=posttau[,1]/sum(posttau[,1])
    weight2=posttau[,2]/sum(posttau[,2])
    temp_sigma=rep(0,TT)
    for (t in 1:TT) {
      temp_sigma[t]=olsbeta[t,]%*%(scov/p)%*%olsbeta[t,]
    }
    sigma_error_pre1=sum(weight1*temp_sigma)
    
    sigma_error_pre2=0
    for (t in 1:TT) {
      if (sum(is.na(olsbeta[t,]))>0) {
        sigma_error_pre2=sigma_error_pre2+mean((y[obs[,t],t]-mean(y[obs[,t],t]))^2)
      } else {
        sigma_error_pre2=sigma_error_pre2+mean((y[obs[,t],t]-cbind(1,x[obs[,t],,drop=FALSE])%*%c(ols_intercept[t],olsbeta[t,]))^2)
      }
    }
    sigma_error_pre=max(min(sigma_error_pre1,sigma_error_pre2/TT,5), sigma_error_pre_lower)
  } else {
    sigma_error_pre=0
    for (t in 1:TT) {
      sigma_error_pre=sigma_error_pre+mean((y[obs[,t],t]-mean(y[obs[,t],t]))^2)
    }
    sigma_error_pre=max(min(sigma_error_pre/TT,5), sigma_error_pre_lower)
  }
  
  sigma_error_initial=sigma_error_pre
  
  if (!is.null(olsbeta)) {
    beta_pre=t(olsbeta)%*%posttau[,2]/sum(posttau[,2])
  } else {
    beta_pre=Inf
  }
  
  if (is.null(eta_pre)) {
    eta_pre=N*sigma_error_pre
  }
  
  posttau_pre=Inf
  prob_pre=Inf
  
  for (i in 1:itmax) {
    tempX=matrix(0,p,p)
    tempy=rep(0,p)
    for (t in 1:TT) {
      temp_invers[obs[,t], obs[,t], t]=solve(sigma_error_pre*I[obs[,t], obs[,t]]+eta_pre*H[obs[,t], obs[,t]])
      tempX = tempX + posttau[t,2]*t(x[obs[,t],])%*%temp_invers[obs[,t], obs[,t], t]%*%x[obs[,t],]
      tempy = tempy + posttau[t,2]*t(x[obs[,t],])%*%temp_invers[obs[,t], obs[,t], t]%*%y[obs[,t],t]
    }
    out <- tryCatch(solve(tempX) , error = function(e) e)
    if (any(class(out) == "error")) {
      beta=beta_pre
    } else {
      beta=solve(tempX)%*%tempy
    }
    
    
    for (t in 1:TT) {
      cov1=sigma_error_pre*I[obs[,t], obs[,t]]
      inverscov1=I[obs[,t], obs[,t]]/sigma_error_pre
      cov2=sigma_error_pre*I[obs[,t], obs[,t]]+eta_pre*H[obs[,t], obs[,t]]
      
      logG[t,1]=logpdfG(y[obs[,t],t], rep(0,obs_num[t]), cov1, inverscov1)
      logG[t,2]=logpdfG(y[obs[,t],t], as.numeric(x[obs[,t],]%*%beta), cov2)
      
      tempr=exp(logG[t,2]-logG[t,1])
      if ((prob[1]<0.1) & (prob[1]<tempr)) {
        posttau[t,1]=(prob[1]/tempr)/(prob[2]+prob[1]/tempr)
      } else {
        posttau[t,1]=1/(1+(prob[2]/prob[1])*tempr)
      }
      posttau[t,2]=1-posttau[t,1]
    }
    if (sum(posttau[,1])==0) {
      prob[1]=0.01
      prob[2]=1-prob[1]
    } else if (sum(posttau[,2])==0) {
      prob[2]=0.01
      prob[1]=1-prob[2]
    } else {
      prob[1]=sum(posttau[,1])/sum(posttau)
      prob[2]=sum(posttau[,2])/sum(posttau)
    }
    
    
    res=matrix(0,N,t)
    upper1=sum(posttau[,1]*vary)
    lower1=0
    for (t in 1:TT) {
      res[obs[,t],t] = y[obs[,t],t] - x[obs[,t],] %*% beta
      upper1 = upper1 + sigma_error_pre^2*posttau[t,2]*t(res[obs[,t],t]) %*% temp_invers[obs[,t], obs[,t], t] %*% temp_invers[obs[,t], obs[,t], t] %*% res[obs[,t],t]
      lower1 = lower1 + sigma_error_pre*posttau[t,2]*sum(diag(temp_invers[obs[,t], obs[,t], t]))
    }
    sigma_error=as.numeric(upper1/(lower1+sum(posttau[,1]*obs_num)))
    
    if (sum(posttau[,2])==0) {
      eta=eta_pre
    } else {
      upper2=0
      lower2=0
      for (t in 1:TT) {
        upper2 = upper2 + eta_pre^2*posttau[t,2]*t(res[obs[,t],t]) %*% temp_invers[obs[,t], obs[,t], t] %*% H[obs[,t], obs[,t]] %*% temp_invers[obs[,t], obs[,t], t] %*% res[obs[,t],t]
        lower2 = lower2 + eta_pre*posttau[t,2]*sum(diag(temp_invers[obs[,t], obs[,t], t]%*%H[obs[,t], obs[,t]]))
      }
      eta=as.numeric(upper2/lower2)
    }
    
    if (sum(posttau[,1])==0) {
      if (sum(is.infinite(logG[,2]))>0) {
        Q=-Inf
      } else {
        Q=sum(posttau[,2]*logG[,2])+log(prob[2])*sum(posttau[,2])
      }
    } else if (sum(posttau[,2])==0) {
      if (sum(is.infinite(logG[,1]))>0) {
        Q=-Inf
      } else {
        Q=sum(posttau[,1]*logG[,1])+log(prob[1])*sum(posttau[,1])
      }
    } else {
      if (sum(is.infinite(logG))>0) {
        Q=-Inf
      } else {
        Q=sum(posttau*logG)+log(prob[1])*sum(posttau[,1])+log(prob[2])*sum(posttau[,2])
      }
    }
    
    diff=(sigma_error-sigma_error_pre)^2+sum((eta-eta_pre)^2)+sum((beta-beta_pre)^2)+sum((posttau-posttau_pre)^2)+sum((prob-prob_pre)^2)
    
    if (abs(Q_pre)<1 | is.infinite(Q_pre)) {
      Q_diff=(Q-Q_pre)^2
    } else {
      Q_diff=((Q-Q_pre)/Q_pre)^2
    }

    if (is.na(Q_diff<eps2 & diff<eps)){
      print(Q_diff<eps2 & diff<eps)
    }
    
    if (Q_diff<eps2 & diff<eps) {
      break
    }
    Q_pre=Q
    
    sigma_error_pre=sigma_error
    eta_pre=eta
    beta_pre=beta
    posttau_pre=posttau
    prob_pre=prob
  }
  if (is.nan(Q)) {
    print('Error')
  }
  converge=TRUE
  if (i==itmax) {
    print('May not converge')
    converge=FALSE
  }
  postbeta=matrix(0, TT, p)
  for (t in 1:TT) {
    postbeta[t,]=posttau[t,2]*solve((t(x)%*%x)/eta+(t(x[obs[,t],])%*%x[obs[,t],])/sigma_error)%*%((t(x)%*%x)%*%beta/eta+(t(x[obs[,t],])%*%y[obs[,t],t])/sigma_error)
    intercept[t]=centery[t]-crossprod(postbeta[t,],as.numeric(centerx))
  }
  BF=exp(logG[,1]-logG[,2])
  returnlist=list('postbeta'=postbeta, 'beta'=beta, 'posttau'=posttau, 'sigma_error'=sigma_error, 
                  'eta'=eta, 'prob'=prob, 'converge'=converge, 'i'=i, 
                  'centerx'=centerx, 'centery'=centery, 'scalex'=scalex, 'scaley'=scaley, 
                  'intercept'=intercept, 'BF'=BF)
  return(returnlist)
}


Gassianpdf <- function (x, mean, cov) {
  k=length(mean)
  pdf=exp(-0.5*(x-mean)%*%solve(cov)%*%(x-mean))/(sqrt((2*pi)^k*det(cov)))
  return(pdf)
}

logpdfG <- function (x, mean, cov, inverscov=NULL) {
  k=length(mean)
  if (is.null(inverscov)){
    inverscov=solve(cov)
  }
  logpdf=-0.5*(x-mean)%*%inverscov%*%(x-mean)-0.5*(k*log(2*pi)+as.numeric(determinant(cov)$modulus))  #log(det(cov)))
  return(logpdf)
}
