
library(MASS)
library(vars)
#library(parcor)


###################################################################
# VAR(p) GLS estimation with constraints
#
# Constraints are given by matrix J of dimension k x kp with 0 or 1 entries 
# depending on whether the elements of [A1 A2 ... Ap] are included or not
#
# A multivariate series y is of dimension k x T
#
###################################################################

VAR.const.gls = function(y, p, J, szInv){
    k <- dim(y)[1]
    if(missing(szInv)){ szInv = diag(k) }

    y <- y - rowMeans(y)
    T <- dim(y)[2]
    T1 <- T-p

# create vector X1 and Y1
    X1 <- mat.or.vec(k*p, T1)
    Y1 <- mat.or.vec(k, T1)
    for(j in 1:T1){
	      id <- seq(from=j+p-1, to=j, by=-1)
	      x <- as.vector(y[,id])
        X1[,j] <- x
        Y1[,j] <- y[,(j+p)]
    }

# create constraint matrix R
    vecA1 <- c(J)
    A1 <- matrix(vecA1, nrow=k)
    k2 <- sum(vecA1)
	  R <- matrix(0, k^2*p, k2)

	  for(i in 1:(k^2*p)){
	      if(vecA1[i] == 1){ 
	        id <- sum(vecA1[1:i])  
	        R[i,id] <- 1 }  
	  }

# Do iteration to estimate A1 (this follows Lutkepohl, pp. 194-) 
    varA1 <- solve(t(R)%*%kronecker(X1%*%t(X1), szInv)%*%R)
    hA1old <- matrix(R%*% varA1 %*%t(R)%*%kronecker(X1, szInv)%*%as.vector(Y1), nrow=k)
    out.sig <- VAR.sigma(y, p, hA1old)
    if(out.sig$LL != Inf){
      szInv <- solve(out.sig$Sigma_z)
      diff <- 100 
      iter <- 1
        while( diff >= .001 & iter < 30){
          varA1 <- solve(t(R)%*%kronecker(X1%*%t(X1), szInv)%*%R)
          hatA1 <- R%*%varA1 %*%t(R)%*%kronecker(X1, szInv)%*%as.vector(Y1)
	
          hA1 <- matrix(hatA1, nrow=k)
          szInvold <- szInv
          out.sig <- VAR.sigma(y, p, hA1)
          szInv <- solve(out.sig$Sigma_z)
          diff <- sum((szInv - szInvold)^2)
          iter <- iter + 1
		    }
	  }

    out <- out.sig
    out$hatA <- hA1
    out$varA <- diag(varA1)
    out$T <- as.vector(hA1[hA1!=0])/ sqrt(diag(varA1))
    return(out)
}


###################################################################
# This estimates the Sigma matrix of the residuals and a few other things
#
# A multivariate series y is of dimension k x T
#
# One supposes that the rows of y have zero mean
#
###################################################################

VAR.sigma = function(y, p, A){
	  k <- dim(y)[1]
	  T <- dim(y)[2]
	  T1 <- T - p
	  m <- sum(A != 0) # counting the num of non-zero coeff
	  eps <- 1.0e-1
	  
    # Calculate residuals
    Resi <- mat.or.vec(k, T1)
    X <- mat.or.vec(k*p, T1)
    Y1 <- mat.or.vec(k, T1)
    for(j in 1:T1){
	    id <- seq(from= j+p-1, to = j, by=-1)
	    x <- as.vector(y[,id])
      Resi[,j] <-  y[,p+j]  - A%*%x
	  }

# Estimate Sigma
    Sigma_z <- Resi%*%t(Resi)/T1
    LL <- T1*log(det(Sigma_z)) 
    BIC <- LL + log(T1)*m

    out <- list()
    out$Resi <- Resi
    out$Sigma_z <- Sigma_z
    out$BIC <- BIC
    out$LL <- LL
    out$m <- m
    out$T1 <- T1
    return(out)
}


######################################
# 
#   Doing variable selection through t-statistics
#
#   p: VAR(p)
#
####################################

sVAR.tstat = function(y, p){

  k <- dim(y)[1]
  # no restriction to start with
  out.var <- VAR.const.gls(y, p=p, J=matrix(1,k,k*p))
  ln <- k*p

  init <- VAR.sigma(y, p=p,out.var$hatA)$Sigma_z
  szInv <- solve(init)

  absT <- abs(out.var$T) # need to fix this
  rt <- rank(-absT) 
  I <- matrix( rep(1:k, ln), nrow=k)
  J <- matrix( rep(1:ln, each=k), nrow=k)
  R1 <- matrix(0,k,k*p)
  kk <- sum(out.var$hatA != 0)
  LL2 <- BIC2 <- numeric(kk)

  for(i in 1:kk){
	  id <- which(rt == i)
    R1[I[id], J[id]] <- 1
	  lse.out <- VAR.const.gls(y, p=p, R1, szInv)
	  LL2[i] <- lse.out$LL
	  BIC2[i] <- lse.out$BIC
  }
  
	minid <- which.min(BIC2)
  R2 <- matrix(0,k,k*p)
	for(j in 1:minid){
	  id <- which(rt == j)
    R2[I[id], J[id]] <- 1
	}

	tstat.out <- VAR.const.gls(y, p=p, R2, szInv)

  tstat.out$BIC2 <- BIC2
  tstat.out$LL2 <- LL2
  tstat.out$order <- p
  return(tstat.out)
  
}


#############################################
# sparse VAR with regular Lasso
##############################################

sVAR.lasso = function(data, p, nf){

if(missing(nf)){ nf = 10 };
if(missing(p)){ p = 1 };

  y = data - rowMeans(data);
  T1 = dim(y)[2]-p;
  k = dim(y)[1];
  # create vector X1 and Y1
  X1 = matrix(0, k*p, T1); 
  Y1 = matrix(0, k, T1);
  for(j in 1:T1){
    # ar term
    id = seq(from= j+p-1, to = j, by=-1);
    x = as.vector(y[,id]);
    X1[,j] = x;
    Y1[,j] = y[,(j+p)];
  }
  ty = t(y); ty = data.frame(ty);
  out = VAR(ty, type="none", p=p);
  Sig= summary(out)$covres;
  sv=svd(Sig)
  hal = sv$u%*%diag(1/sqrt(sv$d))%*%t(sv$v)
  y1 = as.vector(Y1);
  x1 = kronecker(t(X1), diag(1, k));

  diff =100; iter=1;
  while( diff >= .01 & iter < 5){
    adjSig = kronecker(diag(1, T1), hal);
    Y2 = adjSig%*%y1;
    X2 = adjSig%*%x1;
    cvfit = cv.glmnet(X2, Y2, alpha = 1, intercept=TRUE, standardize=FALSE, type.measure = "mse", nfolds=nf);
    cf.cv = coef(cvfit, s = "lambda.min")
    cf.cv = cf.cv[-1];
    hA1 = matrix(cf.cv, nrow=k, byrow=FALSE);
    Signew= VAR.sigma(y, p, hA1)$Sigma_z;
    diff = sum((Sig - Signew)^2);
    Sig = Signew;
    iter = iter +1;
    sv=svd(Sig)
    hal = sv$u%*%diag(1/sqrt(sv$d))%*%t(sv$v)
}

out=list();
out$cv= cvfit;
out$hatA = hA1;
out$lam = cvfit$lambda.min;
out$Sigma_z = Sig;
return(out);
}

###########################################
# sparse VAR with adaptivelasso
###########################################

sVAR.adaplasso = function(y, p, nf){ 
  if(missing(nf)){ nf=10};
  y = y - rowMeans(y);
  k = dim(y)[1];
  T = dim(y)[2];
  T1 = T-p;
  
  # create vector X1 and Y1
  X1 = matrix(0, k*p, T1); 
  Y1 = matrix(0, k, T1);
  for(j in 1:T1){
    # ar term
    id = seq(from= j+p-1, to = j, by=-1);
    x = as.vector(y[,id]);
    X1[,j] = x;
    Y1[,j] = y[,(j+p)];
  }
  
  y1 = as.vector(Y1);
  x1 = kronecker(t(X1), diag(1, k));
  
  hA0 = kronecker(solve(X1%*%t(X1))%*%X1, diag(k))%*%y1;
  hA0 = matrix(hA0, nrow=k);
  Sig = VAR.sigma(y, p, hA0)$Sigma_z;
  sv = svd(Sig)
  hal = sv$u%*%diag(1/sqrt(sv$d))%*%t(sv$v)
  
  hhA0 = hA0;
  diff =100; iter=1;
  while( diff >= .01 & iter < 5){
    
    adjSig = kronecker(diag(1, T1), hal);
    Y2 = adjSig%*%y1;
    X2 = adjSig%*%x1;
    
    pp = adalasso(X2, Y2, k=nf, intercept=FALSE);
    hA2 = matrix(pp$coefficients.lasso, nrow=k);
    hA1 = matrix(pp$coefficients.adalasso, nrow=k);
    diff = sum((hA1 - hA0)^2);
    iter = iter +1;
    hA0 = hA1;
    Sig= VAR.sigma(y, p, hA0)$Sigma_z;
    sv=svd(Sig)
    hal = sv$u%*%diag(1/sqrt(sv$d))%*%t(sv$v)
  }
  
  out=list();
  out$cv= pp$cv.adalasso;
  out$hatA = hA1;
  out$lam = pp$lambda.adalasso;
  out$Sigma_z = Sig;
  ## standard LASSO result
  out$lasA = hA2;
  return(out);
}

##############################################################
# Generating VAR(p) model
##############################################################

VAR.sim = function(T, A, Sigma){
k = dim(A)[1];
p = dim(A)[2]/k;
burn = 500;

inno = mvrnorm(n=T+burn, rep(0, k), Sigma);
init = mvrnorm(n=p, rep(0, k), Sigma);
init = matrix(init, nrow=p);

# Find index for previous observations
	j=1;
	# ar term
	id = seq(from= j+p-1, to = j, by=-1);

Y = matrix(0, (T+burn), k);
for(r in 1:1:(T+burn)){
	Y[r,] = A%*%as.vector(t(init[id,])) + inno[r,];
      init = rbind(init[-1,], Y[r,]);
}

return(t(Y[-(1:burn),])) # Final data is k*T matrix
}

###########################################
# OLS estimation of VAR
###########################################

VAR.lse  = function(y, p){
y = y - rowMeans(y);
T = dim(y)[2];
T1 = T-p;
k = dim(y)[1];

# create vector X1 and Y1
    X1 = matrix(0,k*p, T1); 
    Y1 = matrix(0,k, T1);
    for(j in 1:T1){
	# ar term
	id = seq(from= j+p-1, to = j, by=-1);
	x = as.vector(y[,id]);
      X1[,j] = x;
      Y1[,j] = y[,(j+p)];
}

hatA = Y1%*%t(X1)%*%solve(X1%*%t(X1));

Resi = matrix(0,k, T1); 
# residuals
    for(j in 1:T1){
	id = seq(from= j+p-1, to = j, by=-1);
	x = as.vector(y[,id]);
      Resi[,j] =  y[,p+j]  - hatA%*%x;
	}
     Sigma_z = Resi%*%t(Resi)/T1;
     bic =  T1*log(det(Sigma_z)) + log(T1)* sum(hatA != 0);
return(list(hatA = hatA, Sigma=Sigma_z, bic=bic, p=p))
}


##############################################################
# Forecasting VAR(p) model
##############################################################

VAR.forecast = function(yf, h, A){

T1 = dim(yf)[2];
k = dim(A)[1];
p = dim(A)[2]/k;

# Find index for previous observations
id = seq(from= T1, to = T1-p+1, by=-1);

Y1 = Y = matrix(0, k, h);
for(r in 1:1:h){
	Y[,r] = A%*%as.vector(yf[,id]);
      yf = cbind(yf[,-1], Y[,r]);
	Y1[,r] = Y[,r];
}

return(Y1) # Final data is k*T matrix
}


################################
# VAR(p) best model 
# BIC is calculated based on RSS
################################

VAR.best = function(y, maxorder){
  ## Input data is dim*N
  N = ncol(y);
  dim = nrow(y);
  BIC = numeric(maxorder);
  for(p in 1:maxorder){
    ar1 = VAR(t(y), p=p, type="none")
    BIC[p] =  as.numeric(VAR.sigma(y, p=p, Bcoef(ar1))$BIC);
  }
  p.opt = which.min(BIC);
  best = VAR(t(y), p=p.opt, type="none")
  out = VAR.sigma(y, p=p.opt, Bcoef(best));
  out$hatA = Bcoef(best);
  out$BIC = BIC;
  out$opt = p.opt
  return(out)
}



