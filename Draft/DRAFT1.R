## 20181002
## 1. update A
sylvester=function(A,B,C){
  I_A=diag(1,dim(A)[1],dim(A)[2])
  I_B=diag(1,dim(B)[1],dim(B)[2])
  M=kronecker(I_B,A) + kronecker(t(B),I_A)
  L=matrix(c(C),dim(A)[1]*dim(B)[1],1)
  cx=ginv(M) %*% L
  X=matrix(cx,dim(A)[1],dim(B)[1])
  return(X)
}

L_num=function(n){
  library(tidyr)
  L=matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      L[i,j] = paste(i,',',j)
    }
  }
  l=c(L[upper.tri(L)])
  ll=data.frame(l)
  
  LL=ll %>% separate(l,c('l1','l2'),sep=',')
  LL=data.frame(l1=as.numeric(LL[,'l1']),l2=as.numeric(LL[,'l2']))
  return(LL)
}

X=matrix(1:(20),5,4)

A=X

n=dim(X)[1]; p=dim(X)[2]
eplison_p=L_num(p)
eplison_n=L_num(n)

nu1=1; nu2=1
lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,4,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,5,dim(eplison_p)[1])

#update_A = function(X, nu1, nu2, lambda_1l, lambda_2k, v, z){
library(MASS)
n=dim(X)[1]; p=dim(X)[2]
eplison_p=L_num(p)
Ep=matrix(0,p,p)
for(i in 1:dim(eplison_p)[1]){
  l1=eplison_p[i,'l1']
  l2=eplison_p[i,'l2']
  el1=matrix(rep(0,p),p,1)
  el2=matrix(rep(0,p),p,1)
  el1[l1,1]=1;el2[l2,1]=1
  Ep=Ep+(el1-el2) %*% t(el1-el2)
}
M = diag(1,p,p) + nu1 * Ep

eplison_n=L_num(n)
En=matrix(0,n,n)
for(i in 1:dim(eplison_n)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  el1=matrix(rep(0,n),n,1)
  el2=matrix(rep(0,n),n,1)
  el1[l1,1]=1;el2[l2,1]=1
  En=En+(el1-el2) %*% t(el1-el2)
}

N = nu2 * En 
C2 = matrix(0,n,p) 

for(i in 1:dim(eplison_p)[1]){
  l1=eplison_p[i,'l1']
  l2=eplison_p[i,'l2']
  el1=matrix(rep(0,p),p,1)
  el2=matrix(rep(0,p),p,1)
  el1[l1,1]=1;el2[l2,1]=1
  C2 = C2+(lambda_1[,i] + nu1 * v[,i]) %*% t(el1-el2)
}

C3 = matrix(0,n,p) 

for(i in 1:dim(eplison_n)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  el1=matrix(rep(0,n),n,1)
  el2=matrix(rep(0,n),n,1)
  el1[l1,1]=1;el2[l2,1]=1
  C3 = C3+(el1-el2) %*% t(lambda_2[,i] + nu2 * z[,i])
}

C = X +  C2 + C3  

A = sylvester(N,M,C)

return(A)
}

update_A = function(X, nu1, nu2, lambda_1, lambda_2, v, z){
  library(MASS)
  n=dim(X)[1]; p=dim(X)[2]
  eplison_n=L_num(n)
  eplison_p=L_num(p)
  
  En=matrix(0,n,n)
  Ep=matrix(0,p,p)
  
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    el1=matrix(rep(0,n),n,1)
    el2=matrix(rep(0,n),n,1)
    el1[l1,1]=1;el2[l2,1]=1
    En=En+(el1-el2) %*% t(el1-el2)
  }
  
  for(i in 1:dim(eplison_p)[1]){
    l1=eplison_p[n,'l1']
    l2=eplison_p[n,'l2']
    el1=matrix(rep(0,p),p,1)
    el2=matrix(rep(0,p),p,1)
    el1[l1,1]=1;el2[l2,1]=1
    Ep=Ep+(el1-el2) %*% t(el1-el2)
  }
  
  M = diag(1,n,n) + nu1 * En
  
  N = nu2 * Ep
  
  C2 = matrix(0,n,p) 
  
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    el1=matrix(rep(0,n),n,1)
    el2=matrix(rep(0,n),n,1)
    el1[l1,1]=1;el2[l2,1]=1
    C2 = C2 + (el1-el2) %*% t(lambda_1[,i] + nu1 * v[,i])
  }
  
  C3 = matrix(0,n,p) 
  
  for(i in 1:dim(eplison_p)[1]){
    l1=eplison_p[i,'l1']
    l2=eplison_p[i,'l2']
    el1=matrix(rep(0,p),p,1)
    el2=matrix(rep(0,p),p,1)
    el1[l1,1]=1;el2[l2,1]=1
    C3 = C3+(lambda_2[,i] + nu2 * z[,i]) %*% t(el1-el2)
  }
  
  C = X +  C2 + C3  
  
  A = sylvester(M,N,C)
  
  return(A)
}

A = update_A(X,nu1, nu2, lambda_1, lambda_2, v, z)

# 3 update lambda
update_lambda=function(X, A, nu1, nu2,v, z){
  n = dim(X)[1]; p = dim(X)[2]
  eplison_p = L_num(p)
  eplison_n = L_num(n)
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    a_l1 = matrix(A[l1,],p,1)
    a_l2 = matrix(A[l2,],p,1)
    lambda_1[,i] = lambda_1[,i] + nu1 * (v[,i] - a_l1 + a_l2)
  }
  for(i in 1:dim(eplison_p)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    a_k1 = matrix(A[,l1],n,1)
    a_k2 = matrix(A[,l2],n,1)
    lambda_2[,i] = lambda_2[,i] + nu2 * (z[,i] - a_k1 + a_k2)
  }
  return(lambda=list(lambda_1,lambda_2))
}


lambda = update_lambda(X, A, nu1, nu2,v, z)

lambda_1 = lambda[[1]]
lambda_2 = lambda[[2]]

## 2. update v and z

prox = function(v,sigma,n=2){
  if(n == 2){
    return(max(1-sigma/sum(v*v),0) %*% v)
  }
}

sigma_2k = gamma_2 * u_k/nu2
sigma_1l = gamma_1 * w_l/nu1

for(i in 1:dim(eplison_n)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  
  a_l1 = A[l1,]; a_l2 = A[l2,]
  v_temp = a_l1 - a_l2 + 1/nu1 * lambda_1[,i]
  w_l = exp(-0.5 * (t(X[l1,] - X[l2,]) %*% (X[l1,] - X[l2,])))
  
  sigma_1l = gamma_1 * w_l/nu1
  
  v[,i] = prox(v_temp,sigma_1l)
}
for(i in 1:dim(eplison_p)[1]){
  l1=eplison_p[i,'l1']
  l2=eplison_p[i,'l2']
  
  a_l1 = A[,l1]; a_l2 = A[,l2]
  v_temp = a_l1 - a_l2 + 1/nu2 * lambda_2[,i]
  
  u_k = 1
  sigma_2k = gamma_2 * u_k/nu2
  
  z[,i] = prox(v_temp,sigma_2k)
}




strassenInv <- function(A){
  
  div4 <- function(A, r){
    A <- list(A)
    A11 <- A[[1]][1:(r/2),1:(r/2)]
    A12 <- A[[1]][1:(r/2),(r/2+1):r]
    A21 <- A[[1]][(r/2+1):r,1:(r/2)]
    A22 <- A[[1]][(r/2+1):r,(r/2+1):r]
    A <- list(X11=A11, X12=A12, X21=A21, X22=A22)
    return(A)
  }
  
  if (nrow(A) != ncol(A)) 
  { stop("only square matrices can be inverted") }
  
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  if ( (is.wholenumber(log(nrow(A), 2)) != TRUE) || (is.wholenumber(log(ncol(A), 2)) != TRUE) )
  { stop("only square matrices of dimension 2^k * 2^k can be inverted with Strassen method") }
  
  A <- div4(A, dim(A)[1])
  
  R1 <- solve(A$X11)
  R2 <- A$X21 %*% R1
  R3 <- R1 %*% A$X12
  R4 <- A$X21 %*% R3
  R5 <- R4 - A$X22
  R6 <- solve(R5)
  C12 <- R3 %*% R6
  C21 <- R6 %*% R2
  R7 <- R3 %*% C21
  C11 <- R1 - R7
  C22 <- -R6
  
  C <- rbind(cbind(C11,C12), cbind(C21,C22))
  
  return(C)
}


## combine 