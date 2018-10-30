X = matrix(1:12,4,3)
n=dim(X)[1]; p=dim(X)[2]
eplison_p=L_num(p)
eplison_n=L_num(n)
nu1=1; nu2=1
lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])

library(MASS)
library(Matrix)
A= X
n=dim(X)[1]; p=dim(X)[2]

alk = function(A,n,p){
  count = 0
  el1 = el2 = matrix(0,n,n*(n-1)/2)
  for(i in 1:(n-1)){
    el1[,(count+1):(count+i)] = diag(1,n,i)
    temp = matrix(0,n,i)
    temp[i,] = 1
    el2[,(count+1):(count+i)] = temp
    count = count+i
  }
  el2 = rbind(el2[n,],el2[1:(n-1),])
  al1 = t(A) %*% el1; al2 = t(A) %*% el2
  
  # k
  count = 0
  ek1 = ek2 = matrix(0,p,p*(p-1)/2)
  for(i in 1:(p-1)){
    ek1[,(count+1):(count+i)] = diag(1,p,i)
    temp = matrix(0,p,i)
    temp[i+1,] = 1
    ek2[,(count+1):(count+i)] = temp
    count = count+i
  }
  ak1 = A %*% ek1; ak2 = A %*% ek2
  return(list(al1 = al1, al2 = al2, ak1 = ak1, ak2 = ak2,
              el1 = el1, el2 = el2, ek1 = ek1, ek2 = ek2))
}


update_A = function(X, nu1, nu2, lambda_1, lambda_2, v, z,n,p){
  En = diag(0:(n - 1)) + diag((n - 1):0) - matrix(1, n, n) + diag(1, n, n)
  Ep = diag(0:(p - 1)) + diag((p - 1):0) - matrix(1, p, p) + diag(1, p, p)
  
  M = diag(1,n,n) + nu1 * En
  
  N = nu2 * Ep
  
  alk = alk(X,n,p)
  
  el1=alk$el1
  el2=alk$el2
  ek1=alk$ek1
  ek2=alk$ek2
  
  lv = lambda_1+ nu1 * v
  lz = lambda_2 + nu2 * z
  C2 = 0 -(el2-el1) %*% t(lv)
  C3 = lz %*% t(ek1-ek2)
  C = X +  C2 + C3  
  
  A = sylvester(M,t(N),C)
  return(A)
}