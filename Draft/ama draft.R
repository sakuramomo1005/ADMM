### new update A 
### AMA 
### 11/05


library(MASS)
library(Matrix)

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
elk = function(n,p){
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
  return(list(el1 = el1, el2 = el2, ek1 = ek1, ek2 = ek2))
}
dist_weight <- function(X, phi){
  dist_X <- as.numeric(dist(X))
  return(exp(- phi * dist_X^2))
}
knn_weights <- function(w,k,n) {
  i <- 1
  neighbors <- tri2vec(i,(i+1):n,n)
  keep <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  for (i in 2:(n-1)) {
    group_A <- tri2vec(i,(i+1):n,n)
    group_B <- tri2vec(1:(i-1),i,n)
    neighbors <- c(group_A,group_B)
    knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep <- union(knn,keep)
  }
  i <- n
  neighbors <- tri2vec(1:(i-1),i,n)
  knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep <- union(knn,keep)
  w[-keep] <- 0
  return(w)
}
prox = function(v,sigma){
  return(max(1-sigma/sum(t(v) %*% v),0) * v)
}
tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
}
L_num = function(n){
  t = matrix(0, n*(n-1)/2,2)
  count = 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      count = count + 1
      t[count,1]= i; t[count,2]=j
    }
  }
  LL = data.frame(l1 = t[,1],l2 = t[,2])
  LL = LL[order(LL$l2),]
  rownames(LL) = NULL
  return(LL)
}
bi_ama = function(X, nu1, nu2, gamma_1, gamma_2, niter = 1000, tol = 1e-5){
  
  A = X
  n=dim(X)[1]; p=dim(X)[2]
  eplison_p=L_num(p)
  eplison_n=L_num(n)
  
  lambda_1 = matrix(1,p,dim(eplison_n)[1])
  lambda_2 = matrix(1,n,dim(eplison_p)[1])
  
  w = dist_weight(X / sqrt(p),phi = 0.5)
  w_l = knn_weights(w,k = 1,n)
  u = dist_weight(t(X) / sqrt(n),phi = 0.5)
  u_k = knn_weights(u,k = 1,p)
  
  alk = elk(n,p)
  el1=alk$el1
  el2=alk$el2
  ek1=alk$ek1
  ek2=alk$ek2
  
  for(iter in 1:niter){
    print('iter')
    print(iter)
    A_old = A; lambda_1_old = lambda_1; lambda_2_old = lambda_2
    C2 = 0 -(el2-el1) %*% t(lambda_1)
    C3 = lambda_2 %*% t(ek1-ek2)
    A = X +  C2 + C3  
    
    al1 = t(A) %*% el1; al2 = t(A) %*% el2
    ak1 = A %*% ek1; ak2 = A %*% ek2
    
    lambda1 = lambda_1 - nu1 * (al1 - al2)
    lambda2 = lambda_2 - nu2 * (ak1 - ak2)
    
    for(i in 1:dim(lambda_1)[2]){
      tau = gamma_1 * w_l[i]
      norm_x = sum(lambda1[,i]^2)
      if(norm_x > tau){
        lambda_1[,i] = tau/norm_x * lambda_1[,i]
      }else{
        lambda_1[,i] = lambda_1[,i]
      }
    }
    
    for(i in 1:dim(lambda_2)[2]){
      tau = gamma_2 * u_k[i]
      norm_y = sum(lambda2[,i]^2)
      if(norm_y > tau){
        lambda_2[,i] = tau/norm_y * lambda_2[,i]
      }else{
        lambda_2[,i] = lambda_2[,i]
      }
    }
    
    print(paste('A',sum(abs(A - A_old))))
    print(paste('l',sum(abs(lambda_1 - lambda_1_old))))
    print(paste('2',sum(abs(lambda_2 - lambda_2_old))))
    
    if(sum(abs(A - A_old)) < tol & 
       sum(abs(lambda_1 - lambda_1_old)) < tol & 
       sum(abs(lambda_2 - lambda_2_old)) < tol){
      return(list(A = A, 
                  lambad_1 = lambda1, lambad_2 = lambda2, niter = iter))
      break
    }
    
  }
  
  if(iter == niter){
    print(paste('not converge within',iter, 'times'))
    return(list(A = A, 
                lambad_1 = lambda1, lambad_2 = lambda2, niter = iter))
  }
}

X = matrix(1:18,6,3)
nu1 = nu2 = 1
gamma_1 = 6 
gamma_2 = 0.1

bi_ama(X, nu1, nu2, gamma_1, gamma_2)
