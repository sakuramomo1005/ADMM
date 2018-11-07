### week 10 19

# the work admm parts functions 

# no optimization


### 1. L_num
#The function to generate $l = (l_1, l_2)$, where $l_1 < l_2$
tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
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

### 2. sylvester equation
sylvester = function(A, B, C, eps = 0.0001){
  
  library(MASS)
  library(Matrix)
  
  A1 = Schur(A)
  Q1 = A1$Q; R1 = A1$T
  
  A2 = Schur(B)
  Q2 = A2$Q; R2 = A2$T 
  C = t(Q1) %*% C %*% Q2
  
  Rsq = R1 * R1
  I = diag(dim(A)[1])
  
  k = 1
  n = dim(R2)[1]
  
  while(k < n + 1){
    if(k < n){
      if(abs(R2[k+1, k]) < eps){
        left = R1 + R2[k,k] * I
        right = C[,k]
        temp = matrix(0, dim(X)[1],1)
        if(k == 1){
          temp = temp
        }else{
          for(i in 1:(k-1)){
            temp = temp + X[,i] * R2[i,k]
          }
        }
        temp = matrix(temp, dim(C)[1],1)
        X[,k] = ginv(left) %*% (right - temp)
       # mytry = myTryCatch(solve(left))
        #if(is.null(mytry$error) == 0){er = c(er,tt)}
        k = k+1
      }else{
        r11 = R2[k,k]; r12 = R2[k, k+1]; r21 = R2[k+1, k]; r22 = R2[k+1, k+1]
        temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        if(k == 1){
          temp2 = temp2
          temp3 = temp3
        }else{
          for(i in 1:(k-1)){
            temp2 = temp2 + X[,i] * R2[i,k]
            temp3 = temp3 + X[,i] * R2[i,k+1]}
        }
        temp2 = matrix(temp2, dim(X)[1],1)
        temp3 = matrix(temp3, dim(X)[1],1)
        b1 = C[,k] - temp2 
        b2 = C[,k+1] - temp3
        b1_prime = R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime = R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime = matrix(0, dim(X)[1],2)
        b_prime[,1] = b1_prime; b_prime[,2] = b2_prime
        X[,k:(k+1)] = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                             (r11*r22 - r12*r21) * I) %*% b_prime 
        k = k+2
      }
    }else{
      if(abs(R2[1, k]) > eps){
        left = R1 + R2[k,k] * I
        right = C[,k]
        temp = matrix(0, dim(X)[1],1)
        if(k == 1){
          temp = temp
        }else{
          for(i in 1:(k-1)){
            temp = temp + X[,i] * R2[i,k]
          }
        }
        temp = matrix(temp, dim(C)[1],1)
        X[,k] = ginv(left) %*% (right - temp)
        k = k+1
      }else{
        R22 = R2
        R22 = cbind(R2, rep(0,dim(R2)[1]))
        R22 = rbind(R22,rep(0,dim(R2)[1]+1))
        r11 = R22[k,k]; r12 = R22[k, k+1]; r21 = R22[k+1, k]; r22 = R22[k+1, k+1]
        temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        for(i in 1:(k-1)){
          temp2 = temp2 + X[,i] * R22[i,k]
          temp3 = temp3 + X[,i] * R22[i,k+1]}
        temp2 = matrix(temp2, dim(X)[1],1)
        temp3 = matrix(temp3, dim(X)[1],1)
        b1 = C[,k] - temp2 
        b2 = - temp3
        b1_prime = R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime = R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime = matrix(0, dim(X)[1],2)
        b_prime[,1] = b1_prime; b_prime[,2] = b2_prime
        GOD = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                     (r11*r22 - r12*r21) * I) %*% b_prime 
        X[,k] = GOD[,1]
        k = k+2
      }
    }
  }
  return(Q1 %*% X %*% t(Q2))
}


### 3. update A
update_A = function(X, nu1, nu2, lambda_1, lambda_2, v, z){
  library(MASS)
  library(Matrix)
  n=dim(X)[1]; p=dim(X)[2]
  eplison_n=L_num(n)
  eplison_p=L_num(p)
  
  En=matrix(0,n,n)
  Ep=matrix(0,p,p)
  
  # calculate \sum_{l \in \epsilon_1}(e_{l1} - e_{l_2})(e_{l1} - e_{l_2})^T, which is En
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    el1=matrix(rep(0,n),n,1)
    el2=matrix(rep(0,n),n,1)
    el1[l1,1]=1;el2[l2,1]=1
    En=En+(el1-el2) %*% t(el1-el2)
  }
  
  # calculate \sum_{k \in \epsilon_2} (e_{k1} - e_{k_2})(e_{k1} - e_{k_2})^T, which is Ep
  for(i in 1:dim(eplison_p)[1]){
    l1=eplison_p[i,'l1']
    l2=eplison_p[i,'l2']
    el1=matrix(rep(0,p),p,1)
    el2=matrix(rep(0,p),p,1)
    el1[l1,1]=1;el2[l2,1]=1
    Ep=Ep+(el1-el2) %*% t(el1-el2)
  }
  
  M = diag(1,n,n) + nu1 * En
  
  N = nu2 * Ep
  t3 = Sys.time()
  # C2 is the second part of C, which is \sum_{l \in \epsilon_1} (e_{l1} - e_{l_2})(\lambda_{l1} + \nu_1 v_l)^T 
  C2 = matrix(0,n,p) 
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    el1=matrix(rep(0,n),n,1)
    el2=matrix(rep(0,n),n,1)
    el1[l1,1]=1;el2[l2,1]=1
    C2 = C2 + (el1-el2) %*% t(lambda_1[,i]+ nu1 * v[,i])
  }
  t4= Sys.time()
  # C3 is the third part of C, which is \sum_{k \in \epsilon_2} (\lambda_{k2} + \nu_2 z_k)(e_{k1} - e_{k_2})^T
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
  
  A = sylvester(M,t(N),C)
  
  return(A)
}

# test the function 
X = matrix(rnorm(100),10,10)
n=dim(X)[1]; p=dim(X)[2]
eplison_p=L_num(p)
eplison_n=L_num(n)
nu1=0.001; nu2=0.001
lambda_1 = matrix(rnorm(p*dim(eplison_n)[1],0,1),p,dim(eplison_n)[1])
v = matrix(rnorm(p*dim(eplison_n)[1],1,3),p,dim(eplison_n)[1])
lambda_2 = matrix(rnorm(n*dim(eplison_p)[1],0,1),n,dim(eplison_p)[1])
z = matrix(n*dim(eplison_p)[1],rnorm(-1,2),n,dim(eplison_p)[1])

A=update_A1(X, nu1, nu2, lambda_1, lambda_2, v, z)
A[1:10,1:10]

### 4. update V and Z
#prox = function(v,sigma){
#  return(max(1-sigma/sum(v %*% t(v)),0) %*% v)
#}

prox = function(v,sigma){
  return(max(1-sigma/sum(t(v) %*% v),0) * v)
}
v = 1:5
prox(v, 0.1)

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

update_vz = function(X, A, lambda_1, lambda_2, gamma_1, gamma_2, nu1, nu2){
  
  n=dim(X)[1]; p=dim(X)[2]
  eplison_n=L_num(n)
  eplison_p=L_num(p)
 
  w <- dist_weight(X / sqrt(p),phi = 0.5)
  w_l <- knn_weights(w,k = 1,n)

  u = dist_weight(t(X) / sqrt(n),phi = 0.5)
  u_k <- knn_weights(u,k = 1,p)
  
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    a_l1 = A[l1,]; a_l2 = A[l2,]
    v_temp = a_l1 - a_l2 - 1/nu1 * lambda_1[,i]
    sigma_1l = gamma_1 * w_l[i]/nu1
    v[,i] = prox(v_temp,sigma_1l)
  }
  
  for(i in 1:dim(eplison_p)[1]){
    l1=eplison_p[i,'l1']
    l2=eplison_p[i,'l2']
    a_l1 = A[,l1]; a_l2 = A[,l2]
    v_temp = a_l1 - a_l2 - 1/nu2 * lambda_2[,i]
    #u_k = exp(-0.5 * (t(X[,l1] - X[,l2]) %*% (X[,l1] - X[,l2])))
    sigma_2k = gamma_2 * u_k[i]/nu2
    z[,i] = prox(v_temp,sigma_2k)
  }
  
  return(list(v = v, z = z))
}


### 5. update lambda
update_lambda=function(X, A, nu1, nu2,v, z){
  n = dim(X)[1]; p = dim(X)[2]
  eplison_p = L_num(p)
  eplison_n = L_num(n)
  # update lambda 1
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    a_l1 = matrix(A[l1,],p,1)
    a_l2 = matrix(A[l2,],p,1)
    lambda_1[,i] = lambda_1[,i] + nu1 * (v[,i] - a_l1 + a_l2)
  } 
  # update lambda 2
  for(i in 1:dim(eplison_p)[1]){
    l1=eplison_p[i,'l1']
    l2=eplison_p[i,'l2']
    a_k1 = matrix(A[,l1],n,1)
    a_k2 = matrix(A[,l2],n,1)
    lambda_2[,i] = lambda_2[,i] + nu2 * (z[,i] - a_k1 + a_k2)
  }
  return(lambda=list(lambda_1,lambda_2))
}

### ADMM function
# Combine previous function together to correct ADMM function

Bi_ADMM = function(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, e = 0.001){
  A = X; 
  for(iter in 1: 100){
    A_old = A; v_old = v; z_old = z; lambda_1_old = lambda_1; lambda_2_old = lambda_2
    
    # update A 
    A = update_A(X, nu1, nu2, lambda_1, lambda_2, v, z)
    
    # update V and Z
    vz = update_vz(X, A, lambda_1, lambda_2, gamma_1, gamma_2, nu1, nu2)
    v = vz[[1]]
    z = vz[[2]]
    
    # update lambda_1 and lambda_2 
    lambda = update_lambda(X, A, nu1, nu2,v, z)
    lambda_1 = lambda[[1]]
    lambda_2 = lambda[[2]]
    print(paste('A',sum(abs(A - A_old))))
    print(paste('v',sum(abs(v - v_old))))
    print(paste('z',sum(abs(z -z_old))))
    print(paste('l',sum(abs(lambda_1 - lambda_1_old))))
    print(paste('2',sum(abs(lambda_2 - lambda_2_old))))
    
    # whether coverage
    if(sum(abs(A - A_old)) < e & 
       sum(abs(v - v_old)) < e & 
       sum(abs(z - z_old)) < e & 
       sum(abs(lambda_1 - lambda_1_old)) < e & 
       sum(abs(lambda_2 - lambda_2_old)) < e){
      return(list(A = A))
      break
    }
  }
  if(iter == 100){print('not converage within 100 iters')
    return(list(A = A, v = v, z = z, 
                lambda_1 = lambda_1, lambda_2 = lambda_2))}
}

# test
nu1 = 1; nu2 = 1; 

X = cbind(rnorm(10,3,1), rnorm(10,-3,1),rep(0,20),rep(0,20),rep(0,20),rep(0,20))
X = rbind(X, cbind(rnorm(10,6,1), rnorm(10,-6,1),rep(0,20),rep(0,20),rep(0,20),rep(0,20)))
A = X
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])

gamma_1 = 0.5; gamma_2 = 0; nu1 = nu2 = 0.05

B = Bi_ADMM(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2)
A[1:10,]




Bi_ADMM = function(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, e = 0.001){
  A = X; 
  n = dim(X)[1]; p = dim(X)[2]
  eplison_p = L_num(p)
  eplison_n = L_num(n)
  A_old = matrix(111,dim(A)[1],dim(A)[2])
  v_old = matrix(111,p,dim(eplison_n)[1])
  z_old = matrix(111,n,dim(eplison_p)[1])
  lambda_1_old = matrix(111,p,dim(eplison_n)[1])
  lambda_2_old = matrix(111,n,dim(eplison_p)[1])
  e = 0.1
  for(iter in 1: 100){
    print('iter')
    print(iter)
    
    print(paste('old: '))
    print(lambda_1_old[1:5])
    
    # update A 
    A = update_A(X, nu1, nu2, lambda_1, lambda_2, v, z)
    
    # update V and Z
    vz = update_vz(X, A, lambda_1, lambda_2, gamma_1, gamma_2, nu1, nu2)
    v = vz[[1]]
    z = vz[[2]]
    
    # update lambda_1 and lambda_2 
    lambda = update_lambda(X, A, nu1, nu2,v, z)
    lambda_1 = lambda[[1]]
    lambda_2 = lambda[[2]]
    print(paste('lambda: '))
    print(lambda_1[1:5])
 
    # whether coverage
    if(sum(abs(A - A_old)) < e & 
       sum(abs(v - v_old)) < e & 
       sum(abs(z - z_old)) < e & 
       sum(abs(lambda_1 - lambda_1_old)) < e & 
       sum(abs(lambda_2 - lambda_2_old)) < e){
      return(list(A = A))
      break
    }
    print(paste('A',sum(abs(A - A_old))))
    print(paste('v',sum(abs(v - v_old))))
    print(paste('z',sum(abs(z -z_old))))
    print(paste('l',sum(abs(lambda_1 - lambda_1_old))))
    print(paste('2',sum(abs(lambda_2 - lambda_2_old))))
    A_old = A; v_old = v; z_old = z; lambda_1_old = lambda_1; lambda_2_old = lambda_2
    
   
  }
  if(iter == 100){print('not converage within 100 iters')
    return(list(A = A, v = v, z = z, 
                lambda_1 = lambda_1, lambda_2 = lambda_2))}
}
