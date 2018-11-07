library(cvxclustr)

p <- 10
n <- 20
seed <- 12345
nProbs <- 10
errors <- double(nProbs)

rnd_problem <- create_clustering_problem(p,n,seed=seed,method='admm')
X <- rnd_problem$X
ix <- rnd_problem$ix
M1 <- rnd_problem$M1
M2 <- rnd_problem$M2
s1 <- rnd_problem$s1
s2 <- rnd_problem$s2
w <- rnd_problem$w
nK <- length(w)
Lambda <- matrix(rnorm(p*nK),p,nK)
gamma <- 0.1
nu <- 1
max_iter <- 1e6
tol_abs <- 1e-15
tol_rel <- 1e-15
sol_admm_acc <- cvxclust_admm(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,max_iter=max_iter,
                              tol_abs=tol_abs,tol_rel=tol_rel,accelerate=TRUE)







bi_ADMM=function(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, niter = 1000, imax = 0.1){
  A = X
  n = dim(X)[1]; p = dim(X)[2]
  eplison_p = L_num(p); eplison_n = L_num(n)
  fx=rep(0,niter)
  A_old = A; v_old = v; z_old = z; 
  lambda_1_old = lambda_1; lambda_2_old = lambda_2
  for(iter in 2:niter){
    if(iter %% 10 ==0){print(paste('iter',iter))}
    # update A 
    begin = Sys.time()
    A = update_A(X, nu1, nu2, lambda_1, lambda_2, v, z)
    end = Sys.time()
    end - begin
    # update V and Z
    vz = update_vz(X, A, lambda_1, lambda_2, gamma_1, gamma_2, nu1, nu2)
    v = vz[[1]]
    z = vz[[2]]
    
    # update lambda_1 and lambda_2 
    lambda = update_lambda(X, A, nu1, nu2,v, z)
    lambda1 = lambda[[1]]
    lambda2 = lambda[[2]]
    
    print(paste('A:',round(sum(abs(A - A_old)),2),
                '  v:',round(sum(abs(v - v_old)),2) ,
                '  z:',round(sum(abs(z - z_old)),2) ,
                '  l1:',round(sum(abs(lambda1 - lambda_1_old)),2),
                '  l2:',round(sum(abs(lambda2 - lambda_2_old)),2)))
    
    term1=term2=term3=0
    for(i in 1:n){
      term1 = term1 + t(X[i,] - A[i,]) %*% (X[i,] - A[i,])
    }
    
    for(i in  1:dim(eplison_n)[1]){
      l1=eplison_n[i,'l1']
      l2=eplison_n[i,'l2']
      term2 = term2 + t(v[,i] + 1/nu1 * lambda1[,i]- A[l1,] + A[l2,]) %*% (v[,i] - A[l1,] + A[l2,])
    }
    
    for(i in  1:dim(eplison_p)[1]){
      l1=eplison_p[i,'l1']
      l2=eplison_p[i,'l2']
      term3 = term3 + t(z[,i] + 1/nu2 * lambda2[,i]- A[,l1] + A[,l2]) %*% (z[,i] - A[,l1] + A[,l2])
    }
    
    fx[iter] = 0.5 * term1 + 0.5* nu1 * term2 + 0.5*nu2*term3
    
    if(fx[iter] < fx[iter - 1]){print('small')}
    
    # whether coverage
    if(sum(abs(A - A_old)) < e & 
       sum(abs(v - v_old)) < e & 
       sum(abs(z - z_old)) < e & 
       sum(abs(lambda1 - lambda_1_old)) < e & 
       sum(abs(lambda2 - lambda_2_old)) < e){
      return(list(A = A, v = v, z = z, lambad_1 = lambda1, lambad_2 = lambda2, niter = iter-1))
      break
    }
    A_old = A; v_old = v; z_old = z; 
    lambda_1_old = lambda1; lambda_2_old = lambda2
  }
  if(iter == 1000){
    print(paste('not converge within',iter, 'times'))
    return(list(A = A, v = v, z = z, lambad_1 = lambda1, lambad_2 = lambda2, niter = iter-1))
  }
}




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


# optimization
# change L_num:
L_num = function(n){
  t = matrix(0, n*(n-1)/2,2)
  count = 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      count = count + 1
      t[count,1]= i; t[count,2]=j
    }
  }
  return(LL = data.frame(l1 = t[,1],l2 = t[,2]))
}

### 2.  sylvester equation
sylvester2 = function(A, B, C, eps = 0.0001){
  
  library(MASS)
  library(Matrix)
  
  t1 = Sys.time()
  A1 = Schur(A)
  Q1 = A1$Q; R1 = A1$T
  
  A2 = Schur(B)
  Q2 = A2$Q; R2 = A2$T 
  C = t(Q1) %*% C %*% Q2
  
  Rsq = R1 * R1
  I = diag(dim(A)[1])
  
  k = 1
  n = dim(R2)[1]
  t2 = Sys.time()
  while(k < n + 1){
    if(k < n){
      if(abs(R2[k+1, k]) < eps){
        left = R1 + R2[k,k] * I
        right = C[,k]
        temp = matrix(0, dim(X)[1],1)
        if(k == 1){
          temp = temp
        }else{
          temp = (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
        }
        
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
            temps = X[,1:(k-1)] %*% matrix(R2[1:(k-1),k:(k+1)],k-1,2)
            temp2 = temps[,1]
            temp3 = temps[,2]
        }
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
          temp = (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
        }
        X[,k] = ginv(left) %*% (right - temp)
        k = k+1
      }else{
        R22 = R2
        R22 = cbind(R2, rep(0,dim(R2)[1]))
        R22 = rbind(R22,rep(0,dim(R2)[1]+1))
        r11 = R22[k,k]; r12 = R22[k, k+1]; r21 = R22[k+1, k]; r22 = R22[k+1, k+1]
        if(k == 1){
          temp2 = temp2
          temp3 = temp3
        }else{
          temps = X[,1:(k-1)] %*% matrix(R2[1:(k-1),k:(k+1)],k-1,2)
          temp2 = temps[,1]
          temp3 = temps[,2]
        }
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

A = matrix(rnorm(3600,0,6),60,60)
B = matrix(rnorm(360*360,0,6),360,360)
X = matrix(rnorm(360*60,0,7),60,360)
C = A%*% X + X%*% B
XX = X
t1 =Sys.time()
res = sylvester2(A,B,C)
t2 =Sys.time()
t2-t1
### 3. update A
update_A = function(X, nu1, nu2, lambda_1, lambda_2, v, z){
  library(MASS)
  library(Matrix)
  n=dim(X)[1]; p=dim(X)[2]
 
  En = diag(0:(n - 1)) + diag((n - 1):0) - matrix(1, n, n) + diag(1, n, n)
  Ep = diag(0:(p - 1)) + diag((p - 1):0) - matrix(1, p, p) + diag(1, p, p)
  
  M = diag(1,n,n) + nu1 * En
  
  N = nu2 * Ep
  
  # for C

  
  eplison_n=L_num(n)
  eplison_p=L_num(p)
  t1 = Sys.time()
  c1 = c2 = matrix(0,n,p)
  
  tmp = lambda_1
  tmp2 = nu1 * v
  
  tmp = data.frame(tmp)
  tmp2 = data.frame(tmp2)
  tmp1 = data.frame(tmp)
  tmp22 = data.frame(tmp2)
  
  colnames(tmp) = eplison_n[,1]
  colnames(tmp2) = eplison_n[,1]
  colnames(tmp1) = eplison_n[,2]
  colnames(tmp22) = eplison_n[,2]
  
  for(i in 1:(n-2)){
    j = n+1-i
    c1[i,] = apply(tmp[,colnames(tmp) == i], 1, sum) +
             apply(tmp2[,colnames(tmp2) == i], 1, sum)
    c2[j,] = apply(tmp[,colnames(tmp1) == j], 1, sum)  +
             apply(tmp2[,colnames(tmp22) == j], 1, sum)
  }
  c1[n-1,] = tmp[,colnames(tmp) == i+1] + tmp2[,colnames(tmp2) == i+1]
  c2[2,] = tmp[,1] + tmp2[,1]
  c = c1 - c2 
  t2 = Sys.time()
  
  
  
  C2 = matrix(0,n,p) 
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    el1=matrix(rep(0,n),n,1)
    el2=matrix(rep(0,n),n,1)
    el1[l1,1]=1;el2[l2,1]=1
    C2 = C2 + (el1) %*% t(matrix(1:25,5,10)[,i])
    print(C2)
  }
  
  t2 = Sys.time()
  n=6
  eplison_n = L_num(n)
  nn = dim(eplison_n)[1]
  tt = matrix(0,n,nn)
  colnames(tt) = eplison_n[,1]
  for(i in 1:(n-1)){
    tt[i,which(colnames(tt)==i)]=1
  }
  
  c1=tt %*% t(lambda_1+ nu1 * v)
  
  tt = matrix(0,n,nn)
  colnames(tt) = eplison_n[,2]
  for(i in 1:n){
    print(i)
    tt[i,which(colnames(tt)==i)]=1
  }
  
  c2=tt %*% t(lambda_1+ nu1 * v)
  
  c = c1-c2
  
  t1 = Sys.time()
  
  
  
  C2 = matrix(0,n,p) 
  t = matrix(1:25,5,dim(eplison_n)[1])
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    el1=matrix(rep(0,n),n,1)
    el2=matrix(rep(0,n),n,1)
    el1[l1,1]=1;el2[l2,1]=1
    C2 = C2 + (-el2) %*% t(matrix(1:25,5,10)[,i])
    print(C2)
  }
  
  
  
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


