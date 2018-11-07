### 2018/10/18

### first version, demo
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
    A = update_A(X, nu1, nu2, lambda_1, lambda_2, v, z)
    
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

# test 1
simu_4 = function(n, true_p, p, k, mu, sigma, seed = NULL){
  
  if(! is.null(seed) ){ set.seed(seed)  }
  n_p <- true_p / 2
  n_k <- n / k
  clust.ind <- rep(1:k, each = n_k)
  clust.mat <- rbind( c(rep( mu, n_p), rep(-mu, n_p)), 
                      c(rep(-mu, n_p), rep(-mu, n_p)),
                      c(rep(-mu, n_p), rep( mu, n_p)),
                      c(rep( mu, n_p), rep( mu, n_p))
  )
  
  X = matrix(0,n,p)
  for(i in 1:n){
    mu_mean <- c( clust.mat[clust.ind[i],], rep(0, p - true_p) )
    X[i,] <- rnorm(p, mu_mean,rep(sigma,p))
  }
  
  list(X = X, label = clust.ind, features = c(rep(TRUE, true_p), rep(FALSE, p - true_p)))
  
}
X =simu_4(40,2,10,1,1,1)
X = X$X
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

set.seed(123)
lambda_1 = matrix(rnorm(p*dim(eplison_n)[1],0,1),p,dim(eplison_n)[1])
v = matrix(rnorm(p*dim(eplison_n)[1],1,3),p,dim(eplison_n)[1])
lambda_2 = matrix(rnorm(n*dim(eplison_p)[1],0,1),n,dim(eplison_p)[1])
z = matrix(rnorm(n*dim(eplison_p)[1],-1,2),n,dim(eplison_p)[1])

gamma_1 = 100; gamma_2 = 0; nu1 = 0.05; nu2 = 0.05

res = bi_ADMM(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, niter = 1000, imax = 0.1)

res = bi_ADMM(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, niter = 1000, imax = 0.1)




# test 2

X =simu_4(400,2,100,1,1,1)
X = X$X
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

set.seed(123)
lambda_1 = matrix(rnorm(p*dim(eplison_n)[1],0,1),p,dim(eplison_n)[1])
v = matrix(rnorm(p*dim(eplison_n)[1],1,3),p,dim(eplison_n)[1])
lambda_2 = matrix(rnorm(n*dim(eplison_p)[1],0,1),n,dim(eplison_p)[1])
z = matrix(rnorm(n*dim(eplison_p)[1],-1,2),n,dim(eplison_p)[1])

lambda_1[1:10]
gamma_1 = 100; gamma_2 = 0; nu1 = 0.05; nu2 = 0.05
begin = Sys.time()
res = bi_ADMM(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, niter = 1000, imax = 0.1)
end = Sys.time()
end - begin