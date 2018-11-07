# 2018/10/27

# try to compare

# new ADMM 
# with sparse_cvxclust_path_admm, send code

# generate data
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
X = simu_4(40,2,10,1,1,1)
X = X$X

n=dim(X)[1]; p=dim(X)[2]
eplison_n=L_num(n)
eplison_p=L_num(p)

# calculate weigths outside
w <- dist_weight(X / sqrt(p),phi = 0.5)
w <- knn_weights(w,k = 5,n)

Gamma1 = 100
Gamma2 = 0

gamma_1 = 100; gamma_2 = 0; nu1 = 0.05; nu2 = 0.05

x1 = sparse_cvxclust_path_admm(X,w,Gamma1, Gamma2,nu=1)
res = bi_ADMM(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, niter = 1000, imax = 0.1)
res2 = t(res$A)
mean((res2 - x1$U[[1]])^2) # 3.0874098