# 1106

# try to compare 
# admm
# ama

# different methods from different people


set.seed(123)

# Load scvxclustr package
library(scvxclustr)
library(cvxclustr)
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

### Prediction results for convex clustering
fit_cvxclust <- function(X, gamma1, k_w, phi, method, data_valide, type = 2, dist.type = "euclidean", verbose = TRUE ){
  fit <- list()
  n <- ncol(X)
  p <- nrow(X)
  
  #   dist.type <- c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")
  w <- dist_weight(X / sqrt(p),phi, dist.type = dist.type, p = 2 )
  if(verbose)  hist(w, breaks = 20, main = "Histogram of weight")
  w <- knn_weights(w,k = k_w,n)
  
  nu <- AMA_step_size(w,n)/2  # Largest nu we can have
  
  
  fit <- cvxclust(X = X, w = w, gamma = gamma1, method = method, nu = nu, type = type)
  
  label_est <- list()
  features_est <- list()
  rand_est <- vector()
  for(i in 1:length(gamma1)){
    A <- create_adjacency( round(fit$V[[i]],2), w, n, method = method)
    label_est[[i]] <- find_clusters(A)
    features_est[[i]] <- (! apply(fit$U[[i]], 1, sd) == 0)
    rand_est[i] <- adjustedRandIndex( label_est[[i]]$cluster, data$label)
  }
  
  fit$cluster  <- lapply(label_est, function(x) x$cluster)
  fit$size     <- lapply(label_est, function(x) x$size)
  fit$rand     <- rand_est
  fit$gamma    <- gamma1
  fit$features <- rep(TRUE, p)
  fit$X        <- X
  
  fit$predict <- list()
  for(iter in 1:length(gamma1) ){
    cluster0 = fit$cluster[[iter]]
    features0 = fit$features
    pred = predict_validate(X, cluster = cluster0, features = features0, data_valide = data_valide)
    fit$predict[[iter]] <- c(gamma1 = gamma1[iter], gamma2 = NA, pred)
  }
  fit$predict <- do.call(rbind, fit$predict)
  
  fit
}
n = 100         # Sample size
true_p = 20     # Number of true features
p = 50     # Number of total features (p = 150, 500)
k = 4           # Number of true cluster
mu = 1.2        # mean of normal distribution (mu = 0.6, 0.9)
sigma = 1       # sd of normal distribution
method = "admm"  # Fitted method (method = "ama", "admm")

# Simiulate 4 cluster Gaussian data
data <- simu_4(n = n, true_p = true_p, p = p, k = k, mu = mu, sigma = sigma )

#standardize n by p data matrix
X <- scale(data$X,center=TRUE,scale=FALSE)

# Adaptive Weight (if possible)
g1 <-6
g2 <- 0
Gamma2.weight <- c(rep(0.5, true_p), rep(1,p - true_p) )
k_w <- 5    # Number of nearest neighbors
phi <- 0.5  # scale of the kernel
verbose <- TRUE # show more information
w <- dist_weight( t(X) / sqrt(p),  p = 2 )
w <- knn_weights(w,k = k_w,n)
nu <- AMA_step_size(w,n) /2


## Validate the cvxclust and scvxclust is the same when g2 = 0
# Fit a convex clustering model
x1 = Sys.time()
fit1 <- cvxclust(X = t(X), w = w, gamma = g1, method = "admm", nu = nu, max_iter = 10, tol = 1e-5)
x2 = Sys.time()
# Fit a sparce convex clsutering model
x3 = Sys.time()
fit2 <- scvxclust(X = X, w = w, Gamma1 = g1, Gamma2 = g2, Gamma2_weight = Gamma2.weight, method = method, nu = nu, max_iter = 10, tol_abs = 1e-5)
x4 = Sys.time()
fit1$iters
fit2$iters
diff_U <- as.numeric(fit1$U[[1]] - t(fit2$U[[1]]) )
summary( diff_U )
plot(diff_U)


library(cvxbiclustr)
#X <- lung
X <- X - mean(X)
X <- X/norm(X,'f')
## Create annotation for heatmap
types <- colnames(lung)
ty <- as.numeric(factor(types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)
## Construct weights and edge-incidence matrices
phi <- 0.5; k <- 5
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col
## Connected Components of Row and Column Graphs
wts$nRowComp
wts$nColComp
#### Initialize path parameters and structures
nGamma <- 5
gammaSeq <- 10**seq(0,3,length.out=nGamma)
## Generate solution path
sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq)
ix <- 4
#png('eric.png')
#heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])
#dev.off()
heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA)

rand.index(X,sol$U[[ix]])

ty <- as.numeric(factor(1:20))
cols <- rainbow(4)

heatmap(t(fit1$U[[1]]),col=hmcols,labRow=NA,labCol=NA)
heatmap((fit2$U[[1]]),col=hmcols,labRow=NA,labCol=NA)
png('new function.png')
heatmap(fit3$A,col=hmcols,labRow=NA,labCol=NA,main='new function')
dev.off()

types <- colnames(fit3$X)
ty <- as.numeric(factor(1:20))
cols <- rainbow(4)
heatmap(fit3$A,col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])




#### lalala let's do simulation!

n = 200
p = 200
mu = seq(-6,6,0.5)
theta = 1.5
x = matrix(0, n, p)

# 2 row groups
# 8 column groups
row_group = 2
col_group = 8

mu_rc = matrix(0, row_group,col_group)
for(r in 1:row_group){
  for(c in 1:col_group){
    mu_rc[r,c] = sample(mu)[1]
  }
}

row_assign = c()
col_assign = c()
for(i in 1:n){
  row_assign = c(row_assign, sample(1:2)[1])
}
for(i in 1:p){
  col_assign = c(col_assign, sample(1:8)[1])
}

for(i in 1:n){
  for(j in 1:p){
    r = row_assign[i]
    c = col_assign[j]
    mu_temp = mu_rc[r,c]
    x[i,j] = rnorm(1,mu_temp,theta)
  }
}

rownames(x) = paste('row',row_assign,sep='')
colnames(x) = paste('col',col_assign,sep='')

heatmap(x,col=hmcols,labRow=NA,labCol=NA)

data_gen=function(seed){
  set.seed(seed)
  n = 200
  p = 200
  mu = seq(-6,6,0.5)
  theta = 1.5
  x = matrix(0, n, p)
  
  # 2 row groups
  # 8 column groups
  row_group = 2
  col_group = 8
  
  mu_rc = matrix(0, row_group,col_group)
  for(r in 1:row_group){
    for(c in 1:col_group){
      mu_rc[r,c] = sample(mu)[1]
    }
  }
  
  row_assign = c()
  col_assign = c()
  for(i in 1:n){
    row_assign = c(row_assign, sample(1:2)[1])
  }
  for(i in 1:p){
    col_assign = c(col_assign, sample(1:8)[1])
  }
  
  for(i in 1:n){
    for(j in 1:p){
      r = row_assign[i]
      c = col_assign[j]
      mu_temp = mu_rc[r,c]
      x[i,j] = rnorm(1,mu_temp,theta)
    }
  }

  rownames(x) = paste('row',row_assign,sep='')
  colnames(x) = paste('col',col_assign,sep='')
  
  return(x)
}

## Create annotation for heatmap
types <- colnames(x)
ty <- as.numeric(factor(types))
cols <- rainbow(8)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)
## Construct weights and edge-incidence matrices

eric_method = function(x){
  X <- x - mean(x)
  X <- X/norm(X,'f')
  phi <- 0.5; k <- 5
  wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
  w_row <- wts$w_row
  w_col <- wts$w_col
  E_row <- wts$E_row
  E_col <- wts$E_col
  ## Connected Components of Row and Column Graphs
  wts$nRowComp
  wts$nColComp
  #### Initialize path parameters and structures
  nGamma <- 5
  gammaSeq <- 10**seq(0,3,length.out=nGamma)
  ## Generate solution path
  sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq)
  return(sol)
}

X <- x - mean(x)
X <- X/norm(X,'f')
## Create annotation for heatmap
types <- colnames(x)
ty <- as.numeric(factor(types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)
## Construct weights and edge-incidence matrices
phi <- 0.5; k <- 5
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col
## Connected Components of Row and Column Graphs
wts$nRowComp
wts$nColComp
#### Initialize path parameters and structures
nGamma <- 5
gammaSeq <- 10**seq(0,3,length.out=nGamma)
## Generate solution path
sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq)
ix <- 4
heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])

#heatmap(fit4$A,col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])

res2 = bi_ama(X, 1,1,5,1)

head(round(abs(sol$U[[ix]] - res2$A)/sol$U[[ix]],3))
heatmap(res2$A,col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])

head(sol$U[[ix]])


rand_res = function(x,y){
  rand = c()
  for(i in 1:4){
    for(j in 1:4){
      rand = c(rand,
               rand.index(X[(1+50*(i-1)):(50*i),(1+50*(j-1)):(50*j)],
                          y[(1+50*(i-1)):(50*i),(1+50*(j-1)):(50*j)]))
    }
  }
  return(mean(rand))
}
rand = c()
for(i in 1:4){
  for(j in 1:4){
    rand = c(rand,
             rand.index(X[(1+50*(i-1)):(50*i),(1+50*(j-1)):(50*j)],
                        sol$U[[ix]][(1+50*(i-1)):(50*i),(1+50*(j-1)):(50*j)]))
  }
}



mean(rand)
library(fossil)
rand.index(X[1:50,1:50],sol$U[[ix]][1:50,1:50])
rand.index(X[1:50,1:50],res2$A[1:50,1:50])


fit4 = bi_ADMM1(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2)

a = Sys.time()
bi_ADMM1(X, nu1, nu2, lambda_1, lambda_2, v, z, gamma_1, gamma_2, niter = 1, imax = 0.1)
b = Sys.time()
heatmap(fit5$A,col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])
