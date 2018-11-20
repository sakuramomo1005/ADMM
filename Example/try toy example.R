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

tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
}

vec2tri <- function(k,n) {
  i <- ceiling(0.5*(2*n-1 - sqrt((2*n-1)^2 - 8*k)))
  j <- k - n*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
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

create_adjacency <- function(V,w,n,method='ama') {
  if (!is.null(method) && !(method %in% c("ama","admm")))
    stop("method must be 'ama', 'admm', or NULL.")
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0)
  if (method == 'ama') {
    ix <- vec2tri(which(w>0),n)
  } else {
    ix <- vec2tri(1:(n*(n-1)/2),n)
  }
  i <- ix[connected_ix,1]
  j <- ix[connected_ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}
### Prediction results for sparse convex clustering
fit_sparse <- function(X, gamma1, gamma2, Gamma2_weight, k_w, phi, method = "ama", data_valide, type = 2, dist.type = "euclidean", verbose = TRUE ){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #   dist.type <- c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")
  w <- dist_weight( t(X) / sqrt(p),phi, dist.type = dist.type, p = 2 )
  if(verbose)  hist(w, breaks = 20, main = "Histogram of weight")
  w <- knn_weights(w,k = k_w,n)
  w <- w / sum(w) * sqrt(p) * sqrt(n)   # Standardize Lambda1 weight
  
  nu <- AMA_step_size(w,n)/2  # Largest nu we can have
  
  
  fit <- scvxclustr::scvxclust(X = X, w = w, Gamma1 = gamma1, Gamma2 = gamma2, Gamma2_weight = Gamma2.weight, nu = nu, method = method, tol_abs = 1e-7)
  
  label_est <- list()
  features_est <- list()
  rand_est <- vector()
  for(i in 1:length(gamma1)){
    A <- create_adjacency( round(fit$V[[i]],2), w, n, method = method)
    label_est[[i]] <- find_clusters(A)
    features_est[[i]] <- (! apply(fit$U[[i]], 2, sd) == 0)
    rand_est[i] <- adjustedRandIndex( label_est[[i]]$cluster, data$label)
  }
  
  fit$cluster  <- lapply(label_est, function(x) x$cluster)
  fit$size     <- lapply(label_est, function(x) x$size)
  fit$rand     <- rand_est
  fit$gamma    <- gamma1
  fit$features <- features_est
  fit$X        <- X
  
  
  cluster0 = fit$cluster[[1]]
  features0 = fit$features[[1]]
  pred = predict_validate( t(X), cluster = cluster0, features = features0, data_valide)
  fit$predict = c(gamma1 = gamma1, gamma2 = gamma2, pred)
  
  fit
}

set.seed(123)

# Load scvxclustr package
library(devtools)
library(scvxclustr)


n = 60          # Sample size
true_p = 20     # Number of true features
p = 40        # Number of total features (p = 150, 500)
k = 4           # Number of true cluster
mu = 1.2        # mean of normal distribution (mu = 0.6, 0.9)
sigma = 1       # sd of normal distribution
method = "admm"  # Fitted method (method = "ama", "admm")

# Simiulate 4 cluster Gaussian data
data <- simu_4(n = n, true_p = true_p, p = p, k = k, mu = mu, sigma = sigma )

#standardize n by p data matrix
X <- scale(data$X,center=TRUE,scale=FALSE)

# Validation data
data_valide <- list()
for(i in 1:5){
  data_valide[[i]] <- simu_4(n = n, true_p = true_p, p = p, k = k, mu = mu, sigma = sigma )
}

# Adaptive Weight (if possible)
g1 <-6
g2 <- 0
Gamma2.weight <- c(rep(0.5, true_p), rep(1,p - true_p) )
k_w <- 5    # Number of nearest neighbors
phi <- 0.5  # scale of the kernel
verbose <- TRUE # show more information
w <- dist_weight((X)/sqrt(p),phi, dist.type = "euclidean", p = 2 )
w <- knn_weights(w,k = k_w,n)
nu <- AMA_step_size(w,n) /2

## Validate the cvxclust and scvxclust is the same when g2 = 0
# Fit a convex clustering model
fit1 <- cvxclust(X = t(X), w = w, gamma = g1, method = "ama", nu = nu, max_iter = 10000, tol = 1e-5)

# Fit a sparce convex clsutering model
fit3 <- scvxclust(X = X, w = w, Gamma1 = g1, Gamma2 = g2, Gamma2_weight = Gamma2.weight, method = method, nu = nu, max_iter = 1000, tol_abs = 1e-3)

## Validate the sparse convex clustring model create a correct clustering structure under the tody exmaple.
g1 <- 9
g2 <- 10
fit_predict <- fit_sparse(X = X, gamma1 = g1, gamma2 = g2, Gamma2.weight, k_w, phi, method = method, data_valide = data_valide, verbose = F)
fit_predict$predict
table(data$label, fit_predict$cluster[[1]])

ix <- 4
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

pdf('toy3.pdf')
heatmap(X,col=hmcols,labRow=NA,labCol=NA)
#heatmap(fit1$U[[1]],col=hmcols,labRow=NA,labCol=NA)
heatmap(fit1$U[[1]],col=hmcols,labRow=NA,labCol=NA)
heatmap(fit3$U[[1]],col=hmcols,labRow=NA,labCol=NA)
dev.off()