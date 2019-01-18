### 2019/01/13

#similar simulation with 2018/12 
# the simulation is based on the idea from Eric
# row color bar added

# load functions
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 1115 codes')
source('admm.R')
source('ama.R')

library(fossil)
library(cvxbiclustr)
library(cvxclustr)

# data generation
# data generation
data_gen=function(seed){
  set.seed(seed)
  n = 80
  p = 40
  mu = seq(-10,10,10)
  theta = 5
  x = matrix(0, n, p)
  
  # 4 row groups
  # 4 column groups
  row_group = 4
  col_group = 4
  
  mu_rc = matrix(0, row_group,col_group)
  for(r in 1:row_group){
    for(c in 1:col_group){
      mu_rc[r,c] = sample(mu)[1]
    }
  }
  
  row_assign = c()
  col_assign = c()
  for(i in 1:n){
    row_assign = c(row_assign, sample(1:4)[1])
  }
  for(i in 1:p){
    col_assign = c(col_assign, sample(1:4)[1])
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

X <- data_gen(123)
X <- X - mean(X)
X <- X/norm(X,'f')

# draw the heatmap for the raw data
col_types = colnames(X)
col_ty = as.numeric(factor(col_types))
row_types = rownames(X)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(X,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

# Eric's method 
phi <- 0.5; k <- 5
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col

## Initialize path parameters and structures
nGamma <- 7
gammaSeq <- 10**seq(0,1,length.out=nGamma)

## Generate solution path
sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq,fraction=0.01)

## Plot validation error
verr <- sol$validation_error
#plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'eric method without smooth function')

heatmap(M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'eric method with smooth function')


# ADMM
## same parameter (nu1,nu2) with Eric
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 1

res_admm1 = bi_ADMM_order(X, nu1, nu2, gamma_1, gamma_2, kk=5, 
                          phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(X, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',gammaSeq[ix]))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('ADMM smooth','nu = ',gammaSeq[ix]))

## try another parameter
nu1 = nu2 = 1
gamma_1 = gamma_2 = 1

res_admm1 = bi_ADMM_order(X, nu1, nu2, gamma_1, gamma_2, 
                          kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(X, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',nu1))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('ADMM smooth','nu = ',nu1))

# AMA
nu1 = nu2 = 0.1
res_ama = bi_ama_order(X, nu1, nu2, gamma_1, gamma_2, 
                       kk = 5, niter = 1000, tol = 1e-5, output = 1)
MM2 = cluster_assign(X, result = res_ama, method = 'AMA')

heatmap(res_ama$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('AMA','nu = ',nu1))

heatmap(MM2$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('AMA smooth','nu = ',nu1))

# change parameter
nu1 = nu2 = 1
res_ama = bi_ama_order(X, nu1, nu2, gamma_1, gamma_2, 
                       kk = 5, niter = 1000, tol = 1e-5, output = 1)
MM2 = cluster_assign(X, result = res_ama, method = 'AMA')

heatmap(res_ama$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('AMA','nu = ',nu1))

heatmap(MM2$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('AMA smooth','nu = ',nu1))
        

                
                