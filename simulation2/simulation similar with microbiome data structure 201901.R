### 2019/1/13
### simulation 

# The data generation in this simulation follows
# the same structure with microbiome reall data

# load functions 
library(cvxbiclustr)
library(cvxclustr)

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 1115 codes')
source('admm.R')
source('ama.R')

# generate the data

## two treatment group, same dimension with the real microbiome data
con = matrix(rnorm(30*34,-4,2),34,30)
trl = matrix(rnorm(30*34,4,2),34,30)
dat = cbind(con,trl)

## generate three row groups
sam = sample(1:34)
group1 = sam[1:11]
group2 = sam[12:22]
dat[group1,1:30] = dat[group1,1:30] * 8
dat[group2,31:60] = dat[group2,31:60] * 4

## add row group names and column group names
colnames(dat) = c(rep('control',30),rep('treatment',30))
rown = rep('group3',34)
rown[group1] = 'group1'
rown[group2] = 'group2'
rownames(dat) = rown

## make it looks more like count data
dat = round(dat)
dat = ifelse(dat<0,0,dat)

## each column in dat sums up to 1
for(i in 1:dim(dat)[2]){
  dat[,i] =  dat[,i]/apply(dat,2,sum)[i]
}

length(apply(dat,2,sum))
dim(dat)

## draw the heatmap for the raw data
col_types = colnames(dat)
col_ty = as.numeric(factor(col_types))
row_types = rownames(dat)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(dat,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw simulation data')

# Eric method
X = dat
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
plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method without smooth function')

heatmap(M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method with smooth function')

# ADMM

## same parameter with eric's 
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(dat, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',gammaSeq[ix]))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('ADMM smooth','nu = ',gammaSeq[ix]))
