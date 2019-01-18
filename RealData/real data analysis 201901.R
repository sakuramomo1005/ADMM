### real data analysis
# 2019/01/10

### dimension of the data 34* 60
### column clusters: 30 controls and 30 PATs
### row clusters:
    # three clusters: 11 group1, 11 group2, 12 group3
    # or two clusters: 17 group1, 17 group2

# load data
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/Generate_microbiome_compositional_data')

pac = c("survival","plyr","ggplot2","reshape2","phyloseq",'dirmult','gtools')
.tmp = lapply(pac, library, character.only = T)

load(file = 'mice.data.20181031.RData') 

#nsample = 499;ntaxa = 37
map = sample_data(phy.genus)
otu.tab = otu_table(phy.genus)
tax.tab = tax_table(phy.genus) #the lineage of each taxon

#stimate parameters of the dirichlet multinomial distribution
ind = map$Week==3 &  map$Treatment =='Control' & map$Sex == 'F'
ind2 = map$Week==3 & map$Treatment =='PAT' & map$Sex == 'F'

phy.sub = subset_samples(phy.genus,ind)
phy.sub = prune_taxa(rowSums(otu_table(phy.sub))>1,phy.sub)
phy.sub2 = subset_samples(phy.genus,ind2)
phy.sub2 = prune_taxa(rowSums(otu_table(phy.sub2))>1,phy.sub2)

dim(otu_table(phy.genus))
dim(otu_table(phy.sub))  # 32 bacteria, 19 person
dim(t(otu_table(phy.sub)))
dim(otu_table(phy.sub2))  # 22 bacteria, 18 person
dim(t(otu_table(phy.sub2)))

fit = dirmult(t(otu_table(phy.sub)),epsilon=10^(-3),trace=FALSE) 
fit2 = dirmult(t(otu_table(phy.sub2)),epsilon=10^(-3),trace=FALSE) 

# generate simulation count table
num.samples= 50; mu = 10000; 
#num.samples: sample size; mu: total read count for each sample
sim.otu.tab <- simPop(J=30, n=mu, pi=fit$pi,
                      theta=fit$theta)$data  
sim.otu.tab2 <- simPop(J=30, n=mu, pi=fit2$pi,
                       theta=fit2$theta)$data  
dim(sim.otu.tab) 
dim(sim.otu.tab2)  # 22 bacteria

colnames(sim.otu.tab) =  rownames(otu_table(phy.sub))
colnames(sim.otu.tab2) = rownames(otu_table(phy.sub2))

## bind the two group data
microbiome = smartbind(sim.otu.tab,sim.otu.tab2)

## how many columns that only contain 0
sum(apply(microbiome,2,function(x) sum(is.na(x)))>0)
dim(microbiome)

## deal with missing, if one cell is missing, then it is assigned as 0 
microbiome2 = microbiome
for(i in 1:dim(microbiome)[1]){
  for(j in 1:dim(microbiome)[2]){
    if(is.na(microbiome[i,j])==1){
      microbiome2[i,j] = 0
    }
  }
}
sum(apply(microbiome2,2,function(x) sum(is.na(x)))>0)
dim(microbiome2)

# Results of the simulated microbiome data with three clusters
## group the microbiome, 3 groups
sample_variable = sample(1:dim(microbiome2)[2])
sample1 = sample_variable[1:round(dim(microbiome2)[2]/3)]
sample2 = sample_variable[(1+round(dim(microbiome2)[2]/3)):(
  2*round(dim(microbiome2)[2]/3))]

# Scenario 1: times 10 and times 20
microbiome3 = microbiome2
microbiome3[(1+round(dim(microbiome2)[1]/2)):dim(microbiome2)[1],sample1] = 
  microbiome3[(1+round(dim(microbiome2)[1]/2)):dim(microbiome2)[1],sample1] * 100
microbiome3[1:round(dim(microbiome2)[1]/2),sample2] = 
  microbiome3[1:round(dim(microbiome2)[1]/2),sample2] * 20
microbiome3 = t(microbiome3)

# add column group names and row group names
colnames(microbiome3) =c(rep('control',round(dim(microbiome2)[1]/2)),rep('pta',round(dim(microbiome2)[1]/2)))
micronames = rep('group3',dim(microbiome2)[2])
micronames[sample1] = 'group1'
micronames[sample2] = 'group2'
rownames(microbiome3) = micronames

## make columns values sum to 1
microbiome4 = microbiome3
for(i in 1:dim(microbiome3)[2]){
  microbiome4[,i] =  microbiome3[,i]/apply(microbiome3,2,sum)[i]
}

apply(microbiome4,2,sum)
dim(microbiome4)

# heat map for the origin data
col_types = colnames(microbiome4)
col_ty = as.numeric(factor(col_types))
row_types = rownames(microbiome4)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

# Method:
# Eric's 
X = microbiome4
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
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',gammaSeq[ix]))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('ADMM smooth','nu = ',gammaSeq[ix]))

# ADMM change parameter
nu1 = nu2 = 2.1
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',nu1))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM smooth','nu = ',nu1))



# scenario 2: times 100 and times 200

microbiome3 = microbiome2
microbiome3[(1+round(dim(microbiome2)[1]/2)):dim(microbiome2)[1],sample1] = 
  microbiome3[(1+round(dim(microbiome2)[1]/2)):dim(microbiome2)[1],sample1] * 100
microbiome3[1:round(dim(microbiome2)[1]/2),sample2] = 
  microbiome3[1:round(dim(microbiome2)[1]/2),sample2] * 200
microbiome3 = t(microbiome3)

# add column group names and row group names
colnames(microbiome3) =c(rep('control',round(dim(microbiome2)[1]/2)),rep('pta',round(dim(microbiome2)[1]/2)))
micronames = rep('group3',dim(microbiome2)[2])
micronames[sample1] = 'group1'
micronames[sample2] = 'group2'
rownames(microbiome3) = micronames

### make columns values sum to 1
microbiome4 = microbiome3
for(i in 1:dim(microbiome3)[2]){
  microbiome4[,i] =  microbiome3[,i]/apply(microbiome3,2,sum)[i]
}

apply(microbiome4,2,sum)
dim(microbiome4)

# heat map for the origin data
col_types = colnames(microbiome4)
col_ty = as.numeric(factor(col_types))
row_types = rownames(microbiome4)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

# Method:
# Eric's 
X = microbiome4
phi <- 0.5; k <- 5
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col

#### Initialize path parameters and structures
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
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',gammaSeq[ix]))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('ADMM smooth','nu = ',gammaSeq[ix]))

# ADMM change parameter
nu1 = nu2 = 2.1
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',nu1))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM smooth','nu = ',nu1))


# scenario: times 0.2 and times 2


microbiome3 = microbiome2
microbiome3[(1+round(dim(microbiome2)[1]/2)):dim(microbiome2)[1],sample1] = 
  microbiome3[(1+round(dim(microbiome2)[1]/2)):dim(microbiome2)[1],sample1] * 0.2
microbiome3[1:round(dim(microbiome2)[1]/2),sample2] = 
  microbiome3[1:round(dim(microbiome2)[1]/2),sample2] * 2
microbiome3 = t(microbiome3)

# add column group names and row group names
colnames(microbiome3) =c(rep('control',round(dim(microbiome2)[1]/2)),rep('pta',round(dim(microbiome2)[1]/2)))
micronames = rep('group3',dim(microbiome2)[2])
micronames[sample1] = 'group1'
micronames[sample2] = 'group2'
rownames(microbiome3) = micronames

### make columns values sum to 1
microbiome4 = microbiome3
for(i in 1:dim(microbiome3)[2]){
  microbiome4[,i] =  microbiome3[,i]/apply(microbiome3,2,sum)[i]
}

apply(microbiome4,2,sum)
dim(microbiome4)

# heat map for the origin data
col_types = colnames(microbiome4)
col_ty = as.numeric(factor(col_types))
row_types = rownames(microbiome4)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

# Method:
# Eric's 
X = microbiome4
phi <- 0.5; k <- 5
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col

#### Initialize path parameters and structures
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
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',gammaSeq[ix]))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main =paste('ADMM smooth','nu = ',gammaSeq[ix]))

# ADMM change parameter
nu1 = nu2 = 2.1
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 1000,tol = 0.001,output = 1)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM','nu = ',nu1))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = paste('ADMM smooth','nu = ',nu1))


#######  Results of the simulated microbiome data with two clusters

# Scenario 1:

```{R include = FALSE}

## deal with missing
microbiome3 = microbiome2

microbiome3[1:30,1:17] = 2 * microbiome3[1:30,1:17]
microbiome3 = t(microbiome3)

# add column group names and row group names
colnames(microbiome3) = c(rep('control',30),rep('pta',30))
rownames(microbiome3) = c(rep('g1',17),rep('g2',17))

### make columns values sum to 1
microbiome4 = microbiome3
for(i in 1:dim(microbiome3)[2]){
  microbiome4[,i] =  microbiome3[,i]/apply(microbiome3,2,sum)[i]
}

apply(microbiome4,2,sum)
dim(microbiome4)

# heat map for the origin data
col_types = colnames(microbiome4)
col_ty = as.numeric(factor(col_types))
row_types = rownames(microbiome4)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

### Zero percentage:
sum(abs(microbiome4)<0.001)/(34*60)

## Heatmap of raw data
heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

# Eric's method
X = microbiome4
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
###plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method without smooth')

heatmap(M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method with smooth')
# ADMM

## ADMM with the same parameter (nu1, nu2) with Eric's method 
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 2000,tol = 0.001,output = 0)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM','nu = ',round(gammaSeq[ix],2)))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main =paste('ADMM smooth','nu = ',round(gammaSeq[ix],2)))

## ADMM change parameter
nu1 = nu2 = 2.1
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 2000,tol = 0.001,output = 0)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

## try another parameter
heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM','nu = ',nu1))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM smooth','nu = ',nu1))

# Scenario 2:
microbiome3 = microbiome2

microbiome3[1:30,1:17] = 20 * microbiome3[1:30,1:17]
microbiome3 = t(microbiome3)

# add column group names and row group names
colnames(microbiome3) = c(rep('control',30),rep('pta',30))
rownames(microbiome3) = c(rep('g1',17),rep('g2',17))

### make columns values sum to 1
microbiome4 = microbiome3
for(i in 1:dim(microbiome3)[2]){
microbiome4[,i] =  microbiome3[,i]/apply(microbiome3,2,sum)[i]
}

apply(microbiome4,2,sum)
dim(microbiome4)

# heat map for the origin data
col_types = colnames(microbiome4)
col_ty = as.numeric(factor(col_types))
row_types = rownames(microbiome4)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

### Zero percentage:
sum(abs(microbiome4)<0.001)/(34*60)

## Heatmap of raw data
heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

# Eric's method
X = microbiome4
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
###plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method without smooth')

heatmap(M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method with smooth')

# ADMM
## ADMM with the same parameter (nu1, nu2) with Eric's method 
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 2000,tol = 0.001,output = 0)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
       ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
       main = paste('ADMM','nu = ',round(gammaSeq[ix],2)))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main =paste('ADMM smooth','nu = ',round(gammaSeq[ix],2)))

# ADMM change parameter
nu1 = nu2 = 2.1
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 2000,tol = 0.001,output = 0)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM','nu = ',nu1))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM smooth','nu = ',nu1))

# Scenario 3:
microbiome3 = microbiome2

microbiome3[1:30,1:17] = 200 * microbiome3[1:30,1:17]
microbiome3 = t(microbiome3)

# add column group names and row group names
colnames(microbiome3) = c(rep('control',30),rep('pta',30))
rownames(microbiome3) = c(rep('g1',17),rep('g2',17))

### make columns values sum to 1
microbiome4 = microbiome3
for(i in 1:dim(microbiome3)[2]){
microbiome4[,i] =  microbiome3[,i]/apply(microbiome3,2,sum)[i]
}


# heat map for the origin data
col_types = colnames(microbiome4)
col_ty = as.numeric(factor(col_types))
row_types = rownames(microbiome4)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

### Zero percentage:
sum(abs(microbiome4)<0.001)/(34*60)

## Heatmap of raw data
heatmap(microbiome4,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Raw microbiome data')

# Eric's method
X = microbiome4
phi <- 0.5; k <- 5
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col

#### Initialize path parameters and structures
nGamma <- 7
gammaSeq <- 10**seq(0,1,length.out=nGamma)

## Generate solution path
sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq,fraction=0.01)

## Plot validation error
verr <- sol$validation_error
###plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method without smooth')

heatmap(M,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Eric method with smooth')

# ADMM
## ADMM with the same parameter (nu1, nu2) with Eric's method 
nu1 = nu2 = gammaSeq[ix]
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 2000,tol = 0.001,output = 0)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM','nu = ',round(gammaSeq[ix],2)))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main =paste('ADMM smooth','nu = ',round(gammaSeq[ix],2)))

# ADMM change parameter
nu1 = nu2 = 2.1
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(microbiome4, nu1, nu2, gamma_1, gamma_2, kk=5, phi=0.5,niter = 2000,tol = 0.001,output = 0)
MM = cluster_assign(microbiome4, result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM','nu = ',nu1))

heatmap(MM$M,col=hmcols,labRow=NA,labCol=NA,
ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
main = paste('ADMM smooth','nu = ',nu1))
