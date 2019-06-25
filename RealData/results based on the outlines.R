# 2019-02-10
# real data analysis

set.seed(123)
# 1. load data
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/Generate_microbiome_compositional_data')
pac = c("survival","plyr","ggplot2","reshape2","phyloseq",'dirmult','gtools')
.tmp = lapply(pac, library, character.only = T)
load(file = 'mice.data.20181031.RData') 

# nsample = 499; 499 mice
# ntaxa = 37; bacterial taxonomy
map = sample_data(phy.genus)
rawdata = otu_table(phy.genus)

# covariates of interest:
# 1. Treatment: STAT; Control; PAT
# 2. Sex: F;M
# 3. Week: 3;6;10;13
# 4. Diabetes: 1; 0; NA

# check the number of sub-datasets 
ind1 = map$Week==3 &  map$Treatment =='Control'
ind2 = map$Week==6 &  map$Treatment =='Control'
ind3 = map$Week==10 &  map$Treatment =='Control'
ind4 = map$Week==13 &  map$Treatment =='Control'

ind5 = map$Week==3 &  map$Treatment =='PAT'
ind6 = map$Week==6 &  map$Treatment =='PAT'
ind7 = map$Week==10 &  map$Treatment =='PAT'
ind8 = map$Week==13 &  map$Treatment =='PAT'

phy.sub1 = subset_samples(phy.genus,ind1)
phy.sub2 = subset_samples(phy.genus,ind2)
phy.sub3 = subset_samples(phy.genus,ind3)
phy.sub4 = subset_samples(phy.genus,ind4)
phy.sub5 = subset_samples(phy.genus,ind5)
phy.sub6 = subset_samples(phy.genus,ind6)
phy.sub7 = subset_samples(phy.genus,ind7)
phy.sub8 = subset_samples(phy.genus,ind8)

dat1 = otu_table(phy.sub1)
dat2 = otu_table(phy.sub2)
dat3 = otu_table(phy.sub3)
dat4 = otu_table(phy.sub4)
dat5 = otu_table(phy.sub5)
dat6 = otu_table(phy.sub6)
dat7 = otu_table(phy.sub7)
dat8 = otu_table(phy.sub8)

dim(dat1) + # 47
  dim(dat2) + # 45
  dim(dat3) + # 37
  dim(dat4) + # 36
  dim(dat5) + # 42
  dim(dat6) +# 42
  dim(dat7) +# 31
  dim(dat8) # 32

sum(dat1 == 0) # 694
sum(dat2 == 0) # 641
sum(dat3 == 0) # 444
sum(dat4 == 0) # 411
sum(dat5 == 0) # 1128
sum(dat6 == 0) # 1081
sum(dat7 == 0) # 650
sum(dat8 == 0) # 648

# Here is the number of subjects in different groups at different time.
a = c(dim(dat1)[2],dim(dat2)[2],dim(dat3)[2],dim(dat4)[2])
b = c(dim(dat5)[2],dim(dat6)[2],dim(dat7)[2],dim(dat8)[2])
dms = data.frame(rbind(a,b))
rownames(dms) = c('Control','PAT')
colnames(dms) = c('Week 3','Week 6','Week 10','Week 13')
kable(dms)

# Here is the number of 0s in the counts for different groups at different time.
a = c(sum(dat1 == 0), sum(dat2 == 0),sum(dat3 == 0),sum(dat4 == 0))
b = c(sum(dat5 == 0), sum(dat6 == 0),sum(dat7 == 0),sum(dat8 == 0))
mis = data.frame(rbind(a,b))
rownames(mis) = c('Control','PAT')
colnames(mis) = c('Week 3','Week 6','Week 10','Week 13')
kable(mis)

# i would like to use the week 12 since there is less 0s in dataset 
ind4 = map$Week==13 &  map$Treatment =='Control'
ind8 = map$Week==13 &  map$Treatment =='PAT'
phy.sub4 = subset_samples(phy.genus,ind4)
phy.sub8 = subset_samples(phy.genus,ind8)
# control group
con = otu_table(phy.sub4)
con = data.frame(con)
# pat group
pat = otu_table(phy.sub8)
pat = data.frame(pat)
dim(pat);dim(con)

dat = cbind(pat, con)
# change the group names
colnames(dat) = rep(c('pat','con'),times = c(dim(pat)[2],dim(con)[2]))
# remove the rows and columns with all 0s.
dat = dat[-which(apply(dat,1,sum) == 0),]
dim(dat) # 36 68


## 2. Pure real data (without any change on the real data)
#### 2.1 Whether the control group and the PAT group are different?

micro = rep(rownames(dat),dim(dat)[2]) 
value = unlist(dat)
group = rep(c('pat','con'),
            times = c(sum(colnames(dat) =='pat') * dim(dat)[1],
            sum(colnames(dat) =='con') * dim(dat)[1]))
            
diff1 = data.frame(group = group, micro = micro, value = value)
summary(aov(value ~ micro + group, data = diff1))
summary(aov(value ~ group, data = diff1))

# The counts value in the two groups do not have significant differences

## 2.1.1 Draw heatmap
col_types = colnames(dat)
col_ty = as.numeric(factor(col_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)
heatmap(as.matrix(dat),col=hmcols,
        ColSideCol=cols[col_ty],
        main = 'Fig 1: Raw data (control & PAT)',scale ='none')

dat = dat/ matrix(apply(dat,2,sum),dim(dat)[1],dim(dat)[2],byrow =TRUE)
apply(dat,2,sum)

heatmap(as.matrix(dat),col=hmcols,
        ColSideCol=cols[col_ty],
        main = 'Fig 1: Raw data (control & PAT)')

# 2.2 Let's try the methods:
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 1115 codes')
source('admm.R')
source('ama.R')

library(fossil)
library(cvxbiclustr)
library(cvxclustr)

## 2.2.1. Eric's method
X = as.matrix(dat)
phi = 0.5; k = 5
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
#plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,labCol=NA,
        ColSideCol=cols[col_ty],
        main = 'Fig2: Eric no smooth for raw data')

heatmap(M,col=hmcols,labCol=NA,
        ColSideCol=cols[col_ty],
        main = 'Fig3: Eric with smooth for raw data')

## 2.2.2 ADMM
nu1 = nu2 = 2
gamma_1 = gamma_2 = 5;k_row = k_col = 5

res_admm1 = bi_ADMM_order(as.matrix(dat), nu1, nu2, gamma_1, gamma_2,
                          kk=5, phi=0.5,niter = 1000,
                          tol = 0.00001,output = 1)
MM = cluster_assign(as.matrix(dat), result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,labCol=NA,
        ColSideCol=cols[col_ty],
        main = 'Fig4: ADMM no smooth for raw data')

heatmap(MM$M,col=hmcols,labCol=NA,
        ColSideCol=cols[col_ty],
        main = 'Fig5: ADMM with smooth for raw data')


# B. Semi real data on the subjects.
# Highlight: Select one group, e.g. the control group. First, test whether the distributions of microbiomes are different within the group. If not, select half of the subjects as group 1 and the rest as group 2. Increase or decrease some microbiome values in group 1. Then apply our convex bi-clustering method to test wether we can separate group 1 and group 2.
# 
# More details:
#   1. Within one intervention group, say control group, we randomly selected two subjects to compare their microbiome counts through ANOVA to see whether there are differences within the group.
# 2. Sort the columns of the dataset. Make the values in each column have the similar order.
# 3. Divide the subjects into two groups: group 1, group 2.
# 4. Multiple the part of the microbiomes in group 1 with an arbitrary value larger than 1.
# 5. Multiple the part of the microbiomes in group 2 with an arbitrary value less than 1.
# 6. Apply ANOVA again to check whether those two groups have differences and check whether the microbiomes within one group have differences.
# 7. Apply the convex bi-clustering methods to check whether they can identify different groups of subjects and different groups of microbiomes.

# B.1

# comparision of the microbiome counts for subjects within the control group
which(apply(con,1,sum) == 0)
dim(con)
con[-which(apply(con,1,sum) == 0),]
dim(con[-which(apply(con,1,sum) == 0),])
con = con[-which(apply(con,1,sum) == 0),]
g1 = 1:round(dim(con)[2]/2)
g2 = (round(dim(con)[2]/2)+1): dim(con)[2]
colnames(con) = rep(c('g1','g2'),times = c(length(g1),length(g2)))
dat = con

# pairwise comparison 
paov = c(); pt = c()
n1 = c(); n2 = c()
for(i in 1:(dim(dat)[2]-1)){
  for(j in (i+1):dim(dat)[2]){
    value = c(dat[,i],dat[,j])
    group = c(rep(1,length(dat[,i])),rep(2,length(dat[,j])))
    test = aov(value~group)
    paov = c(paov, summary(test)[[1]]$`Pr(>F)`[1])
    pt = c(pt, t.test(dat[,i],dat[,j])$p.value)
    n1 = c(n1, i); n2 = c(n2, j)
  }
}
pair_compare = data.frame(subj1 = n1, subj2 = n2, pvalue = p)
head(pair_compare)
min(pair_compare$pvalue)
min(p)
 
# group comparsion
micro = rep(rownames(dat),dim(dat)[2]) 
value = unlist(dat)
group = rep(c('g1','g2'),
            times = c(sum(colnames(dat) =='g1') * dim(dat)[1],
                      sum(colnames(dat) =='g2') * dim(dat)[1]))

diff2 = data.frame(group = group, micro = micro, value = value)
summary(aov(value ~ micro + group, data = diff2))
summary(aov(value ~ group, data = diff2))

# Sort the columns of the dataset.
which(apply(con,2,sum) == max(apply(con,2,sum)))
dat2 = con[order(con[,1],con[,2],con[,3],con[,4],con[,5],con[,6],
                 con[,7],con[,8],con[,9],con[,10],con[,11],con[,12],
                 con[,13],con[,14],con[,15],con[,16],con[,17],con[,18]),]

col_types = colnames(dat2)
col_ty = as.numeric(factor(col_types))
heatmap(as.matrix(dat2),col=hmcols,labCol = NA,
        ColSideCol=cols[col_ty],
        main = 'Fig6: Raw data (control)')

order(apply(dat2,1,sum))

dim(dat2)
dat2_temp = dat2
dat2_temp[18:34,1:18] = dat2_temp[18:34,1:18] * 50
#dat2_temp[18:34,19:36] = dat2_temp[1:17,19:36] * 10000
dat2_temp = as.matrix(dat2_temp)
rownames(dat2_temp)[1:18] = c('1')
rownames(dat2_temp)[19:34] = c('2')
row_types = rownames(dat2_temp)
row_ty = as.numeric(factor(row_types))

# check
order(apply(dat2_temp,1,sum))
heatmap(as.matrix(dat2_temp),col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig7: Semi-raw data1 (control)')

# micro = rep(rownames(dat2_temp),dim(dat2_temp)[2])
# value = unlist(dat2_temp)
# group = rep(c('g1','g2'),
#             times = c(sum(colnames(dat2_temp) =='g1') * dim(dat2_temp)[1],
#                       sum(colnames(dat2_temp) =='g2') * dim(dat2_temp)[1]))
# 
# diff3 = data.frame(group = group, micro = micro, value = value)
# summary(aov(value ~ micro + group, data = diff2))
# summary(aov(value ~ group, data = diff2))


## B. 2.2.1. Eric's method
dat2_temp = dat2_temp/ matrix(apply(dat2_temp,2,sum),dim(dat2_temp)[1],
                              dim(dat2_temp)[2],byrow =TRUE)

heatmap(as.matrix(dat2_temp),col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig7: Semi-raw data1 (control)')


apply(dat2_temp,2,sum)
X = as.matrix(dat2_temp)
phi = 0.5; k = 5
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
#plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig8: Eric no smooth for data2')

heatmap(M,col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig9: Eric with smooth for data2')

## B.2.2.2 ADMM

nu1 = nu2 =2
gamma_1 = gamma_2 = 0.5;k_row = k_col = 5
res_admm1 = bi_ADMM_order(as.matrix(dat2_temp), nu1, nu2, gamma_1, gamma_2,
                          kk=5, phi=0.5,niter = 2000,
                          tol = 0.00001,output = 1)
MM = cluster_assign(as.matrix(dat2_temp), result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig10: ADMM no smooth for data 2')

heatmap(MM$M,col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig11: ADMM with smooth for data2')

# C. Semi real data on both subjects and microbiomes.
# We would like to see whether our method can this help us separate microbiome groups when we have already known the treatment groups.
# 
# Steps:
#   1. Combine subjects and microbiomes from PAT group and control group.
# 2. Fill NAs with 0s.
# 3. Sort the columns of the dataset. Make the values in each column have the similar order.
# 4. Divide the subjects into two groups: group 1, group 2. 
# 5. Multiple the part of the microbiomes in group 1 with an arbitrary value larger than 1.
# 6. Multiple the part of the microbiomes in group 2 with an arbitrary value less than 1.
# 7. Apply the convex bi-clustering methods to check whether they can identify different groups of subjects and different groups of microbiomes.

con = otu_table(phy.sub4)
con = data.frame(con)
# pat group
pat = otu_table(phy.sub8)
pat = data.frame(pat)

dat = cbind(pat, con)
# change the group names
colnames(dat) = rep(c('pat','con'),times = c(dim(pat)[2],dim(con)[2]))
# remove the rows and columns with all 0s.
dat = dat[-which(apply(dat,1,sum) == 0),]
dim(dat) # 36 68

dat3 = dat[order(dat[,1],dat[,2],dat[,3],dat[,4],dat[,5],dat[,6],
                 dat[,7],dat[,8],dat[,9],dat[,10],dat[,11],dat[,12],
                 dat[,13],dat[,14],dat[,15],dat[,16],dat[,17],dat[,18]),]

dim(dat3)
dat2_temp = dat3
dat2_temp[18:34,1:34] = dat2_temp[18:34,1:34] * 50
#dat2_temp[18:34,19:36] = dat2_temp[1:17,19:36] * 10000
dat2_temp = as.matrix(dat2_temp)
rownames(dat2_temp)[1:17] = c('1')
rownames(dat2_temp)[18:34] = c('2')
row_types = rownames(dat2_temp)
row_ty = as.numeric(factor(row_types))
col_types = colnames(dat2_temp)
col_ty = as.numeric(factor(col_types))

# check
order(apply(dat2_temp,1,sum))
heatmap(as.matrix(dat2_temp),col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig12: Semi-raw data (control and PAT)')


## B. 2.2.1. Eric's method
dat2_temp = dat2_temp/ matrix(apply(dat2_temp,2,sum),dim(dat2_temp)[1],
                              dim(dat2_temp)[2],byrow =TRUE)
apply(dat2_temp,2,sum)
X = as.matrix(dat2_temp)
phi = 0.5; k = 5
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
#plot(verr)

## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- sol$groups_row[[ix]]
groups_col <- sol$groups_col[[ix]]
M <- biclust_smooth(X,groups_row,groups_col)

heatmap(sol$U[[ix]],col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig13: Eric no smooth for data2')

heatmap(M,col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig14: Eric with smooth for data2')

## B.2.2.2 ADMM
set.seed(1234)
nu1 = nu2 =2
gamma_1 = gamma_2 = 5;k_row = k_col = 5
res_admm1 = bi_ADMM_order(as.matrix(dat2_temp), nu1, nu2, gamma_1, gamma_2,
                          kk=5, phi=0.5,niter = 2000,
                          tol = 0.00001,output = 1)
MM = cluster_assign(as.matrix(dat2_temp), result = res_admm1, method = 'ADMM')

heatmap(res_admm1$A,col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig15: ADMM no smooth for data 2')

heatmap(MM$M,col=hmcols,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Fig16: ADMM with smooth for data2')