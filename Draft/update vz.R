# functions 
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

prox = function(v,sigma){
  return(max(1-sigma/sum(t(v) %*% v),0) * v)
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

# datas 


nu1 = 1; nu2 = 1; gamma_1 = 0.5; gamma_2 = 0
X = matrix(1:30,6,5)
A = X
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])



count = 0
al1 = al2 = matrix(0,n,n*(n-1)/2)
for(i in 1:(n-1)){
  al1[,(count+1):(count+i)] = diag(1,n,i)
  temp = matrix(0,n,i)
  temp[i,] = 1
  al2[,(count+1):(count+i)] = temp
  count = count+i
}
al2 = rbind(al2[n,],al2[1:(n-1),])
al1 = t(A) %*% al1; al2 = t(A) %*% al2

# k
count = 0
ak1 = ak2 = matrix(0,p,p*(p-1)/2)
for(i in 1:(p-1)){
  ak1[,(count+1):(count+i)] = diag(1,p,i)
  temp = matrix(0,p,i)
  temp[i+1,] = 1
  ak2[,(count+1):(count+i)] = temp
  count = count+i
}
ak1 = A %*% ak1; ak2 = A %*% ak2



w = dist_weight(X / sqrt(p),phi = 0.5)
w_l = knn_weights(w,k = 1,n)
u = dist_weight(t(X) / sqrt(n),phi = 0.5)
u_k = knn_weights(u,k = 1,p)

sigma_1 = gamma_1 * w_l/nu1
vtemp = al1 - al2 - 1/nu1 * lambda_1

temp1 = ifelse((1 - sigma_1/apply(vtemp^2,2,sum)) < 0, 0,1 - sigma_1/apply(vtemp^2,2,sum))
temp2 = matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp 

v = temp2

ztemp = ak1 - ak2 - 1/nu2 * lambda_2
sigma_2 = gamma_2 * u_k/nu2

temp3 = ifelse((1 - sigma_2/apply(ztemp^2,2,sum)) < 0, 0,1 - sigma_2/apply(ztemp^2,2,sum))
temp4 = matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp 

z = temp4


update_vz = function(){
  n=dim(X)[1]; p=dim(X)[2]
  eplison_n=L_num(n)
  eplison_p=L_num(p)
  
  w = dist_weight(X / sqrt(p),phi = 0.5)
  w_l = knn_weights(w,k = 1,n)
  u = dist_weight(t(X) / sqrt(n),phi = 0.5)
  u_k = knn_weights(u,k = 1,p)
  
  sigma_1 = gamma_1 * w_l/nu1
  vtemp = al1 - al2 - 1/nu1 * lambda_1
  
  temp1 = ifelse((1 - sigma_1/apply(vtemp^2,2,sum)) < 0, 0,1 - sigma_1/apply(vtemp^2,2,sum))
  temp2 = matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp 
  
  v = temp2
  
  ztemp = ak1 - ak2 - 1/nu2 * lambda_2
  sigma_2 = gamma_2 * u_k/nu2
  
  temp3 = ifelse((1 - sigma_2/apply(ztemp^2,2,sum)) < 0, 0,1 - sigma_2/apply(ztemp^2,2,sum))
  temp4 = matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp 
  
  z = temp4
}






gamma_1 = 0.5; gamma_2 = 0

nu1 = 1; nu2 = 1; 
X = matrix(1:3000,60,50)
A = X
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])


update_vz = function(X, A, lambda_1, lambda_2, gamma_1, gamma_2, nu1, nu2,n,p,eplison_n,eplison_p){
  w = dist_weight(X / sqrt(p),phi = 0.5)
  w_l = knn_weights(w,k = 1,n)
  u = dist_weight(t(X) / sqrt(n),phi = 0.5)
  u_k = knn_weights(u,k = 1,p)
  
  
  sigma_1 = gamma_1 * w_l/nu1
  vtemp = al1 - al2 - 1/nu1 * lambda_1
  
  temp1 = ifelse((1 - sigma_1/apply(vtemp^2,2,sum)) < 0, 0,1 - sigma_1/apply(vtemp^2,2,sum))
  temp2 = matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp 
  
  v = temp2
  
  ztemp = ak1 - ak2 - 1/nu2 * lambda_2
  sigma_2 = gamma_2 * u_k/nu2
  
  temp3 = ifelse((1 - sigma_2/apply(ztemp^2,2,sum)) < 0, 0,1 - sigma_2/apply(ztemp^2,2,sum))
  temp4 = matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp 
  
  z = temp4
  return(list(v=v,z=z))
}

alk = function(A,n,p){
  count = 0
  al1 = al2 = matrix(0,n,n*(n-1)/2)
  for(i in 1:(n-1)){
    al1[,(count+1):(count+i)] = diag(1,n,i)
    temp = matrix(0,n,i)
    temp[i,] = 1
    al2[,(count+1):(count+i)] = temp
    count = count+i
  }
  al2 = rbind(al2[n,],al2[1:(n-1),])
  al1 = t(A) %*% al1; al2 = t(A) %*% al2
  
  # k
  count = 0
  ak1 = ak2 = matrix(0,p,p*(p-1)/2)
  for(i in 1:(p-1)){
    ak1[,(count+1):(count+i)] = diag(1,p,i)
    temp = matrix(0,p,i)
    temp[i+1,] = 1
    ak2[,(count+1):(count+i)] = temp
    count = count+i
  }
  ak1 = A %*% ak1; ak2 = A %*% ak2
}