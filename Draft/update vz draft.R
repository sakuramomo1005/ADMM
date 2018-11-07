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
al1 = al2 = matrix(0,p,p*(p-1)/2)
for(i in 1:(p-1)){
  al1[,(count+1):(count+i)] = diag(1,p,i)
  temp = matrix(0,p,i)
  temp[i+1,] = 1
  al2[,(count+1):(count+i)] = temp
  count = count+i
}
ak1 = A %*% al1; ak2 = A %*% al2

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







nu1 = 1; nu2 = 1; 
X = matrix(1:30,6,5)
A = X
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])



prox0 = function(v,sigma){
  return(max(1-sigma/sum(v %*% t(v)),0) %*% v)
}

prox = function(v,sigma){
  return(max(1-sigma/sum(t(v) %*% v),0) * v)
}

sigma_1l/sum(t(v_temp) %*% v_temp)

v = 1:5
prox(v, 0.1)
t(v_temp) %*% v_temp

v_temp %*% t(v_temp)

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

update_vz = function(X, A, lambda_1, lambda_2, gamma_1, gamma_2, nu1, nu2){
  
  n=dim(X)[1]; p=dim(X)[2]
  eplison_n=L_num(n)
  eplison_p=L_num(p)
  
  w <- dist_weight(X / sqrt(p),phi = 0.5)
  w_l <- knn_weights(w,k = 1,n)
  
  u = dist_weight(t(X) / sqrt(n),phi = 0.5)
  u_k <- knn_weights(u,k = 1,p)
  
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    a_l1 = A[l1,]; a_l2 = A[l2,]
    v_temp = a_l1 - a_l2 - 1/nu1 * lambda_1[,i]
    sigma_1l = gamma_1 * w_l[i]/nu1
    #print(v_temp)
    
    #(max(1-sigma/
      #    print( sum(t(v_temp) %*% (v_temp)) )
           #,0) %*% v)
        print(max(1-sigma_1l/sum(t(v_temp) %*% v_temp),0))
    
    v[,i] = prox(v_temp,sigma_1l)
    #print(v[,i])
  }
  
  sigma_1l = c(); v_temp = c()
  for(i in 1:dim(eplison_n)[1]){
    l1=eplison_n[i,'l1']
    l2=eplison_n[i,'l2']
    a_l1 = A[l1,]; a_l2 = A[l2,]
    v_temp = cbind(v_temp,a_l1 - a_l2 - 1/nu1 * lambda_1[,i])
    sigma_1l = c(sigma_1l,gamma_1 * w_l[i]/nu1)
    v[,i] = prox(v_temp,sigma_1l)
  }
  
  sigma_1 = gamma_1 * w_l/nu1
  vtemp = al1 - al2 - 1/nu1 * lambda_1
  
  sigma = matrix(sigma_1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE)
  
  temp1 = ifelse((1 - sigma_1/apply(vtemp^2,2,sum)) < 0, 0,1 - sigma_1/apply(vtemp^2,2,sum))
  temp2 = matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp 
  prox
  
  
  
  for(i in 1:dim(eplison_p)[1]){
    l1=eplison_p[i,'l1']
    l2=eplison_p[i,'l2']
    a_l1 = A[,l1]; a_l2 = A[,l2]
    v_temp = a_l1 - a_l2 - 1/nu2 * lambda_2[,i]
   # print(a_l1 - a_l2 - 1/nu2 * lambda_2[,i])
    #u_k = exp(-0.5 * (t(X[,l1] - X[,l2]) %*% (X[,l1] - X[,l2])))
    sigma_2k = gamma_2 * u_k[i]/nu2
    print( sigma_2k)
    z[,i] = prox(v_temp,sigma_2k)
  }
  
  return(list(v = v, z = z))
}



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
