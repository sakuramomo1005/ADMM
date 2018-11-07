# update v z 










X = matrix(1:30,6,5)
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])

A =  matrix(1:30,6,5)
# update lambda 1
for(i in 1:dim(eplison_n)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  a_l1 = matrix(A[l1,],p,1)
  a_l2 = matrix(A[l2,],p,1)
  lambda_1[,i] = lambda_1[,i] + nu1 * (v[,i] - a_l1 + a_l2)
} 
# update lambda 2
for(i in 1:dim(eplison_p)[1]){
  l1=eplison_p[i,'l1']
  l2=eplison_p[i,'l2']
  a_k1 = matrix(A[,l1],n,1)
  a_k2 = matrix(A[,l2],n,1)
  lambda_2[,i] = lambda_2[,i] + nu2 * (z[,i] - a_k1 + a_k2)
}







X = matrix(1:30,6,5)
A =  matrix(1:30,6,5)
n = dim(X)[1]; p = dim(X)[2]
eplison_p = L_num(p)
eplison_n = L_num(n)

lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])

nu1 = nu2 = 1

# update lambda 1
for(i in 1:dim(eplison_n)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  a_l1 = matrix(A[l1,],p,1)
  a_l2 = matrix(A[l2,],p,1)
  lambda_1[,i] = lambda_1[,i] + nu1 * (v[,i] - a_l1 + a_l2)
} 
# update lambda 2
for(i in 1:dim(eplison_p)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  a_k1 = matrix(A[,l1],n,1)
  a_k2 = matrix(A[,l2],n,1)
  lambda_2[,i] = lambda_2[,i] + nu2 * (z[,i] - a_k1 + a_k2)
}

















a2 = matrix(0,dim(lambda_1)[1],dim(lambda_1)[2])
for(i in 1:dim(eplison_n)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  
  a_l1 = matrix(A[l1,],p,1)
  a_l2 = matrix(A[l2,],p,1)
  print(a_l1)
  a2[,i] = a_l2
  lambda_1[,i] = lambda_1[,i] + nu1 * (v[,i] - a_l1 + a_l2)
} 

# update lambda 2
a3 = matrix(0,n,dim(eplison_p)[1])
for(i in 1:dim(eplison_p)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  a_k1 = matrix(A[,l1],n,1)
  a_k2 = matrix(A[,l2],n,1)
  a3[,i] = a_k2
  lambda_2[,i] = lambda_2[,i] + nu2 * (z[,i] - a_k1 + a_k2)
}
a2 = matrix(0,dim(lambda_1)[1],dim(lambda_1)[2])
for(i in 1:dim(eplison_n)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  a_l1 = matrix(A[l1,],p,1)
  a_l2 = matrix(A[l2,],p,1)
  print(a_l1)
  a2[,i] = a2[,i] +  a_l2
} 
lambda_1

t(A) %*% cbind(diag(1,4,1), diag(1,4,2), diag(1,4,3))


n=3
al1 = matrix(0,n,n*(n+1)/2)
count = 0
for(i in 1:(n)){
  al1[,(count+1):(count+i)] = diag(1,n,i)
  count = count + i
}

count = 0
al2 = matrix(0,n,n*(n+1)/2)
for(i in 1:(n)){
  print(i)
  temp = matrix(0,n,i)
  temp[i,] = 1
  al2[,(count+1):(count+i)] = temp
  print(count)
  count = count+i
}

lambda_1 = matrix(1,p,dim(eplison_n)[1])
lambda_1 + nu1 * (v - t(A[1:(n-1),]) %*% al1 + t(A[2:n,]) %*% al2)


count = 0
al1 = al2 = matrix(0,n,n*(n+1)/2)
for(i in 1:n){
  al1[,(count+1):(count+i)] = diag(1,n,i)
  print(al1)
  temp = matrix(0,n,i)
  temp[i,] = 1
  al2[,(count+1):(count+i)] = temp
  count = count+i
}


count = 0
al1 = al2 = matrix(0,p,p*(p-1)/2)
for(i in 1:(p-1)){
  al1[,(count+1):(count+i)] = diag(1,p,i)
  temp = matrix(0,p,i)
  temp[i+1,] = 1
  al2[,(count+1):(count+i)] = temp
  count = count+i
}

t(A[,2:p]) %*% al1

lambda_1 + nu1 * (v - t(A[1:(n-1),]) %*% al1 + t(A[2:n,]) %*% al2)
lambda_2 + nu2 * (z - (A) %*% t(al1) + (A) %*% al2)

a2 = matrix(0,dim(lambda_2)[1],dim(lambda_2)[2])
for(i in 1:dim(eplison_p)[1]){
  l1=eplison_n[i,'l1']
  l2=eplison_n[i,'l2']
  a_k1 = matrix(A[,l1],n,1)
  a_k2 = matrix(A[,l2],n,1)
  print(a_k2)
}


lambda_1 = matrix(1,p,dim(eplison_n)[1])
v = matrix(1,p,dim(eplison_n)[1])
lambda_2 = matrix(1,n,dim(eplison_p)[1])
z = matrix(1,n,dim(eplison_p)[1])


# l

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
lambda_1 + nu1 * (v - t(A) %*% al1 + t(A) %*% al2)

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
lambda_2 + nu2 * (z - A %*% al1 + A %*% al2)