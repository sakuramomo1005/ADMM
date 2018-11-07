### 5. update lambda

### 5. update lambda
update_lambda_old=function(X, A, nu1, nu2,v, z){
  n = dim(X)[1]; p = dim(X)[2]
  eplison_p = L_num(p)
  eplison_n = L_num(n)
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
  return(lambda=list(lambda_1,lambda_2))
}

update_lambda=function(X, A, nu1, nu2,v, z){
  n = dim(X)[1]; p = dim(X)[2]
  eplison_p = L_num(p)
  eplison_n = L_num(n)
  # update lambda 1
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
  lambda_1 = lambda_1 + nu1 * (v - t(A) %*% al1 + t(A) %*% al2)
  
  # update lambda 2
  count = 0
  al1 = al2 = matrix(0,p,p*(p-1)/2)
  for(i in 1:(p-1)){
    al1[,(count+1):(count+i)] = diag(1,p,i)
    temp = matrix(0,p,i)
    temp[i+1,] = 1
    al2[,(count+1):(count+i)] = temp
    count = count+i
  }
  lambda_2 = lambda_2 + nu2 * (z - A %*% al1 + A %*% al2)
  return(lambda=list(lambda_1,lambda_2))
}
