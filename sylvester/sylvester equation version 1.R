### sylvester equation

### version 1

sylvester = function(A, B, C, eps = 0.0001){
  
  library(MASS)
  library(Matrix)
  
  A1 = Schur(A)
  Q1 = A1$Q; R1 = A1$T
  
  A2 = Schur(B)
  Q2 = A2$Q; R2 = A2$T 
  C = t(Q1) %*% C %*% Q2
  
  Rsq = R1 * R1
  I = diag(dim(A)[1])
  
  k = 1
  n = dim(R2)[1]
  
  while(k < n + 1){
    if(k < n){
      if(abs(R2[k+1, k]) < eps){
        left = R1 + R2[k,k] * I
        right = C[,k]
        temp = matrix(0, dim(X)[1],1)
        if(k == 1){
          temp = temp
        }else{
          for(i in 1:(k-1)){
            temp = temp + X[,i] * R2[i,k]
          }
        }
        temp = matrix(temp, dim(C)[1],1)
        X[,k] = ginv(left) %*% (right - temp)
        mytry = myTryCatch(solve(left))
        if(is.null(mytry$error) == 0){er = c(er,tt)}
        k = k+1
      }else{
        r11 = R2[k,k]; r12 = R2[k, k+1]; r21 = R2[k+1, k]; r22 = R2[k+1, k+1]
        temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        if(k == 1){
          temp2 = temp2
          temp3 = temp3
        }else{
          for(i in 1:(k-1)){
            temp2 = temp2 + X[,i] * R2[i,k]
            temp3 = temp3 + X[,i] * R2[i,k+1]}
        }
        temp2 = matrix(temp2, dim(X)[1],1)
        temp3 = matrix(temp3, dim(X)[1],1)
        b1 = C[,k] - temp2 
        b2 = C[,k+1] - temp3
        b1_prime = R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime = R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime = matrix(0, dim(X)[1],2)
        b_prime[,1] = b1_prime; b_prime[,2] = b2_prime
        X[,k:(k+1)] = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                             (r11*r22 - r12*r21) * I) %*% b_prime 
        k = k+2
      }
    }else{
      if(abs(R2[1, k]) > eps){
        left = R1 + R2[k,k] * I
        right = C[,k]
        temp = matrix(0, dim(X)[1],1)
        if(k == 1){
          temp = temp
        }else{
          for(i in 1:(k-1)){
            temp = temp + X[,i] * R2[i,k]
          }
        }
        temp = matrix(temp, dim(C)[1],1)
        X[,k] = ginv(left) %*% (right - temp)
        k = k+1
      }else{
        R22 = R2
        R22 = cbind(R2, rep(0,dim(R2)[1]))
        R22 = rbind(R22,rep(0,dim(R2)[1]+1))
        r11 = R22[k,k]; r12 = R22[k, k+1]; r21 = R22[k+1, k]; r22 = R22[k+1, k+1]
        temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        for(i in 1:(k-1)){
          temp2 = temp2 + X[,i] * R22[i,k]
          temp3 = temp3 + X[,i] * R22[i,k+1]}
        temp2 = matrix(temp2, dim(X)[1],1)
        temp3 = matrix(temp3, dim(X)[1],1)
        b1 = C[,k] - temp2 
        b2 = - temp3
        b1_prime = R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime = R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime = matrix(0, dim(X)[1],2)
        b_prime[,1] = b1_prime; b_prime[,2] = b2_prime
        GOD = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                     (r11*r22 - r12*r21) * I) %*% b_prime 
        X[,k] = GOD[,1]
        k = k+2
      }
    }
  }
  return(Q1 %*% X %*% t(Q2))
}

## Test:

A = matrix(rnorm(3600,0,6),60,60)
B = matrix(rnorm(360*360,0,6),360,360)
X = matrix(rnorm(360*60,0,7),60,360)

C = A%*% X + X%*% B
XX = X

begin = Sys.time()
res = sylvester(A,B,C)
print(sum(res - XX))
end = Sys.time()
end - begin

sums= c()
eps = 0.0001
for(tt in 1:20){
  B = matrix(rnorm(360*360,0,6),360,360)
  A = matrix(rnorm(3600,0,6),60,60)
  X = matrix(rnorm(360*60,0,7),60,360)
  C = A%*% X + X%*% B
  XX = X
  res = sylvester(A,B,C)
  sums = c(sums, sum(XX-res) < eps)
}
sums

