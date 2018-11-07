# 1 MATLAB AND R have different result
# 2 when it is ginv, the results not good
# 3 how to determine the k = n? 


A = matrix(c(0.2769  ,  0.0971,
             0.0462   , 0.8235),2,2,byrow = TRUE)
D = matrix(c(0.7547   , 0.6551  ,  0.4984,
             0.2760    ,0.1626   , 0.9597,
             0.6797    ,0.1190   , 0.3404),3,3,byrow = TRUE)
B = t(D)
X = matrix(c(0.5853  ,  0.7513  ,  0.5060,
             0.2238   , 0.2551 ,   0.6991),2,3,byrow=TRUE)

A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T

R2[3,2] = 0.1
B = Q2 %*% R2 %*% t(Q2)

C = A %*% X +X %*% B

XX = X 


A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T

library(Matrix)
library(MASS)

A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T

Q1 %*% R1 %*% t(Q1)






A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T

C = t(Q1) %*% C %*% Q2

Rsq = R1 * R1
I = diag(dim(A)[1])

X = matrix(0, dim(A)[1],dim(B)[1])

k = 1
n = dim(R2)[1]

while(k < n+1){
  k0 = k
  if(k0 == n){
    b = C[,k] - X[,1:(k-1)] %*% matrix(R2[1:(k-1),k],k-1,1)
    left = R1 + R2[k,k] * I
    X[,k] = ginv(left) %*% (b)
    k = k+1
    }
  if( k0 < n){
    if(R2[k+1, k] < eps){
      if(k == 1){
        b = C[,1] - X[,1] * R2[1,1]
      }else{
        b = C[,k] - X[,1:(k-1)] %*% matrix(R2[1:(k-1),k],k-1,1)
      }
      left = R1 + R2[k,k] * I
      X[,k] = ginv(left) %*% (b)
      k = k+1
    }else{
      r11 = R2[k,k]; r12 = R2[k, k+1]; r21 = R2[k+1, k]; r22 = R2[k+1, k+1]
      b = C[,k:(k+1)] - X[,1:(k-1)] %*% R2[1:(k-1),k:(k+1)]
      p1 = R1 %*% b[,1] + r22 * b[,1] - r21 * b[,2]
      p2 = R1 %*% b[,2] + r11 * b[,2] - r12 * b[,1]
      b = cbind(p1,p2)
      X[,k:(k+1)] = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                           (r11*r22 - r12*r21) * I) %*% b 
      k = k+2
  }}
}
Q1 %*% X %*% t(Q2)





B = matrix(rnorm(16,4,6),4,4)
A = matrix(rnorm(9,4,6),3,3)
X = matrix(rnorm(12,6,7),3,4)
C = A%*% X + X%*% B

A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T
  
A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T

R2[4,3] = 0.1 
B = Q2%*%R2%*%t(Q2)  

C = A%*% X + X%*% B
XX = X

C = t(Q1) %*% C %*% Q2
  
  Rsq = R1 * R1
  I = diag(dim(A)[1])
  
  X[,1] = ginv(R1 + R2[1,1] * I) %*% C[,1]
  k = 1
  n = dim(R2)[1]
  
while(k < n + 1){
  if(k < n){
    if(R2[k+1, k] < eps){
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
      X[,k] = ginv(left) %*% (right + temp)
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
  
round(Q1 %*% X %*% t(Q2))

# CORRECT EXAMPLE: 

A_C = A; B_C = B; X_C = XX; Cc = A_C%*% X_C + X_C%*% B_C

A1 = Schur(A_C)
Q1 = A1$Q; R1 = A1$T

A2 = Schur(B_C)
Q2 = A2$Q; R2 = A2$T

XX_C = X_C

C = t(Q1) %*% Cc %*% Q2

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

Q1 %*% X %*% t(Q2)
XX_C




XX = Q1 %*% X %*% t(Q2)

XXX = solve(Q1) %*% XX %*% solve(t(Q2))



R22 = R2
R22 = cbind(R2, rep(0,dim(R2)[1]))
R22 = rbind(R22,rep(0,dim(R2)[1]+1))

r11 = R22[k,k]; r12 = R22[k, k+1]; r21 = R22[k+1, k]; r22 = R22[k+1, k+1]
temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
if(k == 1){
  temp2 = temp2
  temp3 = temp3
}else{
  for(i in 1:(k-1)){
    temp2 = temp2 + X[,i] * R22[i,k]
    temp3 = temp3 + X[,i] * R22[i,k+1]}
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




# 再搞搞还对不对呢
B = matrix(rnorm(16,0,1),8,8)
A = matrix(rnorm(9,0,1),6,6)
X = matrix(rnorm(12,0,1),6,8)
C = A%*% X + X%*% B
XX = X

A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T

A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T

XX= X

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


Q1 %*% X %*% t(Q2)
XX


XXX = solve(Q1) %*% XX %*% solve(t(Q2))




B = matrix(rnorm(25,0,6),5,5)
A = matrix(rnorm(9,0,6),6,6)
X = matrix(rnorm(30,0,7),6,5)
C = A%*% X + X%*% B
XX = X

A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T

A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T 
C = t(Q1) %*% C %*% Q2

Rsq = R1 * R1
I = diag(dim(A)[1])


XXX = solve(Q1) %*% XX %*% solve(t(Q2))

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

XXX
X

sum(XXX -X ) < eps


myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}



er = c()
sums= c()
for(tt in 1:10){
  print(tt)
  B = matrix(rnorm(360*360,0,6),360,360)
  A = matrix(rnorm(3600,0,6),60,60)
  X = matrix(rnorm(360*60,0,7),60,360)
  C = A%*% X + X%*% B
  XX = X
  
  A1 = Schur(A)
  Q1 = A1$Q; R1 = A1$T
  
  A2 = Schur(B)
  Q2 = A2$Q; R2 = A2$T 
  C = t(Q1) %*% C %*% Q2
  
  Rsq = R1 * R1
  I = diag(dim(A)[1])
  
  
  XXX = solve(Q1) %*% XX %*% solve(t(Q2))
  
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
  
  sums = c(sums, sum(XXX -X ) < eps)
}




sylvester = function(A, B, C){
  
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
  return(Q1 %*% X %*% t(Q2))x
}
