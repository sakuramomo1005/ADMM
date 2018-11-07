## 2018 10 14
## solution of sylvester equation
## no optimazatin yet
## version 1

sylvester = function(A,B,C, eps =0.0001){
  
  library(Matrix)
  library(MASS)
  
  A1 = Schur(A)
  Q1 = A1$Q; R1 = A1$T
  
  A2 = Schur(B)
  Q2 = A2$Q; R2 = A2$T
  C = t(Q1) %*% C %*% Q2
  
  Rsq = R1 * R1
  I = diag(dim(A)[1])
  
  X[,1] = ginv(R1 + R2[1,1] * I) %*% C[,1]
  k = 2
  n = dim(R2)[1]
  
  while(k < n){
    if(R2[k+1, k] < eps){
      left = R1 + R2[k,k] * I
      right = C[,k]
      temp = matrix(0, dim(X)[1],1)
      for(i in 1:(k-1)){
        temp = temp + X[,i] * R2[i,k]
      }
      temp = matrix(temp, dim(C)[1],1)
      X[,k] = ginv(left) %*% (right + temp)
      k = k+1
    }else{
      r11 = R2[k,k]; r12 = R2[k, k+1]; r21 = R2[k+1, k]; r22 = R2[k+1, k+1]
      temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
      for(i in 1:(k-1)){
        temp2 = temp2 + X[,i] * R2[i,k]
        temp3 = temp3 + X[,i] * R2[i,k+1]}
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
  }
  
  print(k)
  if(k == n){
    left = R1 + R2[k,k] * I
    right = C[,k]
    matrix(0, dim(X)[1],1)
    for(i in 1:(k-1)){
      temp = temp + X[,i] * R2[i,k]
    }
    temp = matrix(temp, dim(C)[1],1)
    X[,k] = ginv(left) %*% (right + temp)
  } 
  if(k == n+1){print('done')}
  return(round(Q1 %*% X %*% t(Q2)))
}

# test 1: 
A=matrix(c(1,0,2,3,4,1,0,2,0,5,5,6,1,7,9,0),4,4,byrow = TRUE)
B=matrix(c(0,-1,1,0,3,4,5,1,2),3,3,byrow = TRUE)
X=matrix(c(1,0,2,0,0,3,1,1,2,3,1,6),4,3,byrow = TRUE)
XX= X
C = A%*% X + X%*%B

A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T

n = dim(R2)[1]
for(k in 1: (n-1)){
  print(R2[k+1, k])
}
res = sylvester(A,B,C)

abs(res - XX) < 0.001

# test 2: 

A=matrix(c(1,0,2,3,4,1,0,2,0,5,5,6,1,7,9,0),8,8,byrow = TRUE)
B=matrix(c(0,-1,1,0,3,4,5,1,2),4,4,byrow = TRUE)
X=matrix(rnorm(32,3,6),8,4,byrow = TRUE)
XX= X
C = A%*% X + X%*%B

res = sylvester(A,B,C)
 





library(Matrix)
library(MASS)


A=matrix(c(1,0,2,3,4,1,0,2,0,5,5,6,1,7,9,0),8,8,byrow = TRUE)
B=matrix(c(0,-1,1,0,3,4,5,1),4,4,byrow = TRUE)
X=matrix(rnorm(32,3,6),8,4,byrow = TRUE)
XX= X
C = A%*% X + X%*%B

A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T

A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T
C = t(Q1) %*% C %*% Q2

Rsq = R1 * R1
I = diag(dim(A)[1])

X = matrix(0, dim(A)[1],dim(B)[1])

#X[,1] = ginv(R1 + R2[1,1] * I) %*% C[,1]

n = dim(R2)[1]
k = n

while(k>0){
  if(k == 1){
    left = R1 + R2[k,k] * I
    right = C[,k]
    temp = matrix(0, dim(X)[1],1)
    for(i in (k+1):n){
      temp = temp + X[,i] * R2[i,k]
    }
    temp = matrix(temp, dim(C)[1],1)
    X[,k] = ginv(left) %*% (right + temp)
    k = k-1
  }else if(k>1 & R2[k, k-1] < eps){
    left = R1 + R2[k,k] * I
    right = C[,k]
    temp = matrix(0, dim(X)[1],1)
    for(i in 1:(k-1)){
      temp = temp + X[,i] * R2[i,k]
    }
    temp = matrix(temp, dim(C)[1],1)
    X[,k] = ginv(left) %*% (right + temp)
    k = k-1
  }else{
    r11 = R2[k-1,k-1]; r12 = R2[k-1, k]; 
    r21 = R2[k, k-1]; r22 = R2[k, k]
    temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
    for(i in 1:(k-1)){
      temp2 = temp2 + X[,i] * R2[i,k]
      temp3 = temp3 + X[,i] * R2[i,k+1]}
    temp2 = matrix(temp2, dim(X)[1],1)
    temp3 = matrix(temp3, dim(X)[1],1)
    b1 = C[,k-1] - temp2 
    b2 = C[,k] - temp3
    b1_prime = R1 %*% b1 + r22 * b1 - r21 * b2
    b2_prime = R1 %*% b2 + r11 * b2 - r12 * b1
    b_prime = matrix(0, dim(X)[1],2)
    b_prime[,1] = b1_prime; b_prime[,2] = b2_prime
    X[,(k-1):(k)] = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                         (r11*r22 - r12*r21) * I) %*% b_prime 
    k = k-2
  }
}

print(k)
if(k == n){
  left = R1 + R2[k,k] * I
  right = C[,k]
  matrix(0, dim(X)[1],1)
  for(i in 1:(k-1)){
    temp = temp + X[,i] * R2[i,k]
  }
  temp = matrix(temp, dim(C)[1],1)
  X[,k] = ginv(left) %*% (right + temp)
} 
if(k == n+1){print('done')}



tol = 10 * eps
if(norm(E,type='F' < tol)){
  X = matrix(0, dim(E)[1],dim(E)[2])
}

m = dim(E)[1]; n = dim(E)[2]
Y = matrix(0, m, n)

if(is.null(C) & is.null(D)){
  AEQ = TRUE
  C = B
}else{
  AEQ = FALSE
}

if(is.null(C)){
  A1 = Schur(A)$Q
  Z1 = A1$Q
  P = A1$T  
  Q1 = t(Z1)
  S = diag(m)
}
m= 3; n =2
A = matrix(rnorm(m*m,3,5),m,m)
C = A
B = matrix(rnorm(n*n,3,5),n,n)
D = matrix(rnorm(n*n,3,5),n,n)

X = matrix(rnorm(m*n,4,5),m,n)

A = matrix(c(0.2769  ,  0.0971,
     0.0462   , 0.8235),2,2,byrow = TRUE)
B = matrix(c(0.6948   , 0.0344  ,  0.7655,
     0.3171   , 0.4387   , 0.7952,
     0.9502   , 0.3816   , 0.1869),3,3,byrow=TRUE)
C =  matrix(c(0.4898  ,  0.6463,
      0.4456   , 0.7094),2,2,byrow=TRUE)
D = matrix(c(0.7547   , 0.6551  ,  0.4984,
     0.2760    ,0.1626   , 0.9597,
     0.6797    ,0.1190   , 0.3404),3,3,byrow = TRUE)
X = matrix(c(0.5853  ,  0.7513  ,  0.5060,
       0.2238   , 0.2551 ,   0.6991),2,3,byrow=TRUE)


qzbd = qz(D,B)
T = qzbd$S
R = qzbd$T
Q2 = qzbd$Q
Z2 = qzbd$Z

E = A %*% X %*% t(B) + C %*% X %*% t(D) 

XX = X 

A1 = Schur(A)
Z1 = A1$Q
P = A1$T
Q1 = t(Z1)
S = diag(dim(A)[1])

A2 = Schur(D)
Z2 = A2$Q
T = A2$T
Q2 = t(Z2)
R = diag(dim(D)[1])

E = A %*% X +  X %*% t(D) 


F = Q1 %*% E %*% t(Q2) 

m = dim(E)[1];n = dim(E)[2]
PY = matrix(0, m, n)
SY = matrix(0, m, n)

k = n
Y = matrix(0, m, n)


F = matrix(c(   -0.6144,   -1.8989 ,   2.4252,
                0.0834  ,  0.1691 ,  -0.3417),m,n,byrow = TRUE)

T = matrix(c( -0.3825 ,   0.3597  , -0.2319,
              0   , 0.6165  , -0.8465,
              0  ,       0 ,  1.1986),n,n,byrow = TRUE)
R = matrix(c(0.3596   , 0.0144  , -0.2347,
             0    ,0.8534  , -0.9041,
             0     ,    0 ,   1.1612),n,n,byrow=TRUE)

Q2 = matrix(c( 0.8643  , -0.3962  , -0.3100,
               0.3277 ,  -0.0241   , 0.9445,
               -0.3816 ,  -0.9179 ,   0.1090),n,n,byrow=TRUE)

Z2 = matrix(c(0.2109  ,  0.8964  , -0.3899,
              -0.9431  ,  0.0817 ,  -0.3223,
              0.2570   ,-0.4356 ,  -0.8626),n,n,byrow=TRUE)

k=n
while(k > 0){
  if(k ==1){
    jj = (k+1):n
    rhs = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),1) - 
      SY[,jj] %*% matrix(t(T[k,jj]),length(jj),1)
    tmp = (P + (T[k,k]/R[k,k]) * S)
    rhs = rhs/R[k,k]
    Y[,k] = solve(tmp,rhs)
    
    PY[,k] = P %*% Y[,k]
    SY[,k] = S %*% Y[,k]
    k = k-1
  }else if((k< n) & (T[k, k-1] < tol)){
    jj = (k+1):n
    rhs = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),1) - 
      SY[,jj] %*% matrix(t(T[k,jj]),length(jj),1)
    tmp = (P + (T[k,k]/R[k,k]) * S)
    rhs = rhs/R[k,k]
    Y[,k] = solve(tmp,rhs)
    
    PY[,k] = P %*% Y[,k]
    SY[,k] = S %*% Y[,k]
    k = k-1
  }else if((k == n) & (T[k, k-1] < tol)){
    rhs = F[,k]
    tmp = (P + (T[k,k]/R[k,k]) * S)
    rhs = rhs/R[k,k]
    Y[,k] = solve(tmp,rhs)
    
    PY[,k] = P %*% Y[,k]
    SY[,k] = S %*% Y[,k]
    k = k-1
  }else if((k < n) & (T[k, k-1] > tol)){
    jj = (k+1):n
    rhs1 = F[,k-1] - PY[,jj] %*% matrix(t(R[k-1,jj]),length(jj),k-1) - 
      SY[,jj] %*% matrix(t(T[k-1,jj]),length(jj),k-1)
    rhs2 = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),k) - 
      SY[,jj] %*% matrix(t(T[k,jj]),length(jj),k)
    SM1 = cbind(R[k-1,k-1] * P + T[k-1, k-1] * S, R[k-1, k-1] * P + T[k-1,k]*S)
    SM2 = cbind(T[k,k-1] * S, R[k,k] * P + T[k,k] * S)
    SM = rbind(SM1,SM2)    
    idx = c(matrix(1:(2*m),byrow = TRUE,2,m))
    
    rhs = cbind(rhs1, rhs2)
    rhs = matrix(rhs, dim(rhs1)[1],2)
    UM = solve(SM[idx,idx],rhs[idx])
    UM[idx] = UM
    
    Y[,(k-1):k] = matrix(UM,m,2)
    PY[,(k-1):k] = P %*% Y[,(k-1):k]
    SY[,(k-1):k] = S %*% Y[,(k-1):k]
    k = k-2
  }else{
    rhs1 = F[,k-1]
    rhs2 = F[,k] 
    SM1 = cbind(R[k-1,k-1] * P + T[k-1, k-1] * S, R[k-1, k-1] * P + T[k-1,k]*S)
    SM2 = cbind(T[k,k-1] * S, R[k,k] * P + T[k,k] * S)
    SM = rbind(SM1,SM2)    
    idx = c(matrix(1:(2*m),byrow = TRUE,2,m))
    
    rhs = cbind(rhs1, rhs2)
    rhs = matrix(rhs, dim(rhs1)[1],2)
    UM = solve(SM[idx,idx],rhs[idx])
    UM[idx] = UM
    
    Y[,(k-1):k] = matrix(UM,m,2)
    PY[,(k-1):k] = P %*% Y[,(k-1):k]
    SY[,(k-1):k] = S %*% Y[,(k-1):k]
    k = k-2
  }
}

X = Z1 %*% Y %*% t(Z2)



test = function(A,B,C){
  E = C
  D = t(B)
  A1 = Schur(A)
  Z1 = A1$Q
  P = A1$T
  Q1 = t(Z1)
  S = diag(dim(A)[1])
  
  A2 = Schur(D)
  Z2 = A2$Q
  T = A2$T
  Q2 = t(Z2)
  R = diag(dim(D)[1])
  
  F = Q1 %*% E %*% t(Q2) 
  
  m = dim(E)[1];n = dim(E)[2]
  PY = matrix(0, m, n)
  SY = matrix(0, m, n)
  
  k = n
  Y = matrix(0, m, n)
  
  while(k > 0){
    if(k ==1){
      jj = (k+1):n
      rhs = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),1) - 
        SY[,jj] %*% matrix(t(T[k,jj]),length(jj),1)
      tmp = (P + (T[k,k]/R[k,k]) * S)
      rhs = rhs/R[k,k]
      Y[,k] = solve(tmp,rhs)
      
      PY[,k] = P %*% Y[,k]
      SY[,k] = S %*% Y[,k]
      k = k-1
    }else if((k< n) & (T[k, k-1] < tol)){
      jj = (k+1):n
      rhs = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),1) - 
        SY[,jj] %*% matrix(t(T[k,jj]),length(jj),1)
      tmp = (P + (T[k,k]/R[k,k]) * S)
      rhs = rhs/R[k,k]
      Y[,k] = solve(tmp,rhs)
      
      PY[,k] = P %*% Y[,k]
      SY[,k] = S %*% Y[,k]
      k = k-1
    }else if((k == n) & (T[k, k-1] < tol)){
      rhs = F[,k]
      tmp = (P + (T[k,k]/R[k,k]) * S)
      rhs = rhs/R[k,k]
      Y[,k] = solve(tmp,rhs)
      
      PY[,k] = P %*% Y[,k]
      SY[,k] = S %*% Y[,k]
      k = k-1
    }else if((k < n) & (T[k, k-1] > tol)){
      jj = (k+1):n
      rhs1 = F[,k-1] - PY[,jj] %*% matrix(t(R[k-1,jj]),length(jj),k-1) - 
        SY[,jj] %*% matrix(t(T[k-1,jj]),length(jj),k-1)
      rhs2 = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),k) - 
        SY[,jj] %*% matrix(t(T[k,jj]),length(jj),k)
      SM1 = cbind(R[k-1,k-1] * P + T[k-1, k-1] * S, R[k-1, k-1] * P + T[k-1,k]*S)
      SM2 = cbind(T[k,k-1] * S, R[k,k] * P + T[k,k] * S)
      SM = rbind(SM1,SM2)    
      idx = c(matrix(1:(2*m),byrow = TRUE,2,m))
      
      rhs = cbind(rhs1, rhs2)
      #rhs = matrix(rhs, dim(rhs1)[1],2)
      UM = solve(SM[idx,idx],rhs[idx])
      UM[idx] = UM
      
      Y[,(k-1):k] = matrix(UM,m,2)
      PY[,(k-1):k] = P %*% Y[,(k-1):k]
      SY[,(k-1):k] = S %*% Y[,(k-1):k]
      k = k-2
    }else{
      rhs1 = F[,k-1]
      rhs2 = F[,k] 
      SM1 = cbind(R[k-1,k-1] * P + T[k-1, k-1] * S, R[k-1, k-1] * P + T[k-1,k]*S)
      SM2 = cbind(T[k,k-1] * S, R[k,k] * P + T[k,k] * S)
      SM = rbind(SM1,SM2)    
      idx = c(matrix(1:(2*m),byrow = TRUE,2,m))
      
      rhs = cbind(rhs1, rhs2)
      rhs = matrix(rhs, dim(rhs1)[1],2)
      UM = solve(SM[idx,idx],rhs[idx])
      UM[idx] = UM
      
      Y[,(k-1):k] = matrix(UM,m,2)
      PY[,(k-1):k] = P %*% Y[,(k-1):k]
      SY[,(k-1):k] = S %*% Y[,(k-1):k]
      k = k-2
    }
  }
  
  X = Z1 %*% Y %*% t(Z2)
  print(X)
  #print(XX)
}

B = matrix(rnorm(16,0,6),4,4)
A = matrix(rnorm(9,0,6),3,3)
X = matrix(rnorm(12,0,7),3,4)
C = A%*% X + X%*% B
XX = X
XX
test(A,B,C)


A = matrix(c(1:5,5:1,1:4,5,0:4,5:1),5,5,byrow = TRUE)
D = matrix(c(1:4,4:1,2:5,5:2),4,4,byrow = TRUE)
X = matrix(c(1,2,0,1,1,2,1,1,2,2,1,1,2,0,1,1,1,1,1,1),5,4)

E =  A %*% X + X %*% t(D);

A1 = Schur(A)
Z1 = A1$Q
P = A1$T
Q1 = t(Z1)
S = diag(dim(A)[1])

A2 = Schur(D)
Z2 = A2$Q
T = A2$T
Q2 = t(Z2)
R = diag(dim(D)[1])

F = Q1 %*% E %*% t(Q2) 

m = dim(E)[1];n = dim(E)[2]
PY = matrix(0, m, n)
SY = matrix(0, m, n)

k = n
Y = matrix(0, m, n)


while(k > 0){
  print(k)
  k0 = k
  mark = 0
  if(k0 ==1){
    jj = (k+1):n
    rhs = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),1) - 
      SY[,jj] %*% matrix(t(T[k,jj]),length(jj),1)
    tmp = (P + (T[k,k]/R[k,k]) * S)
    rhs = rhs/R[k,k]
    Y[,k] = solve(tmp) %*% rhs
    
    PY[,k] = P %*% Y[,k]
    SY[,k] = S %*% Y[,k]
    k = k-1
    mark = 1
  }
  if(mark == 0){
    if(k0 <n & T[k0, k0-1] < tol){
      jj = (k+1):n
      rhs = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),1) - 
        SY[,jj] %*% matrix(t(T[k,jj]),length(jj),1)
      tmp = (P + (T[k,k]/R[k,k]) * S)
      rhs = rhs/R[k,k]
      Y[,k] = solve(tmp) %*% rhs
      
      PY[,k] = P %*% Y[,k]
      SY[,k] = S %*% Y[,k]
      k = k-1
      mark = 1
    }
  }
  if(mark == 0){
    if(k0==n & T[k0, k0-1] < tol){
      rhs = F[,k]
      tmp = (P + (T[k,k]/R[k,k]) * S)
      rhs = rhs/R[k,k]
      Y[,k] = solve(tmp) %*% rhs
      
      PY[,k] = P %*% Y[,k]
      SY[,k] = S %*% Y[,k]
      k = k-1
      mark = 1
    }
  }
  if(mark == 0){
    if(k0 <n & T[k0, k0-1] >= tol){
      mark = 1
      jj = (k+1):n
      rhs1 = F[,k-1] - PY[,jj] %*% matrix(t(R[k-1,jj]),length(jj),k-1) - 
        SY[,jj] %*% matrix(t(T[k-1,jj]),length(jj),k-1)
      rhs2 = F[,k] - PY[,jj] %*% matrix(t(R[k,jj]),length(jj),k) - 
        SY[,jj] %*% matrix(t(T[k,jj]),length(jj),k)
      SM1 = cbind(R[k-1,k-1] * P + T[k-1, k-1] * S, R[k-1, k-1] * P + T[k-1,k]*S)
      SM2 = cbind(T[k,k-1] * S, R[k,k] * P + T[k,k] * S)
      SM = rbind(SM1,SM2)    
      idx = c(matrix(1:(2*m),byrow = TRUE,2,m))
      
      rhs = cbind(rhs1, rhs2)
      #rhs = matrix(rhs, dim(rhs1)[1],2)
      UM = solve(SM[idx,idx]) %*% rhs[idx]
      UM[idx] = UM
      
      Y[,(k-1):k] = matrix(UM,m,2)
      PY[,(k-1):k] = P %*% Y[,(k-1):k]
      SY[,(k-1):k] = S %*% Y[,(k-1):k]
      k = k-2
    }
  }
  if(mark == 0){
    if(k0==n & T[k0, k0-1] >= tol){
      mark =1
      rhs1 = F[,k-1]
      rhs2 = F[,k] 
      SM1 = cbind(R[k-1,k-1] * P + T[k-1, k-1] * S, R[k-1, k-1] * P + T[k-1,k]*S)
      SM2 = cbind(T[k,k-1] * S, R[k,k] * P + T[k,k] * S)
      SM = rbind(SM1,SM2)    
      idx = c(matrix(1:(2*m),byrow = TRUE,2,m))
      
      rhs = cbind(rhs1, rhs2)
      rhs = matrix(rhs, dim(rhs1)[1],2)
      UM = solve(SM[idx,idx]) %*% rhs[idx]
      UM[idx] = UM
      
      Y[,(k-1):k] = matrix(UM,m,2)
      PY[,(k-1):k] = P %*% Y[,(k-1):k]
      SY[,(k-1):k] = S %*% Y[,(k-1):k]
      k = k-2
    }
  }
}


Z1 %*% Y %*% t(Z2)

solve(Z1) %*% XX %*% solve(t(Z2))


Z1 %*% solve(Z1) %*% XX %*% solve(t(Z2)) %*% t(Z2)