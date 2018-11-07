# 20181002
# solving sylvester equation

# AX+XB = C
# A m*m, x m*n, B n*n, c m*n
library(MASS)
library(Matrix)

#A=matrix(16,2,2)
#B=matrix(9,2,2)
#X=matrix(12,2,2)
#XX = X
#C = A%*% X + X%*%B

A=matrix(c(1,0,2,3,4,1,0,2,0,5,5,6,1,7,9,0),4,4,byrow = TRUE)
B=matrix(c(0,-1,1,0,3,4,5,1,2),3,3,byrow = TRUE)
X=matrix(c(1,0,2,0,0,3,1,1,2,3,1,6),4,3,byrow = TRUE)
XX= X
C = A%*% X + X%*%B

sylvester=function(A,B,C){
  I_A=diag(1,dim(A)[1],dim(A)[2])
  I_B=diag(1,dim(B)[1],dim(B)[2])
  M=kronecker(I_B,A) + kronecker(t(B),I_A)
  L=matrix(c(C),dim(A)[1]*dim(B)[1],1)
  cx=ginv(M) %*% L
  X=matrix(cx,dim(A)[1],dim(B)[1])
  return(X)
}

X=sylvester(A,B,C)
round(X);XX

round(A%*%X + X%*%B)


#A, B, C 

#A = B =matrix(1:9, 3,3)
#X = matrix(9:1 , 3, 3)
#C = A%*% X + X%*%B

A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T

A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T

C = -C
C = t(Q1) %*% C %*% Q2
Rsq = R1 * R1
I = diag(dim(A)[1])

b = - C[,1]
X = matrix(0, dim(A)[1], dim(B)[1])
X[,1] = (solve(R1 +  R2[1,1] * I) %*% b)

n = dim(B)[1]
j=2
while(j < n){
  if(j<n & abs(R2[j+1,j]<0.01)){
    left = R1 + R2[j,j] * I
    if(j == 2){
      t1 = matrix(X[, 1:(j-1)], dim(X)[1],1:(j-1))
      t2 = matrix(R2[1: (j-1), j], 1:(j-1),j)
      right = - C[,j] - t1 %*% t2
    }else{
      t1 = matrix(X[, 1:(j-1)], dim(X)[1],1:(j-1))
      t2 = matrix(R2[1: (j-1), j], 1:(j-1),j)
      right = - C[,j] - t1 %*% t2
    }
    j = j+1
    X[,j] = qr.solve(left) %*% right
  }else{
    r11 = R2[j,j];r12 = R2[j, j+1]
    r21 = R2[j+1,j]; r22 = R2[j+1, j+1]
    t1 = matrix(X[, 1:(j-1)], dim(X)[1],1:(j-1))
    t2 = matrix(R2[1: (j-1), j], 1:(j-1),2)
    b = -C[,j:(j+1)] - t1 %*% t2
    
    b = rbind(R1 %*% b[,1] + r22 * b[,1] - r21 * b[,2], R1 %*% b[,2] + r11 * b[,2] - r12 * b[,1])
    X[, j:(j+1)] = qr.solve(Rsq + (r11 + r22) * R1 + (r11 * r22 - r12 * r21) * I) %*% b
    j = j+2
  }
}




A = matrix( 1:36, 12,12)
B = matrix( 49:1, 14, 14)
X = matrix(1:42, 12, 14)

#A=matrix(c(1,0,2,3,4,1,0,2,0,5,5,6,1,7,9,0),4,4,byrow = TRUE)
#B=matrix(c(0,-1,1,0,3,4,5,1,2),3,3,byrow = TRUE)
#X=matrix(c(1,0,2,0,0,3,1,1,2,3,1,6),4,3,byrow = TRUE)
XX= X
C = A%*% X + X%*%B

A1 = Schur(A)
Q1 = A1$Q; R1 = A1$T

A2 = Schur(B)
Q2 = A2$Q; R2 = A2$T

C = -C
C = t(Q1) %*% C %*% Q2
Rsq = R1 * R1
I = diag(dim(A)[1])
j = 1
X = matrix(0, dim(A)[1], dim(B)[1])
m = dim(A)[1]; n = dim(B)[1]
library(MASS)
while(j<(n+1)){
  if(j < n){
    if(R2[j+1,j]<eps){
      if(j == 1){
        b = -C[, j]
      }else{
        b = -C[, j] - matrix(X[, 1:(j-1)], dim(X)[1],j-1) %*% matrix(R2[1:(j-1),j], j-1,1)
      }
      x = ginv(R1 + R2[j,j] * I) %*% b
      X[, j] = x
      j = j+1
    
    }else{
      r11 = R2[j,j];r12 = R2[j,j+1];r21 = R2[j+1, j]; r22 = R2[j+1, j+1]
      b = -C[,j:(j+1)] - matrix(X[, 1:(j-1)], dim(X)[1],j-1) %*% matrix(R2[1:(j-1),j:(j+1)],j-1,2)
      b = cbind(R1 %*% b[,1] + r22 * b[,1] - r21 * b[,2], R1 %*% b[,2] + r11 * b[,2] - r12 * b[,1])
      x = ginv(Rsq + (r11+r22)*R1 + (r11* r22 - r12 * r21) * I) %*% b
      X[,j:(j+1)] = x
      j = j+2
      }
  }else{
    
    if(j < n){
      r11 = R2[j,j];r12 = R2[j,j+1];r21 = R2[j+1, j]; r22 = R2[j+1, j+1]
    }else{
      r11 = R2[j,j];r12 = R2[j,1];r21 = R2[1, j]; r22 = R2[1, 1]
    }
    if(j == 1){b = -C[,j:(j+1)]
    }else if(j == n){
      b = -C[,c(j,1)] - matrix(X[, 1:(j-1)], dim(X)[1],j-1) %*% matrix(R2[1:(j-1),c(j,1)],j-1,2)
    }else{
      b = -C[,j:(j+1)] - matrix(X[, 1:(j-1)], dim(X)[1],j-1) %*% matrix(R2[1:(j-1),j:(j+1)],j-1,2)
    }
    b = cbind(R1 %*% b[,1] + r22 * b[,1] - r21 * b[,2], R1 %*% b[,2] + r11 * b[,2] - r12 * b[,1])
    x = ginv(Rsq + (r11+r22)*R1 + (r11* r22 - r12 * r21) * I) %*% b
    if(j == n){
      X[,c(j,1)] = x
    }else{
      X[,j:(j+1)] = x
    }
    j = j+2
  }
}
round(Q1 %*% X %*% t(Q2))
XX


j =1
while(j < (n+1)){
  if(j < n){
    if(R2[j+1,j]<eps){
      if(j ==1){
        b = -C[, j]
      }else{
        b = -C[, j] - matrix(X[, 1:(j-1)], dim(X)[1],j-1) %*% matrix(R2[1:(j-1),j], j-1,1)
      }
      x = ginv(R1 + R2[j,j] * I) %*% b
      X[, j] = x
      j = j+1
    }else{
      r11 = R2[j,j];r12 = R2[j,j+1];r21 = R2[j+1, j]; r22 = R2[j+1, j+1]
      if(j ==1){
        b = -C[,j:(j+1)]
      }else{
        b = -C[,j:(j+1)] - matrix(X[, 1:(j-1)], dim(X)[1],j-1) %*% matrix(R2[1:(j-1),j:(j+1)],j-1,2)
      }
      b = cbind(R1 %*% b[,1] + r22 * b[,1] - r21 * b[,2], R1 %*% b[,2] + r11 * b[,2] - r12 * b[,1])
      x = ginv(Rsq + (r11+r22)*R1 + (r11* r22 - r12 * r21) * I) %*% b
      X[,j:(j+1)] = x
      j = j+2
    }
  }else{
    r11 = R2[j,j];r12 = R2[j,1];r21 = R2[1, j]; r22 = R2[1, 1]
    b = -C[,c(j,1)] - matrix(X[, 1:(j-1)], dim(X)[1],j-1) %*% matrix(R2[1:(j-1),c(j,1)],j-1,2)
    b = cbind(R1 %*% b[,1] + r22 * b[,1] - r21 * b[,2], R1 %*% b[,2] + r11 * b[,2] - r12 * b[,1])
    x = ginv(Rsq + (r11+r22)*R1 + (r11* r22 - r12 * r21) * I) %*% b
    X[,c(j)] = x[,1]
    j = j+2
  }
}

round(Q1 %*% X %*% t(Q2))
XX



library(R.matlab)

library(Matrix)
library(MASS)

A=matrix(c(1,0,2,3,4,1,0,2,0,5,5,6,1,7,9,0),4,4,byrow = TRUE)
B=matrix(c(0,-1,1,0,3,4,5,1,2),3,3,byrow = TRUE)
X=matrix(c(1,0,2,0,0,3,1,1,2,3,1,6),4,3,byrow = TRUE)
XX= X
C = A%*% X + X%*%B

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
    matrix(0, dim(X)[1],1)
    for(i in 1:(k-1)){
      temp = temp + X[,i] * R2[i,k]
    }
    temp = matrix(temp, dim(C)[1],1)
    X[,k] = ginv(left) %*% (rigth + temp)
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
  X[,k] = ginv(left) %*% (rigth + temp)
} 
if(k == n+1){print('done')}
round(Q1 %*% X %*% t(Q2))
XX


Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)