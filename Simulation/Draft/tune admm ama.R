library(fossil)
library(cvxbiclustr)
library(cvxclustr)
 

dist_weight = function (X, phi, dist.type, p){
  dist_X <- as.numeric(dist(t(X), method = dist.type))
  exp(-phi * dist_X^p)
}
elk = function(n,p){
  # n,l
  count = 0
  el1 = el2 = matrix(0,n,n*(n-1)/2)
  
  for(i in (n-1):1){
    temp = matrix(0,n,i)
    temp[n-i,] = 1
    el1[,(count+1):(count+i)] = temp
    el2[(n-i+1):n,(count+1):(count+i)] = diag(1,i,i)
    count = count +i
  }  
  
  # p,k
  count = 0
  ek1 = ek2 = matrix(0,p,p*(p-1)/2)
  for(i in (p-1):1){
    temp = matrix(0,p,i)
    temp[p-i,] = 1
    ek1[,(count+1):(count+i)] = temp
    ek2[(p-i+1):p,(count+1):(count+i)] = diag(1,i,i)
    count = count +i
  }
  return(list(el1 = el1, el2 = el2, ek1 = ek1, ek2 = ek2))
}
L_num = function(n){
  t = matrix(0, n*(n-1)/2,2)
  count = 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      count = count + 1
      t[count,1]= i; t[count,2]=j
    }
  }
  LL = data.frame(l1 = t[,1],l2 = t[,2])
  LL = LL[order(LL$l1),]
  rownames(LL) = NULL
  return(LL)
}
sylvester= function(A, B, C, eps = 0.0001){
  
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
          temp = (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
        }
        temp = matrix(temp, dim(C)[1],1)
        X[,k] = ginv(left) %*% (right - temp)
        # mytry = myTryCatch(solve(left))
        #if(is.null(mytry$error) == 0){er = c(er,tt)}
        k = k+1
      }else{
        r11 = R2[k,k]; r12 = R2[k, k+1]; r21 = R2[k+1, k]; r22 = R2[k+1, k+1]
        temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        if(k == 1){
          temp2 = temp2
          temp3 = temp3
        }else{
          temps = X[,1:(k-1)] %*% matrix(R2[1:(k-1),k:(k+1)],k-1,2)
          temp2 = temps[,1]
          temp3 = temps[,2]
        }
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
          temp = (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
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
        
        if(k == 1){
          temp2 = temp2
          temp3 = temp3
        }else{
          temps = X[,1:(k-1)] %*% matrix(R22[1:(k-1),k:(k+1)],k-1,2)
          temp2 = temps[,1]
          temp3 = temps[,2]
        }
        
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

# AMA
bi_ama_order = function(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 1000, tol = 1e-5){
  
  A = X
  
  n=dim(X)[1]; p=dim(X)[2]
  
  n2 = n*(n-1)/2
  p2 = p*(p-1)/2
  
  #eplison_p=L_num(p)
  #eplison_n=L_num(n)
  
  lambda_1 = matrix(1,p,n2)
  lambda_2 = matrix(1,n,p2)
  
  w <- dist_weight(t(X)/sqrt(n),phi, dist.type = "euclidean", p = 2 )
  w_l <- knn_weights(w,k = kk,p)
  u <- dist_weight((X)/sqrt(p),phi, dist.type = "euclidean", p = 2 )
  u_k <- knn_weights(u,k = kk,n)
  
  elks = elk(n,p)
  el1 = elks$el1
  el2 = elks$el2
  ek1 = elks$ek1
  ek2 = elks$ek2
  
  for(iter in 1:niter){
    print('iter')
    print(iter)
    
    A_old = A; lambda_1_old = lambda_1; lambda_2_old = lambda_2
    
    lv = lambda_1
    lz = lambda_2
    
    C2 = (el1-el2) %*% t(lv)
    C3 = lz %*% t(ek1-ek2)
    A = X +  C2 + C3  
    
    al1 = t(A) %*% el1; al2 = t(A) %*% el2
    ak1 = A %*% ek1; ak2 = A %*% ek2
    
    lambda1 = lambda_1 - nu1 * (al1 - al2)
    lambda2 = lambda_2 - nu2 * (ak1 - ak2)
    
    for(i in 1:dim(lambda_1)[2]){
      tau = gamma_1 * w_l[i]
      norm_x = sum(lambda1[,i]^2)
      if(norm_x > tau){
        lambda_1[,i] = tau/norm_x * lambda_1[,i]
      }else{
        lambda_1[,i] = lambda_1[,i]
      }
    }
    
    for(i in 1:dim(lambda_2)[2]){
      tau = gamma_2 * u_k[i]
      norm_y = sum(lambda2[,i]^2)
      if(norm_y > tau){
        lambda_2[,i] = tau/norm_y * lambda_2[,i]
      }else{
        lambda_2[,i] = lambda_2[,i]
      }
    }
    
    print(paste('A',sum(abs(A - A_old))))
    print(paste('l',sum(abs(lambda_1 - lambda_1_old))))
    print(paste('2',sum(abs(lambda_2 - lambda_2_old))))
    
    if(sum(abs(A - A_old)) < tol & 
       sum(abs(lambda_1 - lambda_1_old)) < tol & 
       sum(abs(lambda_2 - lambda_2_old)) < tol){
      return(list(A = A, 
                  lambad_1 = lambda1, lambad_2 = lambda2, niter = iter))
      break
    }
    
  }
  
  if(iter == niter){
    print(paste('not converge within',iter, 'times'))
    return(list(A = A, 
                lambad_1 = lambda1, lambad_2 = lambda2, niter = iter))
  }
}
# ADMA
bi_ADMM_order = function(X, nu1, nu2, gamma_1, gamma_2, kk=5, niter = 1000,tol = 0.1){
  
  library(cvxbiclustr)
  library(MASS)
  library(Matrix)
  
  n=dim(X)[1]; p=dim(X)[2]
  
  n2 = n*(n-1)/2
  p2 = p*(p-1)/2
  
  elks = elk(n,p)
  el1 = elks$el1
  el2 = elks$el2
  ek1 = elks$ek1
  ek2 = elks$ek2
  
  w <- dist_weight(t(X)/sqrt(n),phi, dist.type = "euclidean", p = 2 )
  w_l <- knn_weights(w,k = kk,p)
  
  u <- dist_weight((X)/sqrt(p),phi, dist.type = "euclidean", p = 2 )
  u_k <- knn_weights(u,k = kk,n)
  
  A = matrix(0,n,p)
  v = matrix(0,p,n2)
  z = matrix(0,n,p2)
  lambda_1 = matrix(0,p,n2)
  lambda_2 = matrix(0,n,p2)
  
  
  for(iter in 1: niter){
    print('iter')
    print(iter)
    
    A_old = A; v_old = v; z_old = z; lambda_1_old = lambda_1; lambda_2_old = lambda_2
    
    # update A
    
    En = diag(0:(n - 1)) + diag((n - 1):0) - matrix(1, n, n) + diag(1, n, n)
    Ep = diag(0:(p - 1)) + diag((p - 1):0) - matrix(1, p, p) + diag(1, p, p)
    
    M = diag(1,n,n) + nu1 * En
    
    N = nu2 * Ep
    
    lv = lambda_1+ nu1 * v
    lz = lambda_2 + nu2 * z
    
    C2 = (el1-el2) %*% t(lv)
    C3 = lz %*% t(ek1-ek2)
    C = X +  C2 + C3  
    
    A = sylvester(M,t(N),C)
    
    al1 = t(A) %*% el1; al2 = t(A) %*% el2
    ak1 = A %*% ek1; ak2 = A %*% ek2
    
    
    # update vz
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
    
    # update lambda
    lambda_1 = lambda_1 + nu1 * (v - al1 + al2)
    
    # update lambda 2
    lambda_2 = lambda_2 + nu2 * (z - ak1 + ak2)
    
    print(paste('A',sum(abs(A - A_old))))
    print(paste('v',sum(abs(v - v_old))))
    print(paste('z',sum(abs(z -z_old))))
    print(paste('l',sum(abs(lambda_1 - lambda_1_old))))
    print(paste('2',sum(abs(lambda_2 - lambda_2_old))))
    
    # whether coverage
    if(sum(abs(A - A_old)) < tol & 
       sum(abs(v - v_old)) < 0.01& 
       sum(abs(z - z_old)) < 0.01 & 
       sum(abs(lambda_1 - lambda_1_old)) < 0.01 & 
       sum(abs(lambda_2 - lambda_2_old)) <0.01){
      return(list(A = A, v = v, z = z, 
                  lambad_1 = lambda_1, lambad_2 = lambda_2, niter = iter))
      break
    }
  }
  
  if(iter == niter){
    print(paste('not converge within',iter, 'times'))
    return(list(A = A, v = v, z = z, 
                lambad_1 = lambda_1, lambad_2 = lambda_2, niter = iter))
  }
}

# data generation
data_gen=function(seed){
  set.seed(seed)
  n = 80
  p = 40
  mu = seq(-6,6,0.5)
  theta = 1.5
  x = matrix(0, n, p)
  
  # 2 row groups
  # 8 column groups
  row_group = 4
  col_group = 8
  
  mu_rc = matrix(0, row_group,col_group)
  for(r in 1:row_group){
    for(c in 1:col_group){
      mu_rc[r,c] = sample(mu)[1]
    }
  }
  
  row_assign = c()
  col_assign = c()
  for(i in 1:n){
    row_assign = c(row_assign, sample(1:2)[1])
  }
  for(i in 1:p){
    col_assign = c(col_assign, sample(1:8)[1])
  }
  
  for(i in 1:n){
    for(j in 1:p){
      r = row_assign[i]
      c = col_assign[j]
      mu_temp = mu_rc[r,c]
      x[i,j] = rnorm(1,mu_temp,theta)
    }
  }
  
  rownames(x) = paste('row',row_assign,sep='')
  colnames(x) = paste('col',col_assign,sep='')
  
  return(x)
}


setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 1116/1')

# save results
nu1 = nu2=3
X = lung
X <- X - mean(X)
X <- X/norm(X,'f')

for(gamma_1 in c(0,0.1,0.3,0.5,0.7)){
  for(gamma_2 in c(0,0.1,0.3,0.5,0.7)){
    phi = 0.5;kk=5;
    res1 = bi_ADMM_order(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 10000, tol = 1e-3)
    res2 = bi_ama_order(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 1000, tol = 1e-5)
    ind1 = rand.index(X,res1$A)
    ind2 = rand.index(X,res2$A)
    result = list(res1 = res1, res2 = res2, ind1 = ind1, ind2 = ind2)
    names = paste('lungg1_',gamma_1,'g2_',gamma_2,'.RData',sep='')
    save(result,file = names)
  }
}


# plots
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)


pdf('plots_ama.pdf')
ind1 = c()
ind2= c()
for(gamma_1 in c(1,10,100,1000)){
  for(gamma_2 in c(1,10,100,1000)){
    names = paste('g1_',gamma_1,'g2_',gamma_2,'.RData',sep='')
    load(names)
    ind1 = c(ind1, result$ind1)
    #heatmap(result$res1$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'admm'))
    heatmap(result$res2$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'ama'))
  }
}
dev.off()

pdf('plots_admm.pdf')
ind1 = c()
for(gamma_1 in c(1,10,100,1000)){
  for(gamma_2 in c(1,10,100,1000)){
    names = paste('g1_',gamma_1,'g2_',gamma_2,'.RData',sep='')
    load(names)
    ind1 = c(ind1, result$ind1)
    heatmap(result$res1$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'admm'))
    #heatmap(result$res2$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'ama'))
  }
}
dev.off()

pdf('plots_lung_ama.pdf')
ind2 = c()
for(gamma_1 in c(1,10,100,1000)){
  for(gamma_2 in c(1,10,100,1000)){
    names = paste('lungg1_',gamma_1,'g2_',gamma_2,'.RData',sep='')
    load(names)
    ind2 = c(ind2, result$ind2)
    #heatmap(result$res1$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'admm'))
    heatmap(result$res2$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'ama'))
  }
}
dev.off()

pdf('a.pdf')
ind1 = c()
ind2= c()
for(gamma_1 in c(1,10,10,100)){
  for(gamma_2 in c(1,10,100,1000)){
    names = paste('lungg1_',gamma_1,'g2_',gamma_2,'.RData',sep='')
    load(names)
    ind1 = c(ind1, result$ind1)
    ind2 = c(ind2, result$ind2)
    heatmap(result$res1$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'admm'))
    #heatmap(result$res2$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'ama'))
  }
}
dev.off()



pdf('b.pdf')
ind1 = c()
ind2= c()
for(gamma_1 in c(0,0.1,0.3,0.5,0.7)){
  for(gamma_2 in c(0,0.1,0.3,0.5,0.7)){
    names = paste('lungg1_',gamma_1,'g2_',gamma_2,'.RData',sep='')
    load(names)
    ind1 = c(ind1, result$ind1)
    ind2 = c(ind2, result$ind2)
    heatmap(result$res1$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'admm'))
    heatmap(result$res2$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'ama'))
  }
}
dev.off()

