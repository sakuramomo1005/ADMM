### 11/15
### ama with corrected weigth

library(cvxbiclustr)

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

bi_ama_order = function(X, nu1, nu2, gamma_1, gamma_2, 
                        kk = 5,phi = 0.5, niter = 1000, tol = 1e-5,
                        output = 1){
  
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
    if(output == 1){
      print('iter')
      print(iter)
      
      print(paste('A',sum(abs(A - A_old))))
      print(paste('l',sum(abs(lambda_1 - lambda_1_old))))
      print(paste('2',sum(abs(lambda_2 - lambda_2_old))))
      
    }
     
    if(sum(abs(A - A_old)) < tol & 
       sum(abs(lambda_1 - lambda_1_old)) < tol & 
       sum(abs(lambda_2 - lambda_2_old)) < tol){
      return(list(A = A, 
                  lambad_1 = lambda1, lambad_2 = lambda2, niter = iter))
      break
    }
    
  }
  
  if(iter == niter & output ==1){
    print(paste('not converge within',iter, 'times'))
    return(list(A = A, 
                lambad_1 = lambda1, lambad_2 = lambda2, niter = iter))
  }
}
