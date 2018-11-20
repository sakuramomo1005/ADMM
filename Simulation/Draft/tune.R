setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 1115 codes')
source('admm.R')
source('ama.R')


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

#nu1=1; nu2=1; gamma_1 = 100; gamma_2 = 100
#X = matrix(1:12,4,3)
#phi = 0.5
#kk=1
#res1 = bi_ADMM_order(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 1000, tol = 1e-5)
#res2 = bi_ama_order(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 1000, tol = 1e-5)
#
X = data_gen(123)

nu1 = nu2=3

for(gamma_1 in c(1,10,100,1000)){
  for(gamma_2 in c(1,10,100,1000)){
    phi = 0.5;kk=5;
    res1 = bi_ADMM_order(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 10000, tol = 1e-5)
    res2 = bi_ama_order(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 1000, tol = 1e-5)
    ind1 = rand.index(X,res1$A)
    ind2 = rand.index(X,res2$A)
    result = list(res1 = res1, res2 = res2, ind1 = ind1)
    names = paste('g1_',gamma_1,'g2_',gamma_2,'.RData',sep='')
    save(result,file = names)
  }
}