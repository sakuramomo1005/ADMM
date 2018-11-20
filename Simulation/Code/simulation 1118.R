
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 1115 codes')
source('admm.R')
source('ama.R')

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 1116')

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

# data generation
data_gen=function(seed){
  set.seed(seed)
  n = 80
  p = 40
  mu = seq(-60,60,10)
  theta = 2
  x = matrix(0, n, p)
  
  # 2 row groups
  # 8 column groups
  row_group = 4
  col_group = 4
  
  mu_rc = matrix(0, row_group,col_group)
  for(r in 1:row_group){
    for(c in 1:col_group){
      mu_rc[r,c] = sample(mu)[1]
    }
  }
  
  row_assign = c()
  col_assign = c()
  for(i in 1:n){
    row_assign = c(row_assign, sample(1:4)[1])
  }
  for(i in 1:p){
    col_assign = c(col_assign, sample(1:4)[1])
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


phi = 0.5;kk=5
nu1 = nu2 = 3

ama_rand = c()
corba_rand =c()


for(sim in 1:500){
  print('sim')
  print(sim)
  # generate simulation count table
  num.samples= 50; mu = 10000; 
  #num.samples: sample size; mu: total read count for each sample
  sim.otu.tab <- simPop(J=num.samples, n= mu, pi=fit$pi,
                        theta=fit$theta)$data  
  ##Compositional data:
 X = sim.otu.tab/mu
 # X <- X - mean(X)
  #X <- X/norm(X,'f')
  begin = Sys.time()
  # tune correct for ama
  seq1 = seq2 = c(1,10,100,1000)
  seqs_old = matrix(0,3,3)
  for(tunes in 1:3){
    seqs = c()
    for(gamma_1 in seq1){
      for(gamma_2 in seq2){
        res2 = bi_ama_order(X, nu1, nu2, gamma_1, gamma_2, kk,niter = 1000, tol = 1e-5,output = 0)
        ind2 = rand.index(X,res2$A)
        seqs = rbind(seqs,c(gamma_1,gamma_2,ind2))
      }
    }
    if(max(seqs[,3])>0.999){
      temp = seqs[which(seqs[,3]== max(seqs[,3]))[1],]
      g1 = temp[1];g2 = temp[2]
      ama_res = max(seqs[,3])
      gamma1 = g1; gamma2 = g2
      break
    }
    if(max(seqs_old[,3]) > max(seqs[,3])){
      temp = seqs[which(seqs_old[,3]== max(seqs_old[,3]))[1],]
      g1 = temp[1];g2 = temp[2]
      gamma1 = g1; gamma2 = g2
      ama_res = max(seqs_old[,3])
      break
    }else{
      temp = seqs[which(seqs[,3]== max(seqs[,3]))[1],]
      g1 = temp[1];g2 = temp[2]
      seq1 = seq(g1 - g1 * 0.2,g1 + g1 * 0.2,g1 * 0.1)
      seq2 = seq(g2 - g2 * 0.2,g2 + g2 * 0.2,g2 * 0.1)
      #print(seqs)
      seqs_old = seqs
    }
  }
  
  if(tunes == 3){
    temp = seqs[which(seqs[,3]== max(seqs[,3]))[1],]
    ama_res = max(seqs[,3])
    g1 = temp[1];g2 = temp[2]
    gamma1 = g1; gamma2 = g2
  }
  
  end = Sys.time()
  begin - end
  # save final results
  ama_rand = rbind(ama_rand,c(gamma1,gamma2,ama_res))
  
  
  if(sim %% 10 ==1){
    save(ama_rand,file = paste('sim',sim,'.RData',sep=''))
  }
  
  # cobra
  # tune correct for corba
  begin2 = Sys.time()
  phi <- 0.5; k <- 5
  wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
  w_row <- wts$w_row
  w_col <- wts$w_col
  E_row <- wts$E_row
  E_col <- wts$E_col
  
  ## Connected Components of Row and Column Graphs
  wts$nRowComp
  wts$nColComp
  
  seqs = c(1,10,100,1000)
  inds = c()
  inds_old = 0
  for(tunes in 1:3){
    for(gamma in seqs){
      sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gamma,fraction=0.01)
      inds = c(inds,rand.index(X,sol$U[[1]]))
    }
    if(max(inds) > 0.999){
      print('bbb')
      gamma = seqs[which(inds == max(inds))][1]
      corba_res =  max(inds)
      break
    }
    if(max(inds_old) > max(inds)){
      print('aaa')
      out_g = g
      corba_res = max(inds_old)
      break
    }else{
      g = seqs[which(inds == max(inds))][1]
      seqs = seq(g - g * 0.2,g + g * 0.2,g * 0.1)
      inds_old = inds
    }
  }
  
  if(tunes == 3){
    gamma = g
    corba_res = max(inds_old)
  }
  
  end2 = Sys.time()
  corba_rand = rbind(corba_rand, c(gamma,corba_res))
  
  if(sim %% 10 ==1){
    save(corba_rand,file = paste('simcorba',sim,'.RData',sep=''))
  }
}


cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)
names = corba_rand[,1]

pdf('1111.pdf')
heatmap(sol$U[[1]],col=hmcols,labRow=NA,labCol=NA,main = paste(names,'sol'))
heatmap(res2$A,col=hmcols,labRow=NA,labCol=NA,main = paste(names,'ama'))
dev.off()

