
# To get the input of 'prepars' in DREst() function

LSEst = function(Y,A,Z,W,max.step = 1000,thres = 1e-6,scale = c("add","mult"),
                 psi.start = NULL, alpha.start = NULL){
  pa = dim(Z)[2]
  pb = dim(W)[2]
  if(is.null(psi.start)) psi.start = rep(0, pa)
  if(is.null(alpha.start)) alpha.start = rep(0, pb)
  
  if(scale == "add"){
    
    rss.psi = function(psi){
      p0p1 = getProbRD(Z %*% psi, W %*% alpha) # referred to brm package
      p0    = p0p1[,1];  p1 = p0p1[,2]
      return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
    }
    
    rss.alpha = function(alpha){
      p0p1 = getProbRD(Z %*% psi, W %*% alpha) # referred to brm package
      p0    = p0p1[,1];  p1 = p0p1[,2]
      return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
    }
    
    rss = function(psi,alpha){
      p0p1 = getProbRD(Z %*% psi, W %*% alpha) # referred to brm package
      p0    = p0p1[,1];  p1 = p0p1[,2]
      return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
    }
  }
  
  if(scale == "mult"){
    
    rss.psi = function(psi){
      p0p1 = getProbRR(Z %*% psi, W %*% alpha) # referred to brm package
      p0    = p0p1[,1];  p1 = p0p1[,2]
      return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
    }
    
    rss.alpha = function(alpha){
      p0p1 = getProbRR(Z %*% psi, W %*% alpha) # referred to brm package
      p0    = p0p1[,1];  p1 = p0p1[,2]
      return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
    }
    
    rss = function(psi,alpha){
      p0p1 = getProbRR(Z %*% psi, W %*% alpha) # referred to brm package
      p0    = p0p1[,1];  p1 = p0p1[,2]
      return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
    }
  }
  
  ## Optimization
  Diff = function(x,y) {sum((x-y)^2)/sum(x^2+thres)}
  psi = psi.start; alpha = alpha.start
  diff = thres + 1; step = 0
  while(diff > thres & step < max.step){
    step = step + 1
    
    if(pa==1){
      opt1 = optimize(rss.psi,interval = c(-10,10),
                      tol = thres)
      diff1 = Diff(opt1$minimum,psi)
      psi = opt1$minimum
    }else{
      opt1 = optim(psi,rss.psi,
                   control=list(maxit=max(100,max.step/10)))
      diff1 = Diff(opt1$par,psi)
      psi = opt1$par
    }
    
    if(pb==1){
      opt2 = optimize(rss.alpha,interval = c(-3,3))
      diff  = max(diff1,Diff(opt2$minimum,alpha))
      alpha = opt2$minimum
    }else{
      opt2 = optim(alpha,rss.alpha,
                   control=list(maxit=max(100,max.step/10)))
      diff  = max(diff1,Diff(opt2$par,alpha))
      alpha = opt2$par
    }
    
  }
  
  pars2 = list(psi = psi,alpha=alpha, value = rss(psi,alpha),
               convergence = (step < max.step), step = step, Diff = diff)
  
  return(pars2)
}
