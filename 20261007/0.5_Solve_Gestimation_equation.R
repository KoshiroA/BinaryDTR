
### Perform doubly robust g-estimation

# The INPUT of 'prepars' is obtained by performing LSEst() function in advance

# Packages required
# ---------------------------------------------------- #
if (requireNamespace("brm", quietly = TRUE)) {
  suppressWarnings(library(brm))
} else {
  stop("Required package 'brm' is not installed.")
}

if (requireNamespace("nleqslv", quietly = TRUE)) {
  suppressWarnings(library(nleqslv))
} else {
  install.packages("nleqslv")
}
# ---------------------------------------------------- #

DREst = function(Y,A,Z,W,prepars,PS = NULL,scale = c("add","mult"),
                 max.step = 1000,thres = 1e-6,startpars = NULL){
  
  if(is.null(PS)){
    PS.model = glm(A ~ W,family = binomial(link = "logit")) 
    PS = predict(PS.model,type="response")
  }
  
  if(is.null(startpars) == TRUE){
    psi = rep(0,dim(Z)[2])
  }else{
    psi = prepars$psi
  }
  
  ### Define Estimation Function ###
  if(scale == "add"){
    U.hat = brm::getProbRD(Z %*% prepars$psi,W %*% prepars$alpha)[,1]
    gest = function(psi.dr){
      tmp = t(Z) %*% as.vector((A-PS) * (Y - A* tanh(Z %*% psi.dr) - U.hat))
      return(tmp)
    }
  }
  if(scale == "mult"){
    U.hat = brm::getProbRR(Z %*% prepars$psi,W %*% prepars$alpha)[,1]
    gest = function(psi.dr){
      tmp = t(Z) %*% as.vector((A-PS) * (Y*exp(-A*Z %*% psi.dr) - U.hat))
      return(tmp)
    }
  }
  
  Diff = function(a,b) {sum((a-b)^2)/sum(a^2+thres)}
  diff = thres + 1; step = 0
  while(diff > thres*1000 & step < max.step){
    step = step + 1
    opt  = nleqslv(psi, gest, method = "Newton",
                    control = list(allowSingular=T))
    diff = Diff(opt$x,psi)
    psi = opt$x
  }
  
  opt$alpha = prepars$alpha
  opt$step = prepars$step
  
  if(step < max.step){
    return(opt)
  }else{
    stop("The estimation equation is not properly be solved")
  }
}
