library(deSolve)
library(reshape2)


checkSize = function(param, L, W){
  param_name = deparse(substitute(param))
  param_size = dim(param)
  if(is.null(param_size)){
    print(paste0("Warning: ", param_name," is a scalar value, expanding to size ", L, "x", W, "."))
    param = param*matrix(1,L,W)
    param_size = dim(param)
  }
  else if(param_size[1]== W & param_size[2]==L & W!=L){
    print(paste0("Warning: ", param_name," was given in reverse dimension order, so transposing it before use..."))
    param = t(param)
    param_size = dim(param)
  }
  else if(param_size[1]!=L | param_size[2]!=W){
    print(paste0("Error: Parameter ",param_name," is of size ",param_size[1],"x",param_size[2]," and not ",L, "x",W))
    stop("See above message")
  }
return(param)
}

checkGreaterOrEqual = function(param, bound.l){
  if (sum(param<bound.l)>0) {
    param_name = deparse(substitute(param))
    print(paste0("Error: At least one of the values of ",param_name," is less than ",bound.l))
    stop("See message above")
  }
}

sirODE = function(times,init,params){
  with(as.list(c(params,init)), {
    dXh = nu1-r*(Tr12*Ym+Tr11*Yh)*Xh-mu1*Xh
    dXm = nu2-r*(Tr22*Ym+Tr21*Yh)*Xm-mu2*Xm
    dYh = r*(Tr12*Ym+Tr11*Yh)*Xh-mu1*Yh-gamma1*Yh
    dYm = r*(Tr22*Ym+Tr21*Yh)*Xm-mu2*Ym-gamma2*Ym
    return(list(c(dXh,dXm,dYh,dYm)))
  })
}
