#####################################################################################################
# This script was written by Ibtihal Ferwana, for any inquiries please contact iferwna2@illinois.edu
#####################################################################################################


#Loads required package, installing if needed
for(pack in c("nnls","quadprog", "parallel","Matrix", "pracma"))
  if(!require(pack, character.only = T))
  {
    install.packages(pack)
    require(pack, character.only = T)
  }
library(pracma)

# Optimal Recovery Algorithm
get_ubar <- function(A, u, xx, reg)
  {
  y = as.numeric(u[xx])
  A = as.matrix(A)
  N = nrow(A)
  T = ncol(A)
  I = as.matrix(diag(nrow=T))
  R = t(A)%*%A + I*reg
  Q = pinv(R)
  PHI = t(R[,xx])%*%Q%*%R[,xx]
  w = inv(PHI)%*%y
  ubar = R[,xx]%*%w
  return(list("ubar"=ubar, "PHI"=PHI, "R"=R, "Q"=Q))
}

# function for deriving error bars
error_bars <- function(PHI, xx, R, Q, A, ubar)
{
  PHI_inv = inv(PHI)
  TT = ncol(A)
  N = nrow(A)
  cbar = matrix(0,length(xx),TT)
  cseq = seq(1,TT)
  for(ii in cseq)
  {
    cbar[,ii] = PHI_inv%*%t(R[,xx])%*%Q%*%R[,ii]
  }
  eps = 0.001
  ybar = matrix(0,TT,TT)
  for(ii in cseq)
  {
    ybar[,ii] = R[,c(ii,xx)]%*%c(1,-cbar[,ii])
    if(sum(ybar[,ii])>eps)
    {
      ybar[,ii] = ybar[,ii]/sqrt(t(ybar[,ii])%*%Q%*%ybar[,ii])
    }
  }
  a = abs(rnorm(N,0,1))
  scale = sqrt(abs(t(a)%*%a - t(ubar)%*%Q%*%ubar))
  
  uworst1 = matrix(0,TT,TT)
  uworst2 = matrix(0,TT,TT)
  
  for(ii in cseq)
  {
    uworst1[,ii] = t(ubar)+scale*ybar[,ii]
    uworst2[,ii] = t(ubar)-scale*ybar[,ii]
  }
  max_error = c(scale)*abs(diag(ybar))
  return(list("max_error"=max_error))
}

# Synthetic Control algorithm
synth_control_est <-function(A, u, xx)
{
  full_T = ncol(A)
  T_0 = tail(xx,1)
  A = as.matrix(A)
  
  A_before = A[,xx]
  X = t(A_before)
  y_after = t(A)
  y = as.numeric(u[xx])
  
  
  Dmat = t(X)%*%X
  
  # if(is.positive.definite(Dmat, tol=1e-8)==FALSE){
  #   print("if-else FALSE")
  #   PDmat  = nearPD(Dmat)
  #   Dmat = as.matrix(PDmat$mat)
  # }
  if(isposdef(Dmat, tol=1e-8)==FALSE){
    print("if-else FALSE")
    PDmat  = nearPD(Dmat)
    Dmat = as.matrix(PDmat$mat)
  }
  dvec = t(X)%*%y
  
  Amat = t(rbind(rep(1,ncol(X)),diag(ncol(X))))
  bvec = c(1, rep(0,ncol(X)))
  
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  w = synth_model$solution
  # effects = y_after%*%c(1,-w)
  effects = y_after%*%w
  
  return(list("w" = w, "effects" = effects, "value_min" = (2*synth_model$value +y%*%y) ))
  
}