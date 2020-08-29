# /*
#   Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
# * diag(G) = [diagFirstLastElem,diagElem,...,diagElem,diagFirstLastElem]
# * first sup/sub-diagonal elements are all equal to offDiagElem
# 
# Th routine returns a vector[2*n] L where 
# * L[1:n] is the diagonal of L
# * L[n+1:2*n-1] is the sub-diagonal of L
# * L[2*n] is the log-determinant of L
# */
CIP_TriDiagChol_const1n <- function(n,diagFirstLastElem,diagElem,offDiagElem){
    L <- rep(0.0,2*n)
    LlogDet <- 0.0;
    #// first iteration
    L[1] <- sqrt(diagFirstLastElem);
    LlogDet <- log(L[1]);
    #// iteration 2:n-1
    for ( t in 2:(n-1)){
      L[n+t-1] <- offDiagElem/L[t-1];
      L[t] <- sqrt(diagElem - L[n+t-1]^2);
      LlogDet <- LlogDet + log(L[t]);
    }
    #// last iteration
    L[2*n-1] <- offDiagElem/L[n-1];
    L[n] <- sqrt(diagFirstLastElem - L[2*n-1]^2);
    LlogDet <- LlogDet + log(L[n]);
    #// done Cholesky
    
    L[2*n] = LlogDet;
    return(L);
  }

CIP_TriDiagChol <- function(diagElem, offDiagElem){
  n <- length(diagElem);
  L <- rep(0.0,2*n);
  LlogDet <- 0.0;
  L[1] <- sqrt(diagElem[1]);
  LlogDet <- log(L[1]);
  for ( t in 2:n){
    L[n+t-1] <- offDiagElem[t-1]/L[t-1];
    L[t] <- sqrt(diagElem[t] - L[n+t-1]^2);
    LlogDet <- LlogDet + log(L[t]);
  }
  L[2*n] <- LlogDet;
  return(L);
}



#/*
#  Solves L^T x = b for x when L is the output of one of the tridiagonal
#Cholesky factorizations above
#*/
CIP_TriDiagChol_LT_solve <- function(L,b){
    n <- length(b)
     x <- rep(0.0,n)
    #// first solve
    x[n] <- b[n]/L[n];
    #// remaining solves
    for ( tt in 1:(n-1)){
      x[n-tt] <- (b[n-tt] - x[n-tt+1]*L[2*n-tt])/L[n-tt];
    }
    return(x);
  }

#/*
#  Solves L x = b for x when L is the output of one of the tridiagonal
#Cholesky factorizations above
#*/
CIP_TriDiagChol_L_solve <- function(L,b){
    n <- length(b);
    x <- rep(0.0,n);
    #// first solve
    x[1] <- b[1]/L[1];
    #// remaining solves
    for ( i in 2:n){
      x[i] <- (b[i] - x[i-1]*L[n+i-1])/L[i];
    }
    return(x);
  }

#/*
#  Solves L L^T x = G x = b for x when L is the output of one of the tridiagonal
#Cholesky factorizations above
#*/
CIP_TriDiagChol_LLT_solve <- function(L, b){ 
  return(CIP_TriDiagChol_LT_solve(L,CIP_TriDiagChol_L_solve(L,b)));
}