functions{
/*
		Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
		* diag(G) = [diagFirstLastElem,diagElem,...,diagElem,diagFirstLastElem]
		* first sup/sub-diagonal elements are all equal to offDiagElem
		
		Th routine returns a vector[2*n] L where 
		* L[1:n] is the diagonal of L
		* L[n+1:2*n-1] is the sub-diagonal of L
		* L[2*n] is the log-determinant of L
	*/
	vector CIP_TriDiagChol_const1n(int n, real diagFirstLastElem, real diagElem, real offDiagElem){
		vector[2*n] L;
		real LlogDet;
		// first iteration
		L[1] = sqrt(diagFirstLastElem);
		LlogDet = log(L[1]);
		// iteration 2:n-1
		for ( t in 2:(n-1)){
			L[n+t-1] = offDiagElem/L[t-1];
			L[t] = sqrt(diagElem - pow(L[n+t-1],2));
			LlogDet += log(L[t]);
		}
		// last iteration
		L[2*n-1] = offDiagElem/L[n-1];
		L[n] = sqrt(diagFirstLastElem - pow(L[2*n-1],2));
		LlogDet += log(L[n]);
		// done Cholesky
		
		L[2*n] = LlogDet;
		return(L);
	}
	
	/*
		Solves L^T x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_LT_solve(vector L, vector b){
		int n = rows(b);
		vector[n] x;
		// first solve
		x[n] = b[n]/L[n];
		// remaining solves
		for ( tt in 1:(n-1)){
			x[n-tt] = (b[n-tt] - x[n-tt+1]*L[2*n-tt])/L[n-tt];
		}
		return(x);
	}
	
	/*
		Solves L x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_L_solve(vector L, vector b){
		int n = rows(b);
		vector[n] x;
		// first solve
		x[1] = b[1]/L[1];
		// remaining solves
		for ( i in 2:n){
			x[i] = (b[i] - x[i-1]*L[n+i-1])/L[i];
		}
		return(x);
	}
	
	/*
		Solves L L^T x = G x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_LLT_solve(vector L, vector b){
		return(CIP_TriDiagChol_LT_solve(L,CIP_TriDiagChol_L_solve(L,b)));
	}
	
}

data{
	int T;
	matrix[5*T,5] yy;
	matrix[5*T,5] yinv;
	real ldets;
	int niter;
}

parameters{
	vector[5] mu;
	vector<lower=0>[5] sigmaSq;
	vector<lower=-1,upper=1>[5] delta;
	vector[10] hts;
	matrix[T,5] zs;
	real<lower=6.0> nu;
}

transformed parameters{

	vector[5] sigma;
	matrix[T,5] xs;
	matrix[T,5] ytilde;
	matrix[5,5] H;
	vector[T] lhs;
	vector[T] Gxhx;
	vector[T] h;
	vector[T] grad;
	vector[T-1] devs;
	vector[2*T] L;
	real logdet;
	real tmp;
	vector[5] gradMean;
	
	
	sigma = sqrt(sigmaSq);
	
	// build H-matrix
	for(i in 1:5){
	  H[i,i] = 1.0;
	}
	for(j in 1:4){
	  for(i in (j+1):5){
	    H[j,i] = 0.0;
	  }
	}
	H[2,1] = hts[1];
	H[3,1] = hts[2];
	H[4,1] = hts[3];
	H[5,1] = hts[4];
	
	H[3,2] = hts[5];
	H[4,2] = hts[6];
	H[5,2] = hts[7];
	
	H[4,3] = hts[8];
	H[5,3] = hts[9];
	
	H[5,4] = hts[10];
	
	// compute quasi-observations ytilde
	for(i in 1:5){
	  for(t in 1:T){
	    ytilde[t,i] = sub_col(H,i,i,6-i)' * block(yinv,(t-1)*5+i,i,6-i,6-i) * sub_col(H,i,i,6-i);
	  }
	}
	
	// zeroth iteration
	
	logdet = 0.0;
	
	for(i in 1:5){
	  gradMean[i] = 999.0;
	  tmp = mu[i]*square(1.0-delta[i])/sigmaSq[i];
	  for(t in 2:(T-1)) Gxhx[t] = tmp;
	  Gxhx[1] = mu[i]*(1.0-delta[i])/sigmaSq[i];
	  Gxhx[T] = Gxhx[1];
	  lhs = Gxhx - 0.5*nu*log((1.0/nu)*col(ytilde,i));
	  L = CIP_TriDiagChol_const1n(T,0.5*nu+1.0/sigmaSq[i],0.5*nu+(1.0+square(delta[i]))/sigmaSq[i],-delta[i]/sigmaSq[i]);
	  h = CIP_TriDiagChol_LLT_solve(L,lhs);
	  // further iterations on h here, if niter>0
	  for(iter in 1:niter){
	    devs = (1.0/sigmaSq[i])*(h[2:T] - mu[i] - delta[i]*(h[1:(T-1)]-mu[i]));
	    grad[2:T] = -devs[1:(T-1)] - 0.5*ytilde[2:T,i].*exp(h[2:T]) + 0.5*nu ;
	    grad[1] = -(h[1]-mu[i])*(1.0-square(delta[i]))/sigmaSq[i] - 0.5*ytilde[1,i]*exp(h[1]) + 0.5*nu;
	    grad[1:(T-1)] = grad[1:(T-1)] + delta[i]*devs;
	    gradMean[i] = mean(fabs(grad));
	    h = h + CIP_TriDiagChol_LLT_solve(L,grad);
	  }
	  
	  xs[1:T,i] = h + CIP_TriDiagChol_LT_solve(L,zs[1:T,i]);
	  logdet += L[2*T];
	}
}

model{
  
  target += inv_gamma_lpdf(sigmaSq | 2.0, 0.5);
  target += normal_lpdf(hts| 0.0, 10.0);
  target += normal_lpdf(mu | 0.0, 5.0);
  
	for(i in 1:5) {
	  target += normal_lpdf(xs[1,i] | mu[i], sigma[i]/sqrt(1.0-square(delta[i])));
	  target += normal_lpdf(xs[2:T,i] | mu[i] + delta[i]*(xs[1:(T-1),i]-mu[i]),sigma[i]);
	}
	
	// common factors in the inverse wishart lpdf
	target += -0.5*(nu+6.0)*ldets + 0.5*nu*sum(xs) - 2.5*T*0.6931471805599453*nu;
	for(i in 1:5){
	  target += -T*lgamma(0.5*(nu+1.0-i));
	}
	// trace factor
	for(t in 1:T){
	  target += -0.5*trace(H * diag_pre_multiply(exp(row(xs,t)),H') * block(yinv,(t-1)*5+1,1,5,5));
	  //target += inv_wishart_lpdf(block(yy,(t-1)*5+1,1,5,5) | nu , H * diag_pre_multiply(exp(row(xs,t)),H'));
	}
	
	target += -logdet;
	
}

generated quantities {

	vector[5] firstx;
	vector[5] firstz;
	vector[5] lastz;
	
	firstx = xs[1,1:5]';
	firstz = zs[1,1:5]';
	lastz = zs[T,1:5]';
}

