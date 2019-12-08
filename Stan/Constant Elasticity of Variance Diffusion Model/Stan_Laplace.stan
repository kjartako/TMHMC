functions{
	/*
         Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where
         * diag(G) = diagElem (diagElem has length n)
         * first sup/sub-diagonal elements are in offDiagElem ( offDiagElem has length n-1)

         Th routine returns a vector[2*n] L where
         * L[1:n] is the diagonal of L
         * L[n+1:2*n-1] is the sub-diagonal of L
         * L[2*n] is the log-determinant of L
     */
     vector CIP_TriDiagChol(vector diagElem, vector offDiagElem){
         int n = rows(diagElem);
         vector[2*n] L;
         real LlogDet;
         L[1] = sqrt(diagElem[1]);
         LlogDet = log(L[1]);
         for ( t in 2:n){
             L[n+t-1] = offDiagElem[t-1]/L[t-1];
             L[t] = sqrt(diagElem[t] - pow(L[n+t-1],2));
             LlogDet += log(L[t]);
         }
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
     vector CIP_TriDiagChol_LLT_solve(vector L, vector b){ return(CIP_TriDiagChol_LT_solve(L,CIP_TriDiagChol_L_solve(L,b)));
     }
}

data{ 
 int N;
 vector[N] y;
 real Delta;
 int<lower=0> Nreg;
}

transformed data{
	vector[N] ySq;
	real sqrtDelta;
	ySq = square(y);
	sqrtDelta = sqrt(Delta);
}

parameters{
	real alpha;
	real beta;
	real<lower=0,upper=4> gamma;
	real logSxSq;
	real logSySq;
	vector[N] u; 
}

transformed parameters{
	real sxsq;
	real sysq;
	real sx;
	real sy;
	real gammasq;
	real Delbeta;
	real DelbetaSq;
	
	vector[N] diag;
	vector[N] grad;
	vector[N-1] odiag;
	vector[N] x;
	vector[N] h;
	vector[N] hSq;
	vector[N-1] tmp;
	vector[N-1] tmpSq;
	vector[N-1] tmp2;
	vector[2*N] Lnewton;
	
	// variance parameters
	sxsq = exp(logSxSq);
	sysq = exp(logSySq);
	sx = exp(0.5*logSxSq);
	sy = exp(0.5*logSySq);
	gammasq = pow(gamma,2);
	Delbeta = (Delta*beta-1.0);
	DelbetaSq = pow(Delta*beta-1.0,2);
	
	// eval point
	h=y;
	hSq=ySq;
	
	// Newton iterations

	for(w in 1:Nreg)
	{
		tmp[1:(N-1)] = h[2:N] + Delbeta*h[1:(N-1)] - Delta*alpha;
		tmpSq = square(tmp);
		
		for(i in 1:(N-1))
		{
			tmp2[i] = pow(h[i],2.0*gamma+2.0);
		}
		tmp2 *= sxsq*Delta;
		
		// build Hessian
		// diagonal
		diag[1:(N-1)] = ((DelbetaSq*hSq[1:(N-1)]-4.0*Delbeta*gamma*h[1:(N-1)] .* tmp[1:(N-1)] + (2.0*gammasq+gamma)*tmpSq[1:(N-1)]) ./ tmp2[1:(N-1)])-(gamma ./ hSq[1:(N-1)]);
		diag[N] = 0.0;
		diag[1] += 10000.0;
		diag[2:N] += hSq[1:(N-1)] ./ tmp2[1:(N-1)];
		diag[1:N] += 1.0/sysq;
	
		// first off-diagonal
		odiag = (Delbeta*hSq[1:(N-1)]-2.0*gamma*h[1:(N-1)] .* tmp[1:(N-1)]) ./ tmp2[1:(N-1)];
	
		// Cholesky
		Lnewton = CIP_TriDiagChol(diag,odiag);
		
		// build gradient

		grad[1] = -10000.0*h[1]+10000.0*y[1] + (tmpSq[1]*gamma*h[1] - tmp[1]*Delbeta*hSq[1])/tmp2[1]  - (gamma/h[1]) + ((y[1]-h[1])/sysq);
		
		grad[2:(N-1)] = - ((tmp[1:(N-2)] .* hSq[1:(N-2)]) ./ tmp2[1:(N-2)]) + (gamma*tmpSq[2:(N-1)] .* h[2:(N-1)] - Delbeta*tmp[2:(N-1)] .* hSq[2:(N-1)]) ./ tmp2[2:(N-1)]  - (gamma ./ h[2:(N-1)]) + ((y[2:(N-1)]-h[2:(N-1)])/sysq);
		
		grad[N] = - ((tmp[N-1]*hSq[N-1]) ./ tmp2[N-1]) + ((y[N]-h[N])/sysq);
		
		// update eval point

		h = h + CIP_TriDiagChol_LLT_solve(Lnewton,grad);	
		hSq = square(h);
	}
	
	x = h + CIP_TriDiagChol_LT_solve(Lnewton,u);	
}

model{
alpha ~ normal(0,sqrt(1000));
beta ~ normal(0,sqrt(1000));

// latent process
// t=1
target += normal_lpdf(x[1] | y[1] , 0.01);

// remaining times
for( t in 2:N){
	target += normal_lpdf(x[t] | x[t-1]+Delta*(alpha-beta*x[t-1]),sx*pow(x[t-1],gamma)*sqrtDelta);
}

// observations

target += normal_lpdf(y | x,sy);

// Jacobian correction

target += -Lnewton[2*N];

}

