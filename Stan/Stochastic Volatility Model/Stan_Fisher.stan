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
}

data {
int<lower=0> N; // # time points (equally spaced)
vector[N] y; // return at time t
}

transformed data {
vector[N] ySq;

ySq = square(y);
}

parameters {
real gamma; // mean log volatility
real<lower=0,upper=1> delta_p; // persistence of volatility
real<lower=0> v_p; // white noise shock scale
vector[N] u; 
}

transformed parameters {
vector[N] x;
vector[2*N] L;
real delta;
real vSq;
real deltaSq;

delta = -1.0 + 2.0*delta_p;
deltaSq = pow(delta,2);
vSq = 0.01*10.0/v_p;

L = CIP_TriDiagChol_const1n(N, 0.5+(1.0/vSq), 0.5+(1.0+deltaSq)/vSq, -delta/vSq);

x = CIP_TriDiagChol_LT_solve(L,u);
}

model {

v_p ~ chi_square(10);
delta_p ~ beta(20,1.5);
	
// latent process
// t=1
target += normal_lpdf(x[1] | gamma/(1-delta) , sqrt(vSq/(1-deltaSq)));

// remaining times
for( t in 2:N){
	target += normal_lpdf(x[t] | gamma+delta*x[t-1],sqrt(vSq));
}

// observations

target += normal_lpdf(y | 0,exp(0.5*x));

// Jacobian correction

target += -L[2*N];

}

generated quantities{
real v = sqrt(vSq);
}
