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

data {
int<lower=0> N; // # time points (equally spaced)
vector[N] y; // return at time t
}

parameters {
real tau;
real beta;
real<lower=0,upper=1> delta_p; // persistence of volatility
real<lower=0> v_p; // white noise shock scale
vector[N] u; 
}

transformed parameters {
vector[N] x; 
vector[N] GxGy;
vector[N] h;
vector[N] grad;
vector[N] diag;
vector[2*N] L;
real delta;
real vSq;
real deltaSq;

delta = -1.0 + 2.0*delta_p;
deltaSq = pow(delta,2);
vSq = 0.01*10.0/v_p;

L = CIP_TriDiagChol_const1n(N, (1.0/tau)+(1.0/vSq), (1.0/tau)+(1.0+deltaSq)/vSq, -delta/vSq);

// G_x*h_x + G_y*h_y
GxGy = (1.0/tau)*log(y/beta);

h = CIP_TriDiagChol_LLT_solve(L,GxGy);

x = h+CIP_TriDiagChol_LT_solve(L,u);

}

model {

v_p ~ chi_square(10);
delta_p ~ beta(20,1.5);
	
// latent process
// t=1
target += normal_lpdf(x[1] | 0 , sqrt(vSq/(1-deltaSq)));
//print("target_x1: ", target());

// remaining times
for( t in 2:N){
	target += normal_lpdf(x[t] | delta*x[t-1],sqrt(vSq));
	//print("target_xi: ", target());
}

// observations

target += gamma_lpdf(y | 1/tau,1.0 ./ (tau*exp(x)*beta));
//print("target_y: ", target());

// Jacobian correction

target += -L[2*N];

}

generated quantities{
real v = sqrt(vSq);
}
