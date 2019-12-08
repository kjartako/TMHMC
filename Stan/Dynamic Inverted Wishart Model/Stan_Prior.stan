data{
	int T;
	matrix[5*T,5] yy;
	matrix[5*T,5] yinv;
	real ldets;
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
	matrix[5,T] xs;
	matrix[5,5] H;
	
	for(i in 1:5){
	  sigma[i] = sqrt(sigmaSq[i]);
	  xs[i,1] = mu[i] + sigma[i]/sqrt(1.0-square(delta[i]))*zs[1,i];
	  for(t in 2:T){
	      xs[i,t] = mu[i] + delta[i]*(xs[i,t-1]-mu[i]) + sigma[i]*zs[t,i];
	  }
	}
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
	
}

model{
  
  target += inv_gamma_lpdf(sigmaSq | 2.0, 0.5);
  target += normal_lpdf(hts| 0.0, 10.0);
  target += normal_lpdf(mu | 0.0, 5.0);
  
	// common factors in the inverse wishart lpdf
	target += -0.5*(nu+6.0)*ldets + 0.5*nu*sum(xs) - 2.5*T*0.6931471805599453*nu;
	for(i in 1:5){
	  target += -T*lgamma(0.5*(nu+1.0-i));
	}
	// trace factor
	
	for(t in 1:T){
	  target += normal_lpdf( row(zs,t) | 0.0, 1.0);
	  target += -0.5*trace(H * diag_pre_multiply(exp(col(xs,t)),H') * block(yinv,(t-1)*5+1,1,5,5));
	 
	 // target += inv_wishart_lpdf(block(yy,(t-1)*5+1,1,5,5) | nu , H * diag_pre_multiply(exp(col(xs,t)),H'));
	}
}

generated quantities {

	vector[5] firstx;
	vector[5] firstz;
	vector[5] lastz;
	
	firstx = xs[1,1:5]';
	firstz = zs[1,1:5]';
	lastz = zs[T,1:5]';
}
