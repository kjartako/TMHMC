data {
int<lower=0> N; // # time points (equally spaced)
vector[N] y; // 
}

parameters {
real tau;
real beta;
real<lower=0,upper=1> delta_p; // persistence of volatility
real<lower=0> v_p; // white noise shock scale
vector[N] h_std; // std log volatility time t
}

transformed parameters {
real delta;
real v;
vector[N] h; // log volatility at time t

delta = -1.0 + 2.0*delta_p;
v = sqrt(0.01*10.0/v_p);
h = h_std * v;
h[1] = (h[1] / sqrt(1 - delta * delta));

for (t in 2:N)
h[t] = h[t]+ delta*h[t - 1];
}

model {

v_p ~ chi_square(10);
delta_p ~ beta(20,1.5);
h_std ~ normal(0, 1);

y ~ gamma(1/tau,1 ./ (tau*exp(h)*beta));

}
