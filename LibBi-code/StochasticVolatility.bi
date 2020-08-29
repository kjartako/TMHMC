/**
 * Stochastic volatility model 
 */
model StochasticVolatility {
  param phi_v_prime, sigmaSq_v, gam_v;
  noise eta;
  state v;
  obs y;

  sub parameter {
    phi_v_prime ~ beta(20.0, 1.5);
    sigmaSq_v ~ inverse_gamma(5.0, 0.05);
    gam_v ~ uniform(-1.0, 1.0);
  }

  sub proposal_parameter {
    phi_v_prime ~ truncated_normal(phi_v_prime, 0.0048, 0.00001, 0.99999);
    sigmaSq_v ~ truncated_normal(sigmaSq_v, 0.0084, 0.00001, 100.0);
    gam_v ~ truncated_normal(gam_v, 0.011, -10.0, 10.0);
  }

  sub initial {
    /* stationary distribution for v */
    v ~ normal(gam_v/(1.0-(-1.0+2.0*phi_v_prime)), sqrt(sigmaSq_v/(1.0 - (-1.0+2.0*phi_v_prime)**2)));
  }

  sub transition {
    eta ~ normal();
    v <- gam_v + (-1.0+2.0*phi_v_prime)*v + sqrt(sigmaSq_v)*eta;
  }

  sub observation {
    y ~ normal(0.0, exp(0.5*v));
  }
}
