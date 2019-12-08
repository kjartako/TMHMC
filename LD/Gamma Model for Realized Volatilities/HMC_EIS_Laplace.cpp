// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "../../adept.h"
#include "../../adept_source.h"
#include <math.h>
using adept::adouble;
typedef adept::adouble adtype;

const double c_s0 = 0.01;
const double c_df= 10;
const double c_alpha = 20;
const double c_beta = 1.5;

template <class T>
using Tmat = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

template <class T>
using Tvec = Eigen::Matrix<T,Eigen::Dynamic,1>;

/////////////////////////////

template <class T>
T gammalog(T x)
{ 
  Tvec<T> cof(14);
  cof <<  57.1562356658629235,-59.5979603554754912,14.1360979747417471,-0.491913816097620199,0.339946499848118887e-4,
                               0.465236289270485756e-4,-0.983744753048795646e-4,0.158088703224912494e-3,
                             -0.210264441724104883e-3,0.217439618115212643e-3,-0.164318106536763890e-3,
                              0.844182239838527433e-4,-0.261908384015814087e-4,0.368991826595316234e-5;
  
   if (x <= 0) throw("bad arg in gammln");
   T y=x;
   T tmp = x+5.24218750000000000; // 671/128.
   tmp = (x+0.5)*log(tmp)-tmp;
   T ser = 0.999999999999997092;
   for (int j=0;j<14;j++)
   {
     y = y + 1;
     ser = ser + cof(j)/y; 
   }
   
   return (tmp+log(2.5066282746310005*ser/x));
}


template <class T>
T Tdgammalog(double y, T tau, T beta, T l)
{
  T shape = 1.0/tau;
  T scale = tau*exp(l)*beta;
  T res = -gammalog(shape)-shape*log(scale)+(shape-1)*log(y)-(y/scale);

  return(res);
}

/////////////////////////////

template <class T>
T Tdchilog(T phi, T sigSq, T xPrev, T a1, T a2)
{
  T res = (-0.5*(2.0*a2*pow(phi,2)*pow(xPrev,2)+pow(a1,2)*sigSq+2.0*a1*phi*xPrev)/(2.0*a2*sigSq-1.0))-0.5*log(1.0-2.0*a2*sigSq);
  return(res);
}



/////////////////////////////
template <class T>
T Tdchilog1(T phi, T sigSq, T a1, T a2)
{
  T res = (-0.5*sigSq*pow(a1,2)/(2.0*a2*sigSq+pow(phi,2)-1.0))+0.5*log(1-pow(phi,2))-0.5*log(1.0-2.0*a2*sigSq-pow(phi,2));
  return(res);
}

/////////////////////////////


template <class T>
Tvec<T> cholsolve(Tmat<T> A, Tvec<T> b, int n)
{
  // Cholesky decomposition
  T sum=0;
  for (int i=0;i<n;i++)
  {
    for (int j=i;j<n;j++)
    {
      sum=A(i,j);
      for (int k=i-1;k>=0;k--)
      {
        sum -= A.coeff(i,k)*A.coeff(j,k);
      }
      if (i == j) {
        if (sum <= 0.0) // A, with rounding errors, is not positive-definite.
          throw("Cholesky failed");
        A.coeffRef(i,i)=sqrt(sum);
      } else A.coeffRef(j,i)=sum/A.coeff(i,i);
    }
  }

  Tvec<T> x(n);
  T sum2=0;
  // Solve A*x=b, storing y in x
  for (int i=0;i < n;i++)
  {
    sum2=b(i);
    for (int k=i-1;k>=0;k--)
    {
      sum2 -= A.coeff(i,k)*x.coeff(k);
    }
    x.coeffRef(i)=sum2/A.coeff(i,i);
  }

  // Solve L^T*x=y
  for (int i=n-1;i >= 0;i--)
  {
    sum2=x(i);
    for (int k=i+1;k<n;k++)
    {
      sum2 -= A.coeff(k,i)*x.coeff(k);
    }
    x.coeffRef(i)=sum2/A.coeff(i,i);
  }

  return(x);
}

/////////////////////////////

template <class T>
void LinReg(const Tvec<double> &y, Tmat<T> &A, const Tmat<T> &lambdaM, T tau,T beta,T phi,T sigSq,int n_z, int n_y)
{
  Tvec<T> yreg(n_y);
  Tmat<T> Xreg(n_z,3);

  Tvec<T> l=lambdaM.col(n_y-1);

  for (int j=0; j<n_z; j++)
  {
    T lj = l.coeff(j);
    yreg.coeffRef(j)=Tdgammalog(y.coeff(n_y-1), tau, beta, lj);

    Xreg.coeffRef(j,0)=1;
    Xreg.coeffRef(j,1)=lj;
    Xreg.coeffRef(j,2)=pow(lj,2);
  }

  Tmat<T> XregNorm =Xreg.transpose()*Xreg;
  Tvec<T> yregNorm =Xreg.transpose()*yreg;

  Tvec<T> reg=cholsolve(XregNorm,yregNorm,3);
  A.coeffRef((n_y-1),0)=reg.coeff(1);
  A.coeffRef((n_y-1),1)=reg.coeff(2);


  for (int i=(n_y-2); i>-1; i--)
  {
    Tvec<T> yreg(n_y);
    Tmat<T> Xreg(n_z,3);

    Tvec<T> l=lambdaM.col(i);

    for (int k=0; k<n_z; k++)
    {
      T lk = l.coeff(k);
      yreg(k)=Tdgammalog(y(i), tau, beta,lk)+Tdchilog(phi,sigSq,lk, A(i+1,0),A(i+1,1));

      Xreg(k,0)=1;
      Xreg(k,1)=lk;
      Xreg(k,2)=pow(lk,2);
    }

    Tmat<T> XregNorm =Xreg.transpose()*Xreg;
    Tvec<T> yregNorm =Xreg.transpose()*yreg;

    Tvec<T> reg=cholsolve(XregNorm,yregNorm,3);
    A.coeffRef(i,0)=reg.coeff(1);
    A.coeffRef(i,1)=reg.coeff(2);
  }

}

/////////////////////////////


template <class T>
T loglike_c (Tvec<double> y, Tvec<T> Par,int n_y, int n_z, Tvec<double> z, int n_reg,Tvec<T> u, int n_u){

  T tau = exp(Par(0));
  T beta = exp(Par(1));
  T phi = tanh(Par(2));
  T sigSq = exp(Par(3));

  // Start values for A
  Tmat<T> A(n_y,2);

  for (int i=0; i<n_y; i++)
  {
    A(i,0)=-log(beta/y(i))/tau;
    A(i,1)=-1.0/(2*tau);
  }

  Tmat<T> lambdaM(n_z,n_y);

  for(int j = 0; j < n_reg; j++)
  {
    T m1_var = -sigSq/(2.0*A(0,1)*sigSq+pow(phi,2)-1.0);
    T m1_mean = -(A(0,0)*sigSq)/(2.0*A(0,1)*sigSq+pow(phi,2)-1.0);

    // Sample from m
    for (int i=0;i<n_z;i++)
    {
      lambdaM.coeffRef(i,0) = m1_mean+z.coeff(i*n_y)*sqrt(m1_var);

      for (int k=1;k<n_y;k++)
      {
        T xPrev = lambdaM(i,k-1);
        T m_var = -sigSq/(2.0*A(k,1)*sigSq-1.0);
        T m_mean = -(A(k,0)*sigSq+phi*xPrev)/(2.0*A(k,1)*sigSq-1.0);

        lambdaM.coeffRef(i,k) = m_mean+z.coeff(i*n_y+k)*sqrt(m_var);
      }
    }

      // Regression to update A
      LinReg(y, A, lambdaM, tau,beta,phi,sigSq, n_z, n_y);
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  Tmat<T> lambdaM_u(n_u,n_y);

  T m1_var = -sigSq/(2.0*A(0,1)*sigSq+pow(phi,2)-1.0);
  T m1_mean = -(A(0,0)*sigSq)/(2.0*A(0,1)*sigSq+pow(phi,2)-1.0);

  for (int i=0;i<n_u;i++)
  {
    lambdaM_u.coeffRef(i,0) = m1_mean+u.coeff(i*n_y)*sqrt(m1_var);

    for (int k=1;k<n_y;k++)
    {
      T xPrev = lambdaM(i,k-1);
      T m_var = -sigSq/(2.0*A(k,1)*sigSq-1.0);
      T m_mean = -(A(k,0)*sigSq+phi*xPrev)/(2.0*A(k,1)*sigSq-1.0);

      lambdaM_u.coeffRef(i,k) = m_mean+u.coeff(i*n_y+k)*sqrt(m_var);
    }
  }

  Tvec<T> logsums(n_u);

  for(int i=0;i<n_u;i++)
  {
    Tvec<T> lw(n_y);

    T l = lambdaM_u.coeff(i,n_y-1);
    lw(n_y-1) = Tdgammalog(y.coeff(n_y-1),tau,beta,l) - (A(n_y-1,0)*l+A(n_y-1,1)*pow(l,2));

    for (int j=0;j<(n_y-1);j++)
    {
      T l = lambdaM_u.coeff(i,j);
      lw.coeffRef(j) = Tdgammalog(y.coeff(j), tau, beta,l) + Tdchilog(phi,sigSq, l, A(j+1,0),A(j+1,1)) - (A(j,0)*l+A(j,1)*pow(l,2));
    }
    logsums(i)= lw.sum()+Tdchilog1(phi,sigSq, A(0,0),A(0,1));
  }

  T w = logsums.maxCoeff();

  Tvec<T> lsexp(n_u);

  for (int k=0; k<n_u;k++)
  {
    lsexp.coeffRef(k) = exp(logsums.coeff(k) - w);
  }

  T loglike = w + log(lsexp.mean());
  
  T sig_prior = ((-0.5*c_df*c_s0)/sigSq)-0.5*c_df*log(sigSq);
  T phi_prior = log(0.5*phi+0.5)*(c_alpha-1)+log(0.5-0.5*phi)*(c_beta-1)+log(0.5-0.5*pow(phi,2));

  return(loglike+sig_prior+phi_prior);
}

/////////////////////////////
/////////////////////////////

Rcpp::List loglike_grad_c(Tvec<double> y,  Tvec<double> Par, int n_y, int n_z, Tvec<double> z, int n_reg, Tvec<double> u, int n_u)
{
  adept::Stack stack;

  Tvec<adtype> ad_Par(4);

  for(int k=0;k<4;k++)
  {
    ad_Par.coeffRef(k)=Par.coeff(k);
  }

  Tvec<adtype> ad_u(n_u*n_y);

  for(int i=0;i<(n_u*n_y);i++)
  {
    ad_u.coeffRef(i)=u.coeff(i);
  }

  stack.new_recording();

  adtype res0 = loglike_c<adtype>(y,ad_Par,n_y,n_z,z,n_reg,ad_u,n_u);
  adtype res = res0/1.0;

  res.set_gradient(1.0);
  stack.compute_adjoint();

  Tvec<double> Pargrads(4);

  for (int q=0;q<4;q++)
  {
    Pargrads.coeffRef(q)= ad_Par.coeff(q).get_gradient();
  }

  Tvec<double> ugrads(n_u*n_y);

  for (int j=0;j<(n_u*n_y);j++)
  {
    ugrads(j)= ad_u.coeff(j).get_gradient();
  }

  return Rcpp::List::create(
    Rcpp::Named("Par")  = Pargrads,
    Rcpp::Named("u")  = ugrads);
}

/////////////////////////////

double Hamiltonian(Tvec<double> Par, Tvec<double> y,int n_y,Tvec<double> u,Tvec<double> z,int n_u, int n_z, int n_reg,Tvec<double> p_Par,Tvec<double> p_u, Tmat<double> M_inv)
{
  double U = loglike_c<double>(y,Par,n_y,n_z,z,n_reg,u,n_u);

  double hamil = -U;

  hamil = hamil + 0.5*p_Par.transpose()*M_inv*p_Par;

  for (int i=0;i<(n_u*n_y);i++)
  {
    hamil = hamil + 0.5*(pow(u.coeff(i),2) + pow(p_u.coeff(i),2));
  }

  return(hamil);
}

/////////////////////////////

Rcpp::List HMC_integrator(Tvec<double> y, int n_y, Tvec<double> z,Tvec<double> Par0,Tvec<double> u0, int n_u, int n_z, Tvec<double> p_Par0,Tvec<double> p_u0,double epsilon,int L, double unif, int n_reg, Tmat<double> M_inv)
{
  Tvec<double> Par = Par0;
  Tvec<double> u = u0;
  Tvec<double> p_Par = p_Par0;
  Tvec<double> p_u = p_u0;
  Rcpp::List grad;

  for(int i=0;i<L;i++)
  {
    Tvec<double> p_u_tmp = cos(epsilon/2.0)*p_u - sin(epsilon/2.0)*u;

    //half step for positions
    Par = Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u;
    //
    //full step for momentum
    grad = loglike_grad_c(y,Par,n_y,n_z,z,n_reg,u,n_u);

    Tvec<double> grad_Par = grad[0];
    Tvec<double> grad_u = grad[1];

    p_Par = p_Par + epsilon*(grad_Par);
    Tvec<double> p_u_tmp2 = p_u_tmp + (epsilon*grad_u);
    p_u = cos(epsilon/2.0)*p_u_tmp2 - sin(epsilon/2.0)*u;

    //half step for positions
    Par= Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u_tmp2;
  }

  //Evaluate potential and kinetic energies at start and end of trajectory
  double H_current = Hamiltonian(Par0,y,n_y,u0,z,n_u,n_z,n_reg,p_Par0,p_u0,M_inv);
  double H_prop = Hamiltonian(Par,y,n_y,u,z,n_u,n_z,n_reg,p_Par,p_u,M_inv);

  //Accept or reject the state at end of trajectory
  double accp = exp(H_current-H_prop);
  if (unif < accp)
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par,
      Rcpp::Named("u")  = u,
      Rcpp::Named("acc")  = 1,
      Rcpp::Named("accprob")=accp);
  }
  else
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par0,
      Rcpp::Named("u")  = u0,
      Rcpp::Named("acc")  = 0,
      Rcpp::Named("accprob")=accp);
  }

}

// [[Rcpp::export]]
Rcpp::List HMC_EIS(int N_mcmc,Eigen::VectorXd Par0, Eigen::VectorXd y, int n_y, Eigen::VectorXd u0,Eigen::MatrixXd z,int n_u, int n_z, Eigen::MatrixXd p_Par,Eigen::MatrixXd p_u,Eigen::VectorXd epsilons,Eigen::VectorXd Ls,Eigen::VectorXd unif,int n_reg, Eigen::MatrixXd M_inv)
{
  Tmat<double> Par(N_mcmc+1,4);
  Tmat<double> u((N_mcmc+1),(n_u*n_y));
  Tvec<double> acc(N_mcmc);
  Tvec<double> accp(N_mcmc);

  Par.row(0)=Par0;
  u.row(0)=u0;

  for (int j=0;j<N_mcmc;j++)
  {
    Rcpp::List hmc_val;
    try {hmc_val=HMC_integrator(y,n_y,z.row(j),Par.row(j),u.row(j),n_u,n_z,p_Par.row(j),p_u.row(j),epsilons(j),Ls(j),unif(j),n_reg,M_inv);}
    catch(...){hmc_val=Rcpp::List::create(Rcpp::Named("Par")  = Par.row(j),Rcpp::Named("u")  = u.row(j),Rcpp::Named("acc")= 0,Rcpp::Named("accprob")= 0);}


    Tvec<double> Par_tmp = hmc_val[0];
    Par.row(j+1)=Par_tmp;
    Tvec<double> u_tmp = hmc_val[1];
    u.row(j+1)=u_tmp;
    double acc_tmp = hmc_val[2];
    acc(j)= acc_tmp;
    double accp_tmp = hmc_val[3];
    accp(j)= accp_tmp;
  }

  return Rcpp::List::create(
    Rcpp::Named("Par")  = Par,
    Rcpp::Named("u")  = u,
    Rcpp::Named("acc")= acc,
    Rcpp::Named("accprob")= accp);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where
* diag(G) = [diagFirstLastElem,diagElem,...,diagElem,diagFirstLastElem]
* first sup/sub-diagonal elements are all equal to offDiagElem

Th routine returns a vector[2*n] L where
* L[1:n] is the diagonal of L
* L[n+1:2*n-1] is the sub-diagonal of L
* L[2*n] is the log-determinant of L
*/
template <class T>
Tvec<T> CIP_TriDiagChol_const1n(int N, T diagFirstLastElem, T diagElem, T offDiagElem){
  Tvec<T> L(2*N);
  T LlogDet;
  L(0) = sqrt(diagFirstLastElem);
  LlogDet = log(L(0));
  // iteration 2:n-1
  for(int t=1;t<(N-1);t++)
  {
    L(N+t-1)= offDiagElem/L(t-1);
    L(t) = sqrt(diagElem - pow(L(N+t-1),2));
    LlogDet += log(L(t));
  }
  // last iteration

  L(2*N-2) = offDiagElem/L(N-2);
  L(N-1) = sqrt(diagFirstLastElem - pow(L(2*N-2),2));
  LlogDet += log(L(N-1));
  // done Cholesky

  L(2*N-1) = LlogDet;
  return(L);
}


/*
 Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where
* diag(G) = diagElem (diagElem has length n)
* first sup/sub-diagonal elements are all equal to offDiagElem

Th routine returns a vector[2*n] L where
* L[1:n] is the diagonal of L
* L[n+1:2*n-1] is the sub-diagonal of L
* L[2*n] is the log-determinant of L
*/
template <class T>
Tvec<T> CIP_TriDiagChol_diag_constod(Tvec<T> diagElem, T offDiagElem, int N)
{
  Tvec<T> L(2*N);
  T LlogDet;
  L(0) = sqrt(diagElem(0));
  LlogDet = log(L(0));

  for(int t=1;t<N;t++)
  {
    L(N+t-1)= offDiagElem/L(t-1);
    L(t) = sqrt(diagElem(t) - pow(L(N+t-1),2));
    LlogDet += log(L(t));
  }

  L(2*N-1) = LlogDet;
  return(L);
}

//Solves L^T x = b for x when L is the output of one of the tridiagonal
//Cholesky factorizations above
template <class T>
Tvec<T> CIP_TriDiagChol_LT_solve(Tvec<T> L, Tvec<T> b, int N)
{
  Tvec<T> x(N);
  // first solve
  x(N-1) = b(N-1)/L(N-1);
  // remaining solves
  for(int t=2;t<(N+1);t++)
  {
    x(N-t) = (b(N-t) - x(N-t+1)*L(2*N-t))/L(N-t);
  }
  return(x);
}


//Solves L x = b for x when L is the output of one of the tridiagonal
//Cholesky factorizations above
template <class T>
Tvec<T> CIP_TriDiagChol_L_solve(Tvec<T> L, Tvec<T> b, int N)
{
  Tvec<T> x(N);
  // first solve
  x(0) = b(0)/L(0);
  // remaining solves
  for(int i=1;i<N;i++)
  {
    x(i) = (b(i)-x(i-1)*L(N+i-1))/L(i);
  }
  return(x);
}

//Solves L L^T x = G x = b for x when L is the output of one of the tridiagonal
//Cholesky factorizations above
template <class T>
Tvec<T> CIP_TriDiagChol_LLT_solve(Tvec<T> L, Tvec<T> b, int N){ return(CIP_TriDiagChol_LT_solve(L,CIP_TriDiagChol_L_solve(L,b,N),N));
}

template <class T>
T loglike_Laplace(Tvec<double> y, Tvec<T> Par,int n_y, Tvec<T> u,int Newton_it)
{
  T tau = exp(Par(0));
  T beta = exp(Par(1));
  T delta = tanh(Par(2));
  T vSq = exp(Par(3));

  T tauinv = 1.0/tau;
  T vtb = vSq*tau*beta;
  T deltaSq = pow(delta,2);
  Tvec<T> x(n_y);
  Tvec<T> grad(n_y);
  Tvec<T> diag(n_y);
  Tvec<T> L(2*n_y);
  Tvec<T> GxGy(n_y);
  Tvec<T> h(n_y);

  /////////////////////////////////////////////////////////
  // Build Hessian
  T firstlast = tauinv+(1.0/vSq);
  T diagval = tauinv+(1.0+deltaSq)/vSq;
  T offdiag = -delta/vSq;

  L = CIP_TriDiagChol_const1n(n_y,firstlast,diagval,offdiag);

  // Build mean

  for(int i=0;i<n_y;i++)
  {
    GxGy(i) = tauinv*log(y(i)/beta);
  }

  h = CIP_TriDiagChol_LLT_solve(L,GxGy,n_y);

  // Newton iterations

	for(int w=0;w<Newton_it;w++)
	{
		diag(0) = (vSq*y(0)*exp(-h(0))+tau*beta)/vtb;

		for(int t=1;t<(n_y-1);t++)
		{
			diag(t) = (vSq*y(t)*exp(-h(t))+tau*beta*(1.0+deltaSq))/vtb;
		}

		diag(n_y-1) = (vSq*y(n_y-1)*exp(-h(n_y-1))+tau*beta)/vtb;

		L = CIP_TriDiagChol_diag_constod(diag,offdiag,n_y);

		grad(0) = (beta*(delta*tau*h(1)-tau*h(0)-vSq)+vSq*y(0)*exp(-h(0)))/vtb;

		for(int t=1;t<(n_y-1);t++)
		{
			grad(t) = (beta*delta*tau*(h(t-1)-delta*h(t)+h(t+1))-beta*(tau*h(t)+vSq)+vSq*y(t)*exp(-h(t)))/vtb;
		}

		grad(n_y-1) = (beta*(delta*tau*h(n_y-2)-tau*h(n_y-1)-vSq)+vSq*y(n_y-1)*exp(-h(n_y-1)))/vtb;

		h = h + CIP_TriDiagChol_LLT_solve(L,grad,n_y);
	}
	
  // Transformation
  x = h+CIP_TriDiagChol_LT_solve(L,u,n_y);
	
  //////////////////////////////////////////////////////////

  T var = vSq/(1-deltaSq);
  T loglike= -0.5*pow(x(0),2)/var - 0.5*log(2.0*M_PI*var) + Tdgammalog(y(0), tau, beta, x(0));
  loglike += 0.5*pow(u(0),2);
  
  for (int j=1;j<n_y;j++)
  {
    loglike += -0.5*pow(x(j)-delta*x(j-1),2)/vSq - 0.5*log(2.0*M_PI*vSq) + Tdgammalog(y(j), tau, beta, x(j));
    loglike += 0.5*pow(u(j),2);
  }
  loglike += -L(2*n_y-1);
  
  T v_prior = ((-0.5*c_df*c_s0)/vSq)-0.5*c_df*log(vSq);
  T delta_prior = log(0.5*delta+0.5)*(c_alpha-1)+log(0.5-0.5*delta)*(c_beta-1)+log(0.5-0.5*deltaSq);

  return(loglike+v_prior+delta_prior);
}

/////////////////////////////
/////////////////////////////

Rcpp::List loglike_grad_Laplace(Tvec<double> y,  Tvec<double> Par, int n_y, Tvec<double> u, int Newton_it)
{
  adept::Stack stack;

  Tvec<adtype> ad_Par(4);

  for(int k=0;k<4;k++)
  {
    ad_Par.coeffRef(k)=Par.coeff(k);
  }

  Tvec<adtype> ad_u(n_y);

  for(int i=0;i<n_y;i++)
  {
    ad_u.coeffRef(i)=u.coeff(i);
  }

  stack.new_recording();

  adtype res0 = loglike_Laplace<adtype>(y,ad_Par,n_y,ad_u,Newton_it);
  adtype res = res0/1.0;

  res.set_gradient(1.0);
  stack.compute_adjoint();

  Tvec<double> Pargrads(4);

  for (int q=0;q<4;q++)
  {
    Pargrads.coeffRef(q)= ad_Par.coeff(q).get_gradient();
  }

  Tvec<double> ugrads(n_y);

  for (int j=0;j<n_y;j++)
  {
    ugrads(j)= ad_u.coeff(j).get_gradient();
  }

  return Rcpp::List::create(
    Rcpp::Named("Par")  = Pargrads,
    Rcpp::Named("u")  = ugrads);
}

/////////////////////////////

double Hamiltonian_Laplace(Tvec<double> Par, Tvec<double> y,int n_y,Tvec<double> u,Tvec<double> p_Par,Tvec<double> p_u, Tmat<double> M_inv, int Newton_it)
{
  double U = loglike_Laplace<double>(y,Par,n_y,u,Newton_it);

  double hamil = -U;

  hamil = hamil + 0.5*p_Par.transpose()*M_inv*p_Par;

  for (int i=0;i<n_y;i++)
  {
    hamil = hamil + 0.5*(pow(u.coeff(i),2) + pow(p_u.coeff(i),2));
  }

  return(hamil);
}

/////////////////////////////

Rcpp::List HMC_integrator_Laplace(Tvec<double> y, int n_y, Tvec<double> Par0,Tvec<double> u0, Tvec<double> p_Par0,Tvec<double> p_u0,double epsilon,int L, double unif, Tmat<double> M_inv, int Newton_it)
{
  Tvec<double> Par = Par0;
  Tvec<double> u = u0;
  Tvec<double> p_Par = p_Par0;
  Tvec<double> p_u = p_u0;
  Rcpp::List grad;

  for(int i=0;i<L;i++)
  {
    Tvec<double> p_u_tmp = cos(epsilon/2.0)*p_u - sin(epsilon/2.0)*u;

    //half step for positions
    Par = Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u;
    //
    //full step for momentum
    grad = loglike_grad_Laplace(y,Par,n_y,u,Newton_it);

    Tvec<double> grad_Par = grad[0];
    Tvec<double> grad_u = grad[1];

    p_Par = p_Par + epsilon*(grad_Par);
    Tvec<double> p_u_tmp2 = p_u_tmp + (epsilon*grad_u);
    p_u = cos(epsilon/2.0)*p_u_tmp2 - sin(epsilon/2.0)*u;

    //half step for positions
    Par= Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u_tmp2;
  }

  //Evaluate potential and kinetic energies at start and end of trajectory
  double H_current = Hamiltonian_Laplace(Par0,y,n_y,u0,p_Par0,p_u0,M_inv,Newton_it);
  double H_prop = Hamiltonian_Laplace(Par,y,n_y,u,p_Par,p_u,M_inv,Newton_it);

  //Accept or reject the state at end of trajectory
  double accp = exp(H_current-H_prop);
  if (unif < accp)
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par,
      Rcpp::Named("u")  = u,
      Rcpp::Named("acc")  = 1,
      Rcpp::Named("accprob")=accp);
  }
  else
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par0,
      Rcpp::Named("u")  = u0,
      Rcpp::Named("acc")  = 0,
      Rcpp::Named("accprob")=accp);
  }

}

/////////////////////////////

// [[Rcpp::export]]
Rcpp::List HMC_Laplace(int N_mcmc,Eigen::VectorXd Par0, Eigen::VectorXd y, int n_y, Eigen::VectorXd u0,Eigen::MatrixXd p_Par,Eigen::MatrixXd p_u,Eigen::VectorXd epsilons,Eigen::VectorXd Ls,Eigen::VectorXd unif, Eigen::MatrixXd M_inv, int Newton_it)
{
  Tmat<double> Par(N_mcmc+1,4);
  Tmat<double> u((N_mcmc+1),n_y);
  Tvec<double> acc(N_mcmc);
  Tvec<double> accp(N_mcmc);

  Par.row(0)=Par0;
  u.row(0)=u0;

  for (int j=0;j<N_mcmc;j++)
  {
    Rcpp::List hmc_val;
    try {hmc_val=HMC_integrator_Laplace(y,n_y,Par.row(j),u.row(j),p_Par.row(j),p_u.row(j),epsilons(j),Ls(j),unif(j),M_inv,Newton_it);}
    catch(...){hmc_val=Rcpp::List::create(Rcpp::Named("Par")  = Par.row(j),Rcpp::Named("u")  = u.row(j),Rcpp::Named("acc")= 0,Rcpp::Named("accprob")= 0);}


    Tvec<double> Par_tmp = hmc_val[0];
    Par.row(j+1)=Par_tmp;
    Tvec<double> u_tmp = hmc_val[1];
    u.row(j+1)=u_tmp;
    double acc_tmp = hmc_val[2];
    acc(j)= acc_tmp;
    double accp_tmp = hmc_val[3];
    accp(j)= accp_tmp;
  }

  return Rcpp::List::create(
    Rcpp::Named("Par")  = Par,
    Rcpp::Named("u")  = u,
    Rcpp::Named("acc")= acc,
    Rcpp::Named("accprob")= accp);

}
