// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "../../adept.h"
#include "../../adept_source.h"
#include <math.h>
using adept::adouble;
typedef adept::adouble adtype;

template <class T>
using Tmat = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

template <class T>
using Tvec = Eigen::Matrix<T,Eigen::Dynamic,1>;

/////////////////////////////

template <class T>
T Tdnormlog(double y, T mean, T var)
{
  T res = -0.5*((pow(y-mean,2)/var)+log(2*M_PI*var));

  return(res);
}

template <class T>
T TdnormlogX(T y, T mean, T var)
{
  T res = -0.5*((pow(y-mean,2)/var)+log(2*M_PI*var));

  return(res);
}

template <class T>
T TdnormlogX1(T y, double mean, double var)
{
  T res = -0.5*((pow(y-mean,2)/var)+log(2*M_PI*var));

  return(res);
}

/////////////////////////////

template <class T>
T Tdchilog(double Delta, T alpha, T beta, T gamma, T sigxSq, T m, T a1, T a2)
{
  T tmp = 2.0*a2*sigxSq*pow(m,2*gamma)*Delta-1.0;
  double DeltaSq = pow(Delta,2);
  T res =(-0.5*(2.0*DeltaSq*a2*pow(beta,2)*pow(m,2)-4.0*DeltaSq*a2*alpha*beta*m+Delta*pow(m,2*gamma)*pow(a1,2)*sigxSq+2.0*DeltaSq*a2*pow(alpha,2)-4.0*Delta*a2*beta*pow(m,2)-2.0*Delta*a1*beta*m+4.0*Delta*a2*alpha*m+2.0*Delta*a1*alpha+2.0*a2*pow(m,2)+2.0*a1*m)/tmp)-0.5*log(-tmp);
  return(res);
}

/////////////////////////////
template <class T>
T Tdchilog1(double y1, T a1, T a2)
{
  T res = 4.8309615392+((5000.0*pow(y1,2)*a2+0.25*pow(a1,2)+5000.0*a1*y1)/(5000.0-a2))-0.5*log(M_PI*(5000.0-a2));
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
void LinReg(const Tvec<double> &y, Tmat<T> &A, const Tmat<T> &lambdaM, double Delta, T alpha, T beta, T gamma, T sigxSq,T sigySq,int n_z, int n_y)
{
  Tvec<T> yreg(n_y);
  Tmat<T> Xreg(n_z,3);

  Tvec<T> l=lambdaM.col(n_y-1);

  for (int j=0; j<n_z; j++)
  {
    T lj = l.coeff(j);
    yreg.coeffRef(j)=Tdnormlog(y.coeff(n_y-1), lj, sigySq);

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
      yreg(k)=Tdnormlog(y(i), lk, sigySq)+Tdchilog(Delta, alpha, beta, gamma, sigxSq,lk, A(i+1,0),A(i+1,1));

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
T loglike_c (Tvec<double> y, Tvec<T> Par,int n_y, int n_z, Tvec<double> z, int n_reg,Tvec<T> u, int n_u, double Delta){

  T alpha = Par(0);
  T beta = Par(1);
  T gamma = Par(2);
  T sigxSq = exp(Par(3));
  T sigySq = exp(Par(4));

  // Start values for A
  Tmat<T> A(n_y,2);

  for (int i=0; i<n_y; i++)
  {
    A(i,0)=y(i)/sigySq;
    A(i,1)=-1.0/(2*sigySq);
  }

  Tmat<T> lambdaM(n_z,n_y);

  for(int j = 0; j < n_reg; j++)
  {
    T m1_var = 1.0/(10000.0-2.0*A(0,1));
    T m1_mean = m1_var*(A(0,0)+10000.0*y(0));

    // Sample from m
    for (int i=0;i<n_z;i++)
    {
      lambdaM.coeffRef(i,0) = m1_mean+z.coeff(i*n_y)*sqrt(m1_var);

      for (int k=1;k<n_y;k++)
      {
        T m = lambdaM(i,k-1);
        T m_var = -sigxSq*Delta/(2.0*A(k,1)*sigxSq*Delta-pow(m,-2*gamma));
        T tmp = pow(pow(m,gamma),2);
        T m_mean = -(Delta*tmp*A(k,0)*sigxSq+m*(1-Delta*beta)+Delta*alpha)/(2.0*A(k,1)*sigxSq*tmp*Delta-1.0);

        lambdaM.coeffRef(i,k) = m_mean+z.coeff(i*n_y+k)*sqrt(m_var);
      }
    }

      // Regression to update A
      LinReg(y, A, lambdaM, Delta, alpha, beta, gamma, sigxSq,sigySq, n_z, n_y);
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  Tmat<T> lambdaM_u(n_u,n_y);

  T m1_var = 1.0/(10000.0-2.0*A(0,1));
  T m1_mean = m1_var*(A(0,0)+10000.0*y(0));

  for (int i=0;i<n_u;i++)
  {
    lambdaM_u.coeffRef(i,0) = m1_mean+u.coeff(i*n_y)*sqrt(m1_var);

    for (int k=1;k<n_y;k++)
    {
      T m = lambdaM_u(i,k-1);
      T m_var = -sigxSq*Delta/(2.0*A(k,1)*sigxSq*Delta-pow(m,-2*gamma));
      T tmp = pow(pow(m,gamma),2);
      T m_mean = -(Delta*tmp*A(k,0)*sigxSq+m*(1-Delta*beta)+Delta*alpha)/(2.0*A(k,1)*sigxSq*tmp*Delta-1.0);

      lambdaM_u.coeffRef(i,k) = m_mean+u.coeff(i*n_y+k)*sqrt(m_var);
    }
  }

  Tvec<T> logsums(n_u);

  for(int i=0;i<n_u;i++)
  {
    Tvec<T> lw(n_y);

    T l = lambdaM_u.coeff(i,n_y-1);
    lw(n_y-1) = Tdnormlog(y.coeff(n_y-1), l, sigySq) - (A(n_y-1,0)*l+A(n_y-1,1)*pow(l,2));

    for (int j=0;j<(n_y-1);j++)
    {
      T l = lambdaM_u.coeff(i,j);
      lw.coeffRef(j) = Tdnormlog(y.coeff(j), l, sigySq) + Tdchilog(Delta, alpha, beta, gamma, sigxSq, l, A(j+1,0),A(j+1,1)) - (A(j,0)*l+A(j,1)*pow(l,2));
    }
    logsums(i)= lw.sum()+Tdchilog1(y(0), A(0,0),A(0,1));
  }

  T w = logsums.maxCoeff();

  Tvec<T> lsexp(n_u);

  for (int k=0; k<n_u;k++)
  {
    lsexp.coeffRef(k) = exp(logsums.coeff(k) - w);
  }

  T loglike = w + log(lsexp.mean());

  T alpha_prior = -0.5*pow(alpha,2)/1000.0 - 0.5*log(2.0*M_PI*1000.0);
  T beta_prior = -0.5*pow(beta,2)/1000.0 - 0.5*log(2.0*M_PI*1000.0);
  T gamma_prior = 0;
  if (gamma >= 0 && gamma <=4)
  {
    gamma_prior=log(0.25);
  }

  return(loglike+alpha_prior+beta_prior+gamma_prior);
}

/////////////////////////////
/////////////////////////////

Rcpp::List loglike_grad_c(Tvec<double> y,  Tvec<double> Par, int n_y, int n_z, Tvec<double> z, int n_reg, Tvec<double> u, int n_u,double Delta)
{
  adept::Stack stack;

  Tvec<adtype> ad_Par(5);

  for(int k=0;k<5;k++)
  {
    ad_Par.coeffRef(k)=Par.coeff(k);
  }

  Tvec<adtype> ad_u(n_u*n_y);

  for(int i=0;i<(n_u*n_y);i++)
  {
    ad_u.coeffRef(i)=u.coeff(i);
  }

  stack.new_recording();

  adtype res0 = loglike_c<adtype>(y,ad_Par,n_y,n_z,z,n_reg,ad_u,n_u,Delta);
  adtype res = res0/1.0;

  res.set_gradient(1.0);
  stack.compute_adjoint();

  Tvec<double> Pargrads(5);

  for (int q=0;q<5;q++)
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

double Hamiltonian(Tvec<double> Par, Tvec<double> y,int n_y,Tvec<double> u,Tvec<double> z,int n_u, int n_z, int n_reg,Tvec<double> p_Par,Tvec<double> p_u, Tmat<double> M_inv,double Delta)
{
  double U = loglike_c<double>(y,Par,n_y,n_z,z,n_reg,u,n_u,Delta);

  double hamil = -U;

  hamil = hamil + 0.5*p_Par.transpose()*M_inv*p_Par;

  for (int i=0;i<(n_u*n_y);i++)
  {
    hamil = hamil + 0.5*(pow(u.coeff(i),2) + pow(p_u.coeff(i),2));
  }

  return(hamil);
}

/////////////////////////////

Rcpp::List HMC_integrator(Tvec<double> y, int n_y, Tvec<double> z,Tvec<double> Par0,Tvec<double> u0, int n_u, int n_z, Tvec<double> p_Par0,Tvec<double> p_u0,double epsilon,int L, double unif, int n_reg, Tmat<double> M_inv,double Delta)
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
    grad = loglike_grad_c(y,Par,n_y,n_z,z,n_reg,u,n_u,Delta);

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
  double H_current = Hamiltonian(Par0,y,n_y,u0,z,n_u,n_z,n_reg,p_Par0,p_u0,M_inv,Delta);
  double H_prop = Hamiltonian(Par,y,n_y,u,z,n_u,n_z,n_reg,p_Par,p_u,M_inv,Delta);

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
Rcpp::List HMC_EIS(int N_mcmc,Eigen::VectorXd Par0, Eigen::VectorXd y, int n_y, Eigen::VectorXd u0,Eigen::MatrixXd z,int n_u, int n_z, Eigen::MatrixXd p_Par,Eigen::MatrixXd p_u,Eigen::VectorXd epsilons,Eigen::VectorXd Ls,Eigen::VectorXd unif,int n_reg, Eigen::MatrixXd M_inv, double Delta)
{
  Tmat<double> Par(N_mcmc+1,5);
  Tmat<double> u((N_mcmc+1),(n_u*n_y));
  Tvec<double> acc(N_mcmc);
  Tvec<double> accp(N_mcmc);

  Par.row(0)=Par0;
  u.row(0)=u0;

  for (int j=0;j<N_mcmc;j++)
  {
    Rcpp::List hmc_val;
    try {hmc_val=HMC_integrator(y,n_y,z.row(j),Par.row(j),u.row(j),n_u,n_z,p_Par.row(j),p_u.row(j),epsilons(j),Ls(j),unif(j),n_reg,M_inv,Delta);}
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
* diag(G) = diagElem (diagElem has length n)
* first sup/sub-diagonal elements are in offDiagElem ( offDiagElem has length n-1)

Th routine returns a vector[2*n] L where
* L[1:n] is the diagonal of L
* L[n+1:2*n-1] is the sub-diagonal of L
* L[2*n] is the log-determinant of L
*/
template <class T>
Tvec<T> CIP_TriDiagChol(Tvec<T> diagElem, Tvec<T> offDiagElem, int N)
{
  Tvec<T> L(2*N);
  T LlogDet;
  L(0) = sqrt(diagElem(0));
  LlogDet = log(L(0));

  for(int t=1;t<N;t++)
  {
    L(N+t-1)= offDiagElem(t-1)/L(t-1);
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

/////////////////////////////////
/////////////////////////////////


template <class T>
T loglike_Laplace(Tvec<double> y,Tvec<double> ySq,int n_y, Tvec<T> Par,Tvec<T> u, double Delta, int Newton_it)
{
  T alpha = Par(0);
  T beta = Par(1);
  T gamma = Par(2);
  T sxsq = exp(Par(3));
  T sysq = exp(Par(4));

  T gammaSq = pow(gamma,2);
  T Delbeta = Delta*beta-1.0;
  T DelbetaSq = pow(Delta*beta-1.0,2);
  Tvec<T> tmp(n_y-1);
  Tvec<T> tmp2(n_y-1);
  Tvec<T> diag(n_y);
  Tvec<T> odiag(n_y-1);
  Tvec<T> grad(n_y);
  Tvec<T> L(2*n_y);
  Tvec<T> x(n_y);
  Tvec<T> h(n_y);
  Tvec<T> hSq(n_y);


  /////////////////////////////////////////////////////////
  // Build Hessian

  for(int i=0;i<n_y;i++)
  {
    h(i) = y(i);
    hSq(i) = ySq(i);
  }

  // Newton iterations

  for(int w=0;w<Newton_it;w++)
  {
	  for(int i=0;i<(n_y-1);i++)
    {
      tmp(i) = h(i+1) + Delbeta*h(i) - Delta*alpha;
    }

    for(int j=0;j<(n_y-1);j++)
    {
      tmp2(j) = sxsq*Delta*pow(h(j),2.0*gamma+2.0);
    }

    diag(0) = 10000.0 + 1.0/sysq -(gamma / hSq(0)) +(DelbetaSq*hSq(0)-4.0*tmp(0)*Delbeta*gamma*h(0)+pow(tmp(0),2)*(2.0*gammaSq+gamma))/tmp2(0);

    for(int k=1;k<(n_y-1);k++)
    {
      diag(k) = 1.0/sysq -(gamma / hSq(k))+ hSq(k-1)/tmp2(k-1) +(DelbetaSq*hSq(k)-4.0*tmp(k)*Delbeta*gamma*h(k)+pow(tmp(k),2)*(2.0*gammaSq+gamma))/tmp2(k);
    }

    diag(n_y-1) = 1.0/sysq + hSq(n_y-2)/tmp2(n_y-2);

    // First off-diagonal
    for(int q=0;q<(n_y-1);q++)
    {
      odiag(q) = (Delbeta*hSq(q)/tmp2(q))-(2.0*tmp(q)*h(q)*gamma)/tmp2(q);
    }

    // Cholesky
    L = CIP_TriDiagChol(diag,odiag,n_y);

    grad(0) = -10000.0*h(0) + 10000.0*y(0) + (pow(tmp(0),2)*gamma*h(0)-tmp(0)*Delbeta*hSq(0))/tmp2(0) - (gamma/h(0)) + ((y(0)-h(0))/sysq);

    for(int k=1;k<(n_y-1);k++)
    {
      grad(k) = -(tmp(k-1)*hSq(k-1))/tmp2(k-1) + (pow(tmp(k),2)*gamma*h(k)-tmp(k)*Delbeta*hSq(k))/tmp2(k) - (gamma/h(k)) + ((y(k)-h(k))/sysq);
    }

    grad(n_y-1) = -(tmp(n_y-2)*hSq(n_y-2))/tmp2(n_y-2) + ((y(n_y-1)-h(n_y-1))/sysq);

    // update eval point
    h = h + CIP_TriDiagChol_LLT_solve(L,grad,n_y);

    for(int j=0;j<n_y;j++)
    {
      hSq(j) = pow(h(j),2);
    }
  }

  // Transformation
  x = h + CIP_TriDiagChol_LT_solve(L,u,n_y);

  //////////////////////////////////////////////////////////

  T loglike= -5000.0*pow(x(0)-y(0),2) - 0.5*log(M_PI/5000.0) -0.5*pow(y(0)-x(0),2)/sysq - 0.5*log(2.0*M_PI*sysq);
  loglike += 0.5*pow(u(0),2);

  for (int j=1;j<n_y;j++)
  {
    T var = sxsq*Delta*pow(x(j-1),2.0*gamma);
    loglike += -0.5*pow(x(j)-x(j-1)-Delta*(alpha-beta*x(j-1)),2)/var - 0.5*log(2.0*M_PI*var) -0.5*pow(y(j)-x(j),2)/sysq - 0.5*log(2.0*M_PI*sysq);
    loglike += 0.5*pow(u(j),2);
  }
  loglike += -L(2*n_y-1);

  T alpha_prior = -0.5*pow(alpha,2)/1000.0 - 0.5*log(2.0*M_PI*1000.0);
  T beta_prior = -0.5*pow(beta,2)/1000.0 - 0.5*log(2.0*M_PI*1000.0);
  T gamma_prior = 0;
  if (gamma >= 0 && gamma <=4)
  {
    gamma_prior=log(0.25);
  }

  return(loglike+alpha_prior+beta_prior+gamma_prior);
}

/////////////////////////////
/////////////////////////////
Rcpp::List loglike_Laplace_grad(Tvec<double> y,Tvec<double> ySq, int n_y, Tvec<double> Par,Tvec<double> u, double Delta, int Newton_it)
{
  adept::Stack stack;

  Tvec<adtype> ad_Par(5);

  for(int k=0;k<5;k++)
  {
    ad_Par.coeffRef(k)=Par.coeff(k);
  }

  Tvec<adtype> ad_u(n_y);

  for(int i=0;i<n_y;i++)
  {
    ad_u.coeffRef(i)=u.coeff(i);
  }

  stack.new_recording();

  adtype res0 = loglike_Laplace<adtype>(y,ySq,n_y,ad_Par,ad_u,Delta,Newton_it);
  adtype res = res0/1.0;

  res.set_gradient(1.0);
  stack.compute_adjoint();

  Tvec<double> Pargrads(5);

  for (int q=0;q<5;q++)
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
double Hamiltonian_Laplace(Tvec<double> Par, Tvec<double> y,Tvec<double> ySq,int n_y,Tvec<double> u, Tvec<double> p_Par,Tvec<double> p_u, Tmat<double> M_inv, double Delta, int Newton_it)
{
  double U = loglike_Laplace<double>(y,ySq,n_y,Par,u,Delta,Newton_it);

  double hamil = -U;

  hamil = hamil + 0.5*p_Par.transpose()*M_inv*p_Par;

  for (int i=0;i<n_y;i++)
  {
    hamil = hamil + 0.5*(pow(u.coeff(i),2) + pow(p_u.coeff(i),2));
  }

  return(hamil);
}

/////////////////////////////

Rcpp::List HMC_integrator_Laplace(Tvec<double> y, Tvec<double> ySq, int n_y, Tvec<double> Par0,Tvec<double> u0, Tvec<double> p_Par0,Tvec<double> p_u0,double epsilon,int L, double unif, Tmat<double> M_inv,double Delta, int Newton_it)
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

    //full step for momentum
    grad = loglike_Laplace_grad(y,ySq,n_y,Par,u,Delta,Newton_it);

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
  double H_current = Hamiltonian_Laplace(Par0,y,ySq,n_y,u0,p_Par0,p_u0,M_inv,Delta,Newton_it);
  double H_prop = Hamiltonian_Laplace(Par,y,ySq,n_y,u,p_Par,p_u,M_inv,Delta,Newton_it);

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
Rcpp::List HMC_Laplace(int N_mcmc,Eigen::VectorXd Par0, Eigen::VectorXd y, Eigen::VectorXd ySq, int n_y, Eigen::VectorXd u0, Eigen::MatrixXd p_Par,Eigen::MatrixXd p_u,Eigen::VectorXd epsilons,Eigen::VectorXd Ls,Eigen::VectorXd unif, Eigen::MatrixXd M_inv, double Delta, int Newton_it)
{
  Tmat<double> Par(N_mcmc+1,5);
  Tmat<double> u((N_mcmc+1),n_y);
  Tvec<double> acc(N_mcmc);
  Tvec<double> accp(N_mcmc);

  Par.row(0)=Par0;
  u.row(0)=u0;

  for (int j=0;j<N_mcmc;j++)
  {
    Rcpp::List hmc_val;
    try {hmc_val=HMC_integrator_Laplace(y,ySq,n_y,Par.row(j),u.row(j),p_Par.row(j),p_u.row(j),epsilons(j),Ls(j),unif(j),M_inv,Delta,Newton_it);}
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
