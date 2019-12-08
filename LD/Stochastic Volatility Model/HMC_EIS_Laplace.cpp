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
T Tdnormlog(double r, T var)
{
  T res = -0.5*((pow(r,2)/var)+log(2*M_PI*var));

  return(res);
}

/////////////////////////////

template <class T>
T Tdchilog(T tmp1, T tmp2, T tmp3, T tmp4, T tmp5, T m)
{
  T res=0.5*((m*tmp1+pow(m,2)*tmp2+tmp3)/tmp4)-tmp5;
  return(res);
}

/////////////////////////////
template <class T>
T Tdchilog1(T vSq, T delta, T gamma, Tvec<T> A)
{
  T res=(0.5*(-pow(A(0),2)*delta*vSq+pow(A(0),2)*vSq-2*A(0)*pow(delta,2)*gamma+2*A(1)*delta*pow(gamma,2)+2*A(1)*pow(gamma,2)+2*A(0)*gamma)/((-1+delta)*(2*A(1)*vSq+pow(delta,2)-1)))+0.5*log(1.0-pow(delta,2))-0.5*log(1.0-2.0*A(1)*vSq-pow(delta,2));
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
void LinReg(const Tvec<double> &y, Tmat<T> &A, const Tmat<T> &lambdaM, T vSq, T delta, T gamma, int n_z, int n_y)
{
  Tvec<T> yreg(n_y);
  Tmat<T> Xreg(n_z,3);

  Tvec<T> l=lambdaM.col(n_y-1);

  for (int j=0; j<n_z; j++)
  {
	T lj = l.coeff(j);
    T var = exp(lj);
    yreg.coeffRef(j)=Tdnormlog(y.coeff(n_y-1), var);

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
	  Tvec<T> Arow = A.row(i+1);

	  T tmp1 = 4*Arow(1)*delta*gamma + 2*Arow(0)*delta;
	  T tmp2 = 2*Arow(1)*pow(delta,2);
	  T tmp3 = pow(Arow(0),2)*vSq+2*Arow(1)*pow(gamma,2)+2*Arow(0)*gamma;
	  T tmp4 = 1-2*Arow(1)*vSq;
	  T tmp5 = 0.5*log(tmp4);

    for (int k=0; k<n_z; k++)
    {
	  T lk = l.coeff(k);
      T var = exp(lk);

      yreg.coeffRef(k)=Tdnormlog(y.coeff(i), var)+Tdchilog(tmp1,tmp2,tmp3,tmp4,tmp5,lk);

      Xreg.coeffRef(k,0)=1;
      Xreg.coeffRef(k,1)=lk;
      Xreg.coeffRef(k,2)=pow(lk,2);
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
T loglike_c (Tvec<double> y, Tvec<T> Par,int n_y, int n_z, Tvec<double> z, int n_yeg,Tvec<T> u, int n_u, double s0, double df, double alpha, double beta){

  T gamma = Par(0);
  T delta = tanh(Par(1));
  T vSq = exp(Par(2));

  // Start values for A
  Tmat<T> A(n_y,2);

  for (int i=0; i<n_y; i++)
  {
    A.coeffRef(i,0)=0.5*log(pow(y.coeff(i),2));
    A.coeffRef(i,1)=-0.25;
  }

  Tmat<T> lambdaM(n_z,n_y);

  for(int j = 0; j < n_yeg; j++)
  {
    T m1_var = vSq/(1.0-pow(delta,2)-2.0*vSq*A.coeff(0,1));
    T m1_mean = m1_var*(A.coeff(0,0)+((gamma+delta*gamma)/vSq));

    for (int i=0;i<n_z;i++)
    {
      lambdaM.coeffRef(i,0) = m1_mean+z.coeff(i)*sqrt(m1_var);
    }

    for (int k=1;k<n_y;k++)
    {
      T m_var = vSq/(1.0-2.0*vSq*A.coeff(k,1));

      for (int i=0;i<n_z;i++)
      {
        T m_mean = m_var*(A.coeff(k,0)+(gamma+delta*lambdaM.coeff(i,k-1))/vSq);
        lambdaM.coeffRef(i,k) = m_mean+z.coeff(i+k*n_z)*sqrt(m_var);
      }
    }

      // Regression to update A
      LinReg(y, A, lambdaM, vSq, delta, gamma, n_z, n_y);
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  Tmat<T> lambdaM_u(n_u,n_y);

  T m1_var = vSq/(1.0-pow(delta,2)-2.0*vSq*A.coeff(0,1));
  T m1_mean = m1_var*(A.coeff(0,0)+((gamma+delta*gamma)/vSq));

  for (int i=0;i<n_u;i++)
  {
    lambdaM_u.coeffRef(i,0) = m1_mean+u.coeff(i)*sqrt(m1_var);
  }

  for (int k=1;k<n_y;k++)
  {
    T m_var = vSq/(1.0-2.0*vSq*A.coeff(k,1));

    for (int i=0;i<n_u;i++)
    {
      T m_mean = m_var*(A.coeff(k,0)+(gamma+delta*lambdaM_u.coeff(i,k-1))/vSq);
      lambdaM_u.coeffRef(i,k) = m_mean+u.coeff(i+k*n_u)*sqrt(m_var);
    }
  }

  Tvec<T> logsums(n_u);
  logsums.setZero();

  Tvec<T> Arow0 = A.row(0);
  T chilog1 = Tdchilog1(vSq, delta, gamma, Arow0);

  for(int i=0;i<n_u;i++)
  {
    T l = lambdaM_u.coeff(i,n_y-1);
    T var = exp(l);
    logsums.coeffRef(i)=logsums.coeff(i)+Tdnormlog(y.coeff(n_y-1), var) - (A.coeff(n_y-1,0)*l+A.coeff(n_y-1,1)*pow(l,2))+chilog1;
  }

  for (int j=0;j<(n_y-1);j++)
  {
    Tvec<T> Arow = A.row(j+1);

    T tmp1 = 4*Arow(1)*delta*gamma + 2*Arow(0)*delta;
    T tmp2 = 2*Arow(1)*pow(delta,2);
    T tmp3 = pow(Arow(0),2)*vSq+2*Arow(1)*pow(gamma,2)+2*Arow(0)*gamma;
    T tmp4 = 1-2*Arow(1)*vSq;
    T tmp5 = 0.5*log(tmp4);

    for(int i=0;i<n_u;i++)
    {
      T l = lambdaM_u.coeff(i,j);
      T var = exp(l);
      logsums.coeffRef(i)=logsums.coeff(i)+Tdnormlog(y.coeff(j), var) + Tdchilog(tmp1,tmp2,tmp3,tmp4,tmp5,l) - (A.coeff(j,0)*l+A.coeff(j,1)*pow(l,2));
    }
  }

  ////////////////////

  T w = logsums.maxCoeff();

  Tvec<T> lsexp(n_u);

  for (int k=0; k<n_u;k++)
  {
    lsexp.coeffRef(k) = exp(logsums.coeff(k) - w);
  }

  T loglike = w + log(lsexp.mean());

  T v_prior = ((-0.5*df*s0)/vSq)-0.5*df*log(vSq);
  T delta_prior = log(0.5*delta+0.5)*(alpha-1)+log(0.5-0.5*delta)*(beta-1)+log(0.5-0.5*pow(delta,2));

  return(loglike+v_prior+delta_prior);
}

/////////////////////////////
/////////////////////////////

Tvec<double> loglike_grad_c(Tvec<double> y,  Tvec<double> Par, int n_y, int n_z, Tvec<double> z, int n_yeg, Tvec<double> u, int n_u,double s0, double df, double alpha, double beta)
{
  adept::Stack stack;

  Tvec<adtype> ad_Par(3);

  for(int k=0;k<3;k++)
  {
    ad_Par.coeffRef(k)=Par.coeff(k);
  }

  Tvec<adtype> ad_u(n_u*n_y);

  for(int i=0;i<(n_u*n_y);i++)
  {
    ad_u.coeffRef(i)=u.coeff(i);
  }

  stack.new_recording();

  adtype res0 = loglike_c<adtype>(y,ad_Par,n_y,n_z,z,n_yeg,ad_u,n_u,s0,df,alpha,beta);
  adtype res = res0/1.0;

  res.set_gradient(1.0);
  stack.compute_adjoint();

  Tvec<double> grads(3+n_u*n_y);

  for (int q=0;q<3;q++)
  {
    grads.coeffRef(q)= ad_Par.coeff(q).get_gradient();
  }

  for (int j=0;j<(n_u*n_y);j++)
  {
    grads.coeffRef(j+3)= ad_u.coeff(j).get_gradient();
  }

  return(grads);
}

/////////////////////////////

double Hamiltonian(Tvec<double> Par, Tvec<double> y,int n_y,Tvec<double> u,Tvec<double> z,int n_u, int n_z, int n_yeg,Tvec<double> p_Par,Tvec<double> p_u, Tmat<double> M_inv,double s0, double df, double alpha, double beta)
{
  double U = loglike_c<double>(y,Par,n_y,n_z,z,n_yeg,u,n_u,s0,df,alpha,beta);

  double hamil = -U;

  hamil = hamil + 0.5*p_Par.transpose()*M_inv*p_Par;

  for (int i=0;i<(n_u*n_y);i++)
  {
    hamil = hamil + 0.5*(pow(u.coeff(i),2) + pow(p_u.coeff(i),2));
  }

  return(hamil);
}

/////////////////////////////

Rcpp::List HMC_integrator(Tvec<double> y, int n_y, Tvec<double> z,Tvec<double> Par0,Tvec<double> u0, int n_u, int n_z, Tvec<double> p_Par0,Tvec<double> p_u0,double epsilon,int L, double unif, int n_yeg, Tmat<double> M_inv,double s0, double df, double alpha, double beta)
{
  Tvec<double> Par = Par0;
  Tvec<double> u = u0;
  Tvec<double> p_Par = p_Par0;
  Tvec<double> p_u = p_u0;

  for(int i=0;i<L;i++)
  {
    Tvec<double> p_u_tmp = cos(epsilon/2.0)*p_u - sin(epsilon/2.0)*u;

    //half step for positions
    Par = Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u;
    //
    //full step for momentum
    Tvec<double> grad = loglike_grad_c(y,Par,n_y,n_z,z,n_yeg,u,n_u,s0,df,alpha,beta);

    Tvec<double> grad_Par = grad.segment(0,3);
    Tvec<double> grad_u = grad.segment(3,n_u*n_y);

    p_Par = p_Par + epsilon*(grad_Par);
    Tvec<double> p_u_tmp2 = p_u_tmp + (epsilon*grad_u);
    p_u = cos(epsilon/2.0)*p_u_tmp2 - sin(epsilon/2.0)*u;

    //half step for positions
    Par= Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u_tmp2;
  }

  //Evaluate potential and kinetic energies at start and end of trajectory
  double H_current = Hamiltonian(Par0,y,n_y,u0,z,n_u,n_z,n_yeg,p_Par0,p_u0,M_inv,s0,df,alpha,beta);
  double H_prop = Hamiltonian(Par,y,n_y,u,z,n_u,n_z,n_yeg,p_Par,p_u,M_inv,s0,df,alpha,beta);

  //Accept or reject the state at end of trajectory
  double accp = exp(H_current-H_prop);

  if (unif < accp)
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par,
      Rcpp::Named("u")  = u,
      Rcpp::Named("acc")  = 1,
      Rcpp::Named("accprob")=accp,
      Rcpp::Named("Par1")=Par,
      Rcpp::Named("u1")=u,
      Rcpp::Named("p_Par1")=p_Par,
      Rcpp::Named("p_u1")=p_u);
  }
  else
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par0,
      Rcpp::Named("u")  = u0,
      Rcpp::Named("acc")  = 0,
      Rcpp::Named("accprob")=accp,
      Rcpp::Named("Par1")=Par,
      Rcpp::Named("u1")=u,
      Rcpp::Named("p_Par1")=p_Par,
      Rcpp::Named("p_u1")=p_u);
  }

}

/////////////////////////////
Rcpp::List HMC_integrator_leapfrog(Tvec<double> y, int n_y, Tvec<double> z,Tvec<double> Par0,Tvec<double> u0, int n_u, int n_z, Tvec<double> p_Par0,Tvec<double> p_u0,double epsilon,int L, double unif, int n_yeg, Tmat<double> M_inv,double s0, double df, double alpha, double beta)
{
  Tvec<double> Par = Par0;
  Tvec<double> u = u0;
  Tvec<double> p_Par = p_Par0;
  Tvec<double> p_u = p_u0;

  Tvec<double> grad = loglike_grad_c(y,Par,n_y,n_z,z,n_yeg,u,n_u,s0,df,alpha,beta);
  Tvec<double> grad_Par = grad.segment(0,3);
  Tvec<double> grad_u = grad.segment(3,n_u*n_y);

  Rcpp::Rcout << "grad_Par = " << grad_Par << std::endl;

  for(int i=0;i<L;i++)
  {
    //half step for momentums
    p_Par = p_Par + (epsilon/2.0)*(grad_Par);
    p_u = p_u + (epsilon/2.0)*(grad_u);

    //full step for positions
    Par = Par + epsilon*M_inv*p_Par;
    u = u + epsilon*p_u;

    Rcpp::Rcout << "Par = " << Par <<  std::endl;

    //half step for momentums
    Tvec<double> grad = loglike_grad_c(y,Par,n_y,n_z,z,n_yeg,u,n_u,s0,df,alpha,beta);
    Tvec<double> grad_Par = grad.segment(0,3);
    Tvec<double> grad_u = grad.segment(3,n_u*n_y);

    p_Par = p_Par + (epsilon/2.0)*(grad_Par);
    p_u = p_u + (epsilon/2.0)*(grad_u);
  }

  Rcpp::Rcout << "grad_Par = " << grad_Par << std::endl;

  //Evaluate potential and kinetic energies at start and end of trajectory
  double H_current = Hamiltonian(Par0,y,n_y,u0,z,n_u,n_z,n_yeg,p_Par0,p_u0,M_inv,s0,df,alpha,beta);
  double H_prop = Hamiltonian(Par,y,n_y,u,z,n_u,n_z,n_yeg,p_Par,p_u,M_inv,s0,df,alpha,beta);

  //Accept or reject the state at end of trajectory
  double accp = exp(H_current-H_prop);

  if (unif < accp)
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par,
      Rcpp::Named("u")  = u,
      Rcpp::Named("acc")  = 1,
      Rcpp::Named("accprob")=accp,
      Rcpp::Named("Par1")=Par,
      Rcpp::Named("u1")=u,
      Rcpp::Named("p_Par1")=p_Par,
      Rcpp::Named("p_u1")=p_u);
  }
  else
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par0,
      Rcpp::Named("u")  = u0,
      Rcpp::Named("acc")  = 0,
      Rcpp::Named("accprob")=accp,
      Rcpp::Named("Par1")=Par,
      Rcpp::Named("u1")=u,
      Rcpp::Named("p_Par1")=p_Par,
      Rcpp::Named("p_u1")=p_u);
  }

}


// [[Rcpp::export]]
Rcpp::List HMC_EIS(int N_mcmc,Eigen::VectorXd Par0, Eigen::VectorXd y, int n_y, Eigen::VectorXd u0,Eigen::MatrixXd z,int n_u, int n_z, Eigen::MatrixXd p_Par,Eigen::MatrixXd p_u,Eigen::VectorXd epsilons,Eigen::VectorXd Ls,Eigen::VectorXd unif,int n_yeg, Eigen::MatrixXd M_inv, double s0, double df, double alpha, double beta)
{
  Tmat<double> Par(N_mcmc+1,3);
  Tmat<double> u((N_mcmc+1),(n_u*n_y));
  Tvec<double> acc(N_mcmc);
  Tvec<double> accp(N_mcmc);

  Par.row(0)=Par0;
  u.row(0)=u0;

  for (int j=0;j<N_mcmc;j++)
  {
    Rcpp::List hmc_val;
    try {hmc_val=HMC_integrator(y,n_y,z.row(j),Par.row(j),u.row(j),n_u,n_z,p_Par.row(j),p_u.row(j),epsilons(j),Ls(j),unif(j),n_yeg,M_inv,s0,df,alpha,beta);}
    catch(...){hmc_val=Rcpp::List::create(Rcpp::Named("Par")  = Par.row(j),Rcpp::Named("u")  = u.row(j),Rcpp::Named("acc")= 0,Rcpp::Named("accprob")= -1);}

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
T loglike_laplace(Tvec<double> y,Tvec<double> ySq, Tvec<T> Par, int n_y, Tvec<T> u, double s0, double df, double alpha, double beta,int Newton_it)
{
  T gamma = Par(0);
  T delta = tanh(Par(1));
  T vSq = exp(Par(2));

  T deltaSq = pow(delta,2);
  Tvec<T> x(n_y);
  Tvec<T> grad(n_y);
  Tvec<T> diag(n_y);
  Tvec<T> L(2*n_y);
  Tvec<T> GxGy(n_y);
  Tvec<T> h(n_y);

  /////////////////////////////////////////////////////////
  // Build Hessian
  T firstlast = 0.5+(1.0/vSq);
  T diagval = 0.5+(1.0+deltaSq)/vSq;
  T offdiag = -delta/vSq;

  L = CIP_TriDiagChol_const1n(n_y,firstlast,diagval,offdiag);

  // Build mean

  GxGy(0) = log(fabs(y(0)))+(gamma/vSq);
  for(int i=1;i<(n_y-1);i++)
  {
    GxGy(i) = log(fabs(y(i))) + gamma*(1.0-delta)/vSq;
  }
  GxGy(n_y-1) = log(fabs(y(n_y-1)))+(gamma/vSq);

  h = CIP_TriDiagChol_LLT_solve(L,GxGy,n_y);

  // Newton iterations

	for(int w=0;w<Newton_it;w++)
	{
		diag(0) = 0.5*ySq(0)*exp(-h(0)) + 1.0/vSq;

		for(int t=1;t<(n_y-1);t++)
		{
			diag(t) = 0.5*ySq(t)*exp(-h(t)) + (1.0+deltaSq)/vSq;
		}

		diag(n_y-1) = 0.5*ySq(n_y-1)*exp(-h(n_y-1)) + 1.0/vSq;

		L = CIP_TriDiagChol_diag_constod(diag,offdiag,n_y);

		grad(0) = 0.5*ySq(0)*exp(-h(0)) - 0.5 + (delta*h(1)+gamma-h(0))/vSq;

		for(int t=1;t<(n_y-1);t++)
		{
			grad(t) = 0.5*ySq(t)*exp(-h(t)) - 0.5 + (delta*(h(t-1)+h(t+1)-delta*h(t)-gamma)+gamma-h(t))/vSq;
		}

		grad(n_y-1) = 0.5*ySq(n_y-1)*exp(-h(n_y-1)) - 0.5 + (delta*h(n_y-2)+gamma-h(n_y-1))/vSq;

		h = h + CIP_TriDiagChol_LLT_solve(L,grad,n_y);
	}

  // Transformation
  x = h+CIP_TriDiagChol_LT_solve(L,u,n_y);

  //////////////////////////////////////////////////////////

	T var = vSq/(1-deltaSq);
	T loglike=-0.5*pow(x(0)-gamma/(1-delta),2)/var - 0.5*log(2.0*M_PI*var);
	var = exp(x(0));
	loglike += -0.5*pow(y(0),2)/var - 0.5*log(2.0*M_PI*var);
	loglike += 0.5*pow(u(0),2);

  for (int j=1;j<n_y;j++)
  {
    T var = exp(x(j));
    loglike += -0.5*pow(x(j)-gamma-delta*x(j-1),2)/vSq - 0.5*log(2.0*M_PI*vSq);
    loglike += -0.5*pow(y(j),2)/var - 0.5*log(2.0*M_PI*var);
    loglike += 0.5*pow(u(j),2);
  }
	loglike += -L(2*n_y-1);


  T v_prior = ((-0.5*df*s0)/vSq)-0.5*df*log(vSq);
  T delta_prior = log(0.5*delta+0.5)*(alpha-1)+log(0.5-0.5*delta)*(beta-1)+log(0.5-0.5*pow(delta,2));

  return(loglike+v_prior+delta_prior);
}

/////////////////////////////
/////////////////////////////

Tvec<double> loglike_grad_laplace(Tvec<double> y, Tvec<double> ySq,  Tvec<double> Par, int n_y, Tvec<double> u, double s0, double df, double alpha, double beta, int Newton_it)
{
  adept::Stack stack;

  Tvec<adtype> ad_Par(3);

  for(int k=0;k<3;k++)
  {
    ad_Par.coeffRef(k)=Par.coeff(k);
  }

  Tvec<adtype> ad_u(n_y);

  for(int i=0;i<(n_y);i++)
  {
    ad_u.coeffRef(i)=u.coeff(i);
  }

  stack.new_recording();

  adtype res0 = loglike_laplace<adtype>(y,ySq,ad_Par,n_y,ad_u,s0,df,alpha,beta,Newton_it);
  adtype res = res0/1.0;

  res.set_gradient(1.0);
  stack.compute_adjoint();

  Tvec<double> grads(3+n_y);

  for (int q=0;q<3;q++)
  {
    grads.coeffRef(q)= ad_Par.coeff(q).get_gradient();
  }

  for (int j=0;j<(n_y);j++)
  {
    grads.coeffRef(j+3)= ad_u.coeff(j).get_gradient();
  }

  return(grads);
}

/////////////////////////////

Rcpp::List Hamiltonian_laplace(Tvec<double> Par, Tvec<double> y,Tvec<double> ySq,int n_y,Tvec<double> u,Tvec<double> p_Par,Tvec<double> p_u, Tmat<double> M_inv,double s0, double df, double alpha, double beta, int Newton_it)
{
  double U = loglike_laplace<double>(y,ySq,Par,n_y,u,s0,df,alpha,beta,Newton_it);

  double hamil = -U;

  hamil = hamil + 0.5*p_Par.transpose()*M_inv*p_Par;

  for (int i=0;i<(n_y);i++)
  {
    hamil = hamil + 0.5*(pow(u.coeff(i),2) + pow(p_u.coeff(i),2));
  }

  return Rcpp::List::create(
    Rcpp::Named("Loglike")  = U,
    Rcpp::Named("Hamiltonian")  = hamil);
}

/////////////////////////////

Rcpp::List HMC_integrator_laplace(Tvec<double> y, Tvec<double> ySq, int n_y, Tvec<double> Par0,Tvec<double> u0, Tvec<double> p_Par0,Tvec<double> p_u0,double epsilon,int L, double unif, Tmat<double> M_inv,double s0, double df, double alpha, double beta, int Newton_it)
{
  Tvec<double> Par = Par0;
  Tvec<double> u = u0;
  Tvec<double> p_Par = p_Par0;
  Tvec<double> p_u = p_u0;

  for(int i=0;i<L;i++)
  {
    Tvec<double> p_u_tmp = cos(epsilon/2.0)*p_u - sin(epsilon/2.0)*u;

    //half step for positions
    Par = Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u;
    //
    //full step for momentum
    Tvec<double> grad = loglike_grad_laplace(y,ySq,Par,n_y,u,s0,df,alpha,beta,Newton_it);

    Tvec<double> grad_Par = grad.segment(0,3);
    Tvec<double> grad_u = grad.segment(3,n_y);

    p_Par = p_Par + epsilon*(grad_Par);
    Tvec<double> p_u_tmp2 = p_u_tmp + (epsilon*grad_u);
    p_u = cos(epsilon/2.0)*p_u_tmp2 - sin(epsilon/2.0)*u;

    //half step for positions
    Par= Par + (epsilon/2.0)*M_inv*p_Par;
    u = cos(epsilon/2.0)*u + sin(epsilon/2.0)*p_u_tmp2;
  }

  //Evaluate potential and kinetic energies at start and end of trajectory
  Rcpp::List H_current = Hamiltonian_laplace(Par0,y,ySq,n_y,u0,p_Par0,p_u0,M_inv,s0,df,alpha,beta,Newton_it);
  Rcpp::List H_prop = Hamiltonian_laplace(Par,y,ySq,n_y,u,p_Par,p_u,M_inv,s0,df,alpha,beta,Newton_it);

  double h1 = H_current[1];
  double h2 = H_prop[1];
  double ll = H_current[0];
  double llp = H_prop[0];
  //Accept or reject the state at end of trajectory
  double accp = exp(h1-h2);

  if (unif < accp)
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par,
      Rcpp::Named("u")  = u,
      Rcpp::Named("acc")  = 1,
      Rcpp::Named("accprob")=accp,
      Rcpp::Named("loglike")=ll,
      Rcpp::Named("loglike_prop")=llp);
  }
  else
  {
    return Rcpp::List::create(
      Rcpp::Named("Par")  = Par0,
      Rcpp::Named("u")  = u0,
      Rcpp::Named("acc")  = 0,
      Rcpp::Named("accprob")=accp,
      Rcpp::Named("loglike")=ll,
      Rcpp::Named("loglike_prop")=llp);
  }

}

/////////////////////////////

// [[Rcpp::export]]
Rcpp::List HMC_Laplace(int N_mcmc,Eigen::VectorXd Par0, Eigen::VectorXd y, Eigen::VectorXd ySq, int n_y, Eigen::VectorXd u0, Eigen::MatrixXd p_Par,Eigen::MatrixXd p_u,Eigen::VectorXd epsilons,Eigen::VectorXd Ls,Eigen::VectorXd unif, Eigen::MatrixXd M_inv, double s0, double df, double alpha, double beta, int Newton_it)
{
  Tmat<double> Par(N_mcmc+1,3);
  Tvec<double> u(n_y);
  Tvec<double> acc(N_mcmc);
  Tvec<double> accp(N_mcmc);
  Tvec<double> loglike(N_mcmc);
  Tvec<double> loglike_prop(N_mcmc);

  Par.row(0)=Par0;
  u=u0;

  for (int j=0;j<N_mcmc;j++)
  {
    Rcpp::List hmc_val;
    try {hmc_val=HMC_integrator_laplace(y,ySq,n_y,Par.row(j),u,p_Par.row(j),p_u.row(j),epsilons(j),Ls(j),unif(j),M_inv,s0,df,alpha,beta,Newton_it);}
    catch(...){hmc_val=Rcpp::List::create(Rcpp::Named("Par")  = Par.row(j),Rcpp::Named("u")  = u,Rcpp::Named("acc")= 0,Rcpp::Named("accprob")= -1);}

    Tvec<double> Par_tmp = hmc_val[0];
    Par.row(j+1)=Par_tmp;
    Tvec<double> u_tmp = hmc_val[1];
    u=u_tmp;
    double acc_tmp = hmc_val[2];
    acc(j)= acc_tmp;
    double accp_tmp = hmc_val[3];
    accp(j)= accp_tmp;
    double ll = hmc_val[4];
    loglike(j)= ll;
    double llp = hmc_val[5];
    loglike_prop(j)= llp;
  }

  return Rcpp::List::create(
    Rcpp::Named("Par")  = Par,
    Rcpp::Named("acc")= acc,
    Rcpp::Named("accprob")= accp,
    Rcpp::Named("loglike")=loglike,
    Rcpp::Named("loglike_prop")=loglike_prop);

}
