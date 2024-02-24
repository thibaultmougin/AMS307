#include "xlife++.h"
using namespace xlifepp;

// data on sigma-
Complex gp(const Point& P, Parameters& pa = defaultParameters)
{
  Real y=P(2);
  Number n =0;
  Real h=pa ( "h" ) ;
  Real k=pa ( "k" ) ;
  Complex betan = sqrt(Complex(k*k-n*n*pi_*pi_/(h*h)));
  //return -i_*betan*sqrt(2./h)*cos(n*pi_*y/h);
  return -i_*k;
}

// spectral functions
Real cosny(const Point& P, Parameters& pa = defaultParameters)
{
  Real y=P(2), h=pa("h");            // get the parameter h (user definition)
  Number n=pa("basis index")-1;         // get the index of function to compute
  if (n==0) return sqrt(1./h);
  else      return sqrt(2./h)*cos(n*pi_*y/h);
}

Complex alp(const Point& P, Parameters& pa=defaultParameters){
  Real x = P(1), a = pa("a");
  Complex alp = 1-1*i_;
  if (x<a){
      return 1.;
  }
  else{
    return alp;
  }
}

int main(int argc, char** argv)
{
  init(argc, argv, _lang=en); // mandatory initialization of xlifepp

  Real a =2,r=0.1, b =6.;

  Parameters params ;
  Real h=1. , k=5.;


  params << Parameter (h , "h")<< Parameter (k , "k" ) << Parameter(a,"a") << Parameter(b,"b");

  Number ny=30, na=Number(ny*a/h),nd=30;

  Function alpha(alp, "alpha", params);

  Rectangle Ra(_xmin =0, _xmax = b, _ymin = 0, _ymax = h, 
               _nnodes=Numbers(na,ny), 
               _domain_name = "Omega", _side_names=Strings ( "Gamma_a" , "Sigma_a" , "Gamma_a" , "Sigma_0" ));
  
  Disk D(_center=Point(0.5, 0.5), _radius =r, _nnodes=nd);
  Mesh mail(Ra, _triangle, 1, _gmsh);
  Mesh mesh2d(Ra-D, _triangle, 1, _gmsh);
  Domain omega=mesh2d.domain("Omega");
  Domain sigmaP=mesh2d.domain("Sigma_a"), sigmaM=mesh2d.domain("Sigma_0");

  Space V(_domain=omega, _interpolation=P1, _name="V");
  Unknown u(V, "u"); TestFunction v(u, "v");

  Number N=5;
  Space Sp(_domain=sigmaP, _basis=Function(cosny, params), _dim=N, _name="cos(n*pi*y)");
  Unknown phiP(Sp, "phiP");
  Vector<Complex> lambda(N);
  for (Number n=0; n<N; n++) lambda[n]=sqrt(Complex(k*k-n*n*pi_*pi_/(h*h)));
  TensorKernel tkp(phiP, lambda);

  //BilinearForm auv = intg(omega, grad(u)|grad(v)) - k*k*intg(omega, u*v) - i_*intg(sigmaP, sigmaP, u*tkp*v);
  
  BilinearForm auv = intg(omega, grad(u)|grad(v)) - k*k*intg(omega, u*v) - i_*intg(sigmaP, sigmaP, u*tkp*v);

  //BilinearForm auv = intg(omega, alpha*grad(u)|grad(v)) - k*k*intg(omega, u*v);

  LinearForm fv=intg(sigmaM, Function(gp, params)*v);
  TermMatrix A(auv, "A");
  TermVector B(fv, "B");

  TermVector U = directSolve(A, B);
  saveToFile("U", U, vtu);

/*
  TermMatrix M( intg (Omega, u*v ) ) ;
  TermVector Uex(u , Omega, uex ) ;
  TermVector E=U-Uex ;
  Real er=sqrt ( abs ((M*E|E ) ) ) ;
  std::cout<<"L2 error = "<<er<<eol ;
  */

  return 0;
}
