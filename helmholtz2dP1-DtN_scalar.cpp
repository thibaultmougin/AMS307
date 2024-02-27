#include "xlife++.h"
using namespace xlifepp;

//inputs[1].PointData["real_part_term_u"]-inputs[0].PointData["real_part_term_u"]

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



Matrix<Complex> alpM(const Point& P, Parameters& pa=defaultParameters){
  Real x = P(1), a = pa("a");
  Complex alpha = pa("alpha");

  Vector<Complex> diago(2,1.);

  if (x<a){

      Matrix<Complex> M(diago);

      return M;
  }
  else{
    diago(1)=alpha;
    diago(2)=1./alpha;
    Matrix<Complex> M(diago);
    return M;
  }
}

Complex alpS(const Point& P, Parameters& pa=defaultParameters){
  Real x = P(1), a = pa("a");
  Complex alpha = pa("alpha");

  if (x<a){

      return 1;
  }
  else{
    return 1./alpha;
  }
}



int main(int argc, char** argv)
{
  init(argc, argv, _lang=en); // mandatory initialization of xlifepp

  Options opts ; 
  opts.add ( "a", 2. );

  opts.add ( "b", 3. );
  opts.add ( "k", 10.) ; opts.add ( "ny", 30);
  opts.add("N",5);
  opts.add("mod_alpha",0.25);
  opts.parse ( "data.txt" ) ;

  Real a =opts("a"),r=0.1, b =opts("b");

  Parameters params ;
  Real h=1. , k=opts("k");

  Complex alpha = opts("mod_alpha")*(1-i_);

  params << Parameter (h , "h")<< Parameter (k , "k" ) << Parameter(a,"a") << Parameter(b,"b") << Parameter (alpha , "alpha" );

  Number ny = opts("ny");
  Number na=Number(ny*a/h),nd=30;

  Function alphaM(alpM, "alphaM", params);
  Function alphaS(alpS,"alphaS",params);

  Rectangle Ra(_xmin =0, _xmax = b, _ymin = 0, _ymax = h, 
               _nnodes=Numbers(na,ny), 
               _domain_name = "Omega", _side_names=Strings ( "Gamma_a" , "Sigma_a" , "Gamma_a" , "Sigma_0" ));
  
  Disk D(_center=Point(0.5, 0.5), _radius =r, _nnodes=nd);
  Mesh mail(Ra, _triangle, 1, _gmsh);
  Mesh mesh2d(Ra-D, _triangle, 1, _gmsh);
  Domain omega=mesh2d.domain("Omega");
  Domain sigmaP=mesh2d.domain("Sigma_a"), sigmaM=mesh2d.domain("Sigma_0");
  
  cpuTime("mesh");

  std::cout<< std::endl ; 
  
  Space V(_domain=omega, _interpolation=P1, _name="V");
  Unknown u(V, "u"); TestFunction v(u, "v");

  Number N=opts("N");

  std::cout << std::endl << "k = " << k << ", N = " << N << ", a = " << a << ", b = " << b << ", ny = " << ny << ", mod_alpha = " << opts("mod_alpha") << std::endl;
  Space Sp(_domain=sigmaP, _basis=Function(cosny, params), _dim=N, _name="cos(n*pi*y)");


  
  Unknown phiP(Sp, "phiP");
  Vector<Complex> lambda(N);

  for (Number n=0; n<N; n++) lambda[n]=sqrt(Complex(k*k-n*n*pi_*pi_/(h*h)));
  /*
  for (Number n=0; n<N; n++) {
    Complex betan = sqrt(Complex(k*k-n*n*pi_*pi_/(h*h)));
    lambda[n]=betan*(exp(i_*betan*(b-a)));
    }
    */

  TensorKernel tkp(phiP, lambda);

  //BilinearForm auv = intg(omega, grad(u)|grad(v)) - k*k*intg(omega, u*v) - i_*k*intg(sigmaP, sigmaP, u*v);

  //BilinearForm auv = intg(omega, grad(u)|grad(v)) - k*k*intg(omega, u*v) - i_*intg(sigmaP, sigmaP, u*tkp*v);
  
  BilinearForm auv = intg(omega, (alphaM*grad(u))|grad(v)) - k*k*intg(omega, alphaS*u*v);

  LinearForm fv=intg(sigmaM, Function(gp, params)*v);
  TermMatrix A(auv, "A");
  TermVector B(fv, "B");
  
  cpuTime("assemblage");
  std::cout<< std::endl ; 

  TermVector U = directSolve(A, B);

  cpuTime("solver");
  std::cout<< std::endl ; 

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
