Include "/home/tmougin/xlifepp-sources-v2.3-2022-04-22/etc/gmsh/xlifepp_macros.geo";

h0=0.10000000000000001;
Call xlifepp_init;
l=0;
x1=0; y1=0; z1=0;
x2=10; y2=0; z2=0;
x3=10; y3=1; z3=0;
x4=0; y4=1; z4=0;
h1=h0; h2=h0; h3=h0; h4=h0;

Call xlifepp_Quadrangle;

Transfinite Line {L_1} = 60;
Transfinite Line {L_2} = 30;
Transfinite Line {L_3} = 60;
Transfinite Line {L_4} = 30;
domain_1={L_1[],L_3[]};
domain_2={L_4[]};
domain_3={L_2[]};


Call xlifepp_init;
l=1;
x1=0.5; y1=0.5; z1=0;
x2=0.59999999999999998; y2=0.5; z2=0;
x3=0.5; y3=0.59999999999999998; z3=0;
x4=0.40000000000000002; y4=0.5; z4=0;
x5=0.5; y5=0.40000000000000002; z5=0;
h1=h0; h2=h0; h3=h0; h4=h0; h5=h0;

Call xlifepp_Ellipse;

Transfinite Line {E_1} = 30;
Transfinite Line {E_2} = 30;
Transfinite Line {E_3} = 30;
Transfinite Line {E_4} = 30;

Plane Surface(loops_0)={loops_0[],loops_1[]};

domain_4={loops_0[]};

Physical Line("Gamma_a")= domain_1[];
Physical Line("Sigma_0")= domain_2[];
Physical Line("Sigma_a")= domain_3[];
Physical Surface("Omega")= domain_4[];


Mesh.ElementOrder=1;
Mesh.MshFileVersion = 2.2;
