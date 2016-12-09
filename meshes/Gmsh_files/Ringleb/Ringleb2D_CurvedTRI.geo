// Modifiable Parameters

Refine = 0;

lc = 1.0/2.0^Refine;

NnIn = 1e2;


// Geometry Specification

Point(1) = { 0.328869654659162, 2.273162711063530,0.0,lc};
Point(2) = {-1.485914825006039, 1.429069867868291,0.0,lc};
Point(3) = { 0.328869654659162,-2.273162711063530,0.0,lc};
Point(4) = {-1.485914825006039,-1.429069867868291,0.0,lc};
Point(5) = { 1.384672639155810, 0.0              ,0.0,lc};
Point(6) = {-0.060963132143099, 0.0              ,0.0,lc};


GAMMA = 1.4;
GM1 = GAMMA-1;

Q0   = 0.5;
KMin = 0.7;
KMax = 1.5;


qMin = Q0;
qMax = Q0;
kMin = KMin;
kMax = KMax;


For i In {1:NnIn-1}
	q = qMin;
	k = kMin+i/NnIn*(kMax-kMin);

    a = Sqrt(1.0-GM1/2.0*q*q);
    rho = a^(2.0/GM1);
    J = 1.0/a+1.0/(3.0*a^3)+1.0/(5.0*a^5.0)-0.5*Log((1.0+a)/(1.0-a));

    x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-J/2.0;
    y = 1.0/(k*rho*q)*Sqrt(1.0-(q/k)^2);

	pListBIf[i] = newp;
	Point(pListBIf[i]) = {x,y,0,lc};

    pListBOf[i] = newp;
    Point(pListBOf[i]) = {x,-y,0,lc};
EndFor



For i In {1:NnIn-1}
    q = Q0-i/NnIn*(Q0-kMax);
    k = kMax;

    a = Sqrt(1.0-GM1/2.0*q*q);
    rho = a^(2.0/GM1);
    J = 1.0/a+1.0/(3.0*a^3)+1.0/(5.0*a^5.0)-0.5*Log((1.0+a)/(1.0-a));

    x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-J/2.0;
    y = 1.0/(k*rho*q)*Sqrt(1.0-(q/k)^2);

    pListBIw1[i] = newp;
    Point(pListBIw1[i]) = {x,y,0,lc};


    q = Q0-i/NnIn*(Q0-kMin);
    k = kMin;

    a = Sqrt(1.0-GM1/2.0*q*q);
    rho = a^(2.0/GM1);
    J = 1.0/a+1.0/(3.0*a^3)+1.0/(5.0*a^5.0)-0.5*Log((1.0+a)/(1.0-a));

    x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-J/2.0;
    y = 1.0/(k*rho*q)*Sqrt(1.0-(q/k)^2);

    pListBOw1[i] = newp;
    Point(pListBOw1[i]) = {x,y,0,lc};
EndFor

For i In {1:NnIn-1}
    q = Q0+i/NnIn*(kMax-Q0);
    k = kMax;

    a = Sqrt(1.0-GM1/2.0*q*q);
    rho = a^(2.0/GM1);
    J = 1.0/a+1.0/(3.0*a^3)+1.0/(5.0*a^5.0)-0.5*Log((1.0+a)/(1.0-a));

    x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-J/2.0;
    y = -1.0/(k*rho*q)*Sqrt(1.0-(q/k)^2);

    pListBIw2[i] = newp;
    Point(pListBIw2[i]) = {x,y,0,lc};


    q = Q0+i/NnIn*(kMin-Q0);
    k = kMin;

    a = Sqrt(1.0-GM1/2.0*q*q);
    rho = a^(2.0/GM1);
    J = 1.0/a+1.0/(3.0*a^3)+1.0/(5.0*a^5.0)-0.5*Log((1.0+a)/(1.0-a));

    x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-J/2.0;
    y = -1.0/(k*rho*q)*Sqrt(1.0-(q/k)^2);

    pListBOw2[i] = newp;
    Point(pListBOw2[i]) = {x,y,0,lc};
EndFor

pListBIf[0] = 1; pListBIf[NnIn] = 2;
pListBOf[0] = 3; pListBOf[NnIn] = 4;

pListBIw1[0] = 2; pListBIw1[NnIn] = 6;
pListBOw1[0] = 1; pListBOw1[NnIn] = 5;

pListBIw2[0] = 4; pListBIw2[NnIn] = 6;
pListBOw2[0] = 3; pListBOw2[NnIn] = 5;


Spline(1001) = pListBIf[];
Spline(1002) = pListBOf[];
Spline(1003) = pListBIw1[];
Spline(1004) = pListBOw1[];
Spline(1005) = pListBIw2[];
Spline(1006) = pListBOw2[];
Line(1007)   = {5,6};

Line Loop (4001) = {1001,1003,-1007,-1004};
Line Loop (4002) = {1002,1005,-1007,-1006};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

//Recombine Surface{4002};



// Physical Parameters for '.msh' file

Physical Line(20011) = {1001:1002}; // Curved Dirichlet
Physical Line(20012) = {1003:1006}; // Curved Neumann

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;