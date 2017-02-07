// Modifiable Parameters

Refine = 0;

lc = 0.5/2.0^Refine;


// Geometry Specification
c  = 1.0;
lL = 1.0;
lR = 1.0;
lO = 0.4;
h  = 1.2;
r  = 0.25;
t  = c*Sqrt(r)/1.1019;

a0 =  0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 =  0.2843;
a4 = -0.1036;

N = 1e2;


Point(0) = {-lL,0,0,lc};
Point(1) = {0,0,0,lc};
Point(2) = {c,0,0,lc};
Point(3) = {c+lR,0,0,lc};
Point(4) = {-lL,h,0,lc};
Point(5) = {0-lO,h,0,lc};
Point(6) = {c+lO,h,0,lc};
Point(7) = {c+lR,h,0,lc};

x = 0.5;
y = 5*c*t*(a0*Sqrt(x/c)+a1*(x/c)^1+a2*(x/c)^2+a3*(x/c)^3+a4*(x/c)^4);

Point(8) = {x,y,0,lc};
Point(9) = {x,h,0,lc};


pListBL[0] = 1;
For i In {1:N-1}
	x = c/2*(i/N);
	y = 5*c*t*(a0*Sqrt(x/c)+a1*(x/c)^1+a2*(x/c)^2+a3*(x/c)^3+a4*(x/c)^4);
	
	pListBL[i] = newp;
	Point(pListBL[i]) = {x,y,0,lc};
EndFor
pListBL[N] = 8;

pListBR[0] = 8;
For i In {1:N-1}
	x = c/2*(1+i/N);
	y = 5*c*t*(a0*Sqrt(x/c)+a1*(x/c)^1+a2*(x/c)^2+a3*(x/c)^3+a4*(x/c)^4);
	
	pListBR[i] = newp;
	Point(pListBR[i]) = {x,y,0,lc};
EndFor
pListBR[N] = 2;


Line(1001)   = {0,1};
Spline(1002) = pListBL[];
Spline(1003) = pListBR[];
Line(1004)   = {2,3};
Line(1005)   = {4,5};
Line(1006)   = {5,9};
Line(1007)   = {9,6};

Line(1008)   = {6,7};
Line(1009)   = {0,4};
Line(1010)   = {1,5};
Line(1011)   = {2,6};
Line(1012)   = {3,7};
Line(1013)   = {8,9};


Transfinite Line {1001:1008} = 1*2^(Refine)+1 Using Progression 1;
Transfinite Line {1009:1013} = 1*2^(Refine)+1 Using Progression 1;

Line Loop (4001) = {1001,1010,-1005,-1009};
Line Loop (4002) = {1002,1013,-1006,-1010};
Line Loop (4003) = {1003,1011,-1007,-1013};
Line Loop (4004) = {1004,1012,-1008,-1011};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};
Plane Surface(4004) = {4004};

Transfinite Surface{4001:4004};
//Recombine Surface{4001:4004};




// Physical Parameters for '.msh' file

Physical Line(10001) = {1009,1012};           // Straight Riemann
Physical Line(10002) = {1001,1004,1005:1008}; // Straight SlipWall
Physical Line(20002) = {1002,1003};           // Curved SlipWall

Physical Surface(9401) = {4001:4004};



// Visualization in gmsh

Color Black{ Surface{4001:4004}; }
Geometry.Color.Points = Black;
Coherence;
