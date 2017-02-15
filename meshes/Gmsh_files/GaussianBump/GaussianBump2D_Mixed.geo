// Modifiable Parameters

Refine = 0;

lc = 0.5/2.0^Refine;

// RefType Options: 0 (Std. transfinite), 1 (Std.)
RefType = 0;

BumpExp = 2;
N = 1e2;



// Geometry Specification
a = 0.0625;
b = 0.0;
c = 0.2/2^BumpExp;

l = 0.5;
h = 4.0*a;


xL = -l; x = xL; yL = a*Exp(-(x-b)^2/(2*c^2));
xR =  l; x = xR; yR = a*Exp(-(x-b)^2/(2*c^2));
xM =  0; x = xM; yM = a*Exp(-(x-b)^2/(2*c^2));


Point(1) = {xL,yL,0,lc};
Point(2) = {xR,yR,0,lc};
Point(3) = {xL,h,0,lc};
Point(4) = {xR,h,0,lc};
Point(5) = {xM,yM,0,lc};
Point(6) = {xM,h,0,lc};

pListBL[0] = 1;
For i In {1:N-1}
	x = -l + i/N*l;
	y = a*Exp(-(x-b)^2/(2*c^2));
	
	pListBL[i] = newp;
	Point(pListBL[i]) = {x,y,0,lc};
EndFor
pListBL[N] = 5;

pListBR[0] = 5;
For i In {1:N-1}
	x = 0 + i/N*l;
	y = a*Exp(-(x-b)^2/(2*c^2));
	
	pListBR[i] = newp;
	Point(pListBR[i]) = {x,y,0,lc};
EndFor
pListBR[N] = 2;

Spline(1001) = pListBL[];
Spline(1002) = pListBR[];
Line(1003) = {3,6};
Line(1004) = {6,4};
Line(1005) = {1,3};
Line(1006) = {2,4};
Line(1007) = {5,6};


If (RefType == 0)
	Transfinite Line {1001:1004} = 1*2^(Refine)+1 Using Progression 1;
	Transfinite Line {1005:1007} = 1*2^(Refine)+1 Using Progression 1;
ElseIf (RefType == 1)
	Transfinite Line {1001:1004} = 4*2^(Refine)   Using Progression 1;
	Transfinite Line {1005:1007} = 1*2^(Refine)+1 Using Progression 1;
EndIf

Line Loop (4001) = {1001,1007,-1003,-1005};
Line Loop (4002) = {1002,1006,-1004,-1007};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

If (RefType == 0 || RefType == 1)
	Transfinite Surface{4001};
	Transfinite Surface{4002} Right;
EndIf

Recombine Surface{4001};



// Physical Parameters for '.msh' file

Physical Line(10011) = {1003,1004,1001}; // Straight Dirichlet
Physical Line(10012) = {1005,1006,1002}; // Straight Neumann
//Physical Line(20011) = {1001};      // Curved Dirichlet
//Physical Line(20012) = {1002};      // Curved Neumann

Physical Surface(9401) = {4001:4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
Coherence;
