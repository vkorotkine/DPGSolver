// Modifiable Parameters

Refine = 0;

MeshType = 0; // Options: 0 (Transfinite), 1 (Refined TE)

lc = 5/2.0^Refine;


// Geometry Specification
a = 1.0;
l = a/2.25;



g  = Pi; 
xL = a*2*((Cos(g)^2 - 2*Cos(g) + 1)*l^3 + (2*Cos(g)^2 - 3*Cos(g) + 1)*l^2 + (Cos(g)^2 - 2*Cos(g) + 1)*l - Cos(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);

g  = 0; 
xR = a*2*((Cos(g)^2 - 2*Cos(g) + 1)*l^3 + (2*Cos(g)^2 - 3*Cos(g) + 1)*l^2 + (Cos(g)^2 - 2*Cos(g) + 1)*l - Cos(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);

g  = 0.5*Pi; 
xM = a*2*((Cos(g)^2 - 2*Cos(g) + 1)*l^3 + (2*Cos(g)^2 - 3*Cos(g) + 1)*l^2 + (Cos(g)^2 - 2*Cos(g) + 1)*l - Cos(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);
yM = a*2*(l^3*(Cos(g) - 1)*Sin(g) + 2*l^2*(Cos(g) - 1)*Sin(g) + l*(Cos(g) - 1)*Sin(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);


lL = -xL+2*a;
lU = -xL+1*a;
h  = yM+2*a;


N = 1e1;


Point(0) = {-lL,0,0,lc};
Point(1) = {xL,0,0,lc};
Point(2) = {xR,0,0,lc/20};
Point(3) = {lL,0,0,lc};

Point(4) = {-lL,h,0,lc};
Point(5) = {-lU,h,0,lc};
Point(6) = {lU,h,0,lc};
Point(7) = {lL,h,0,lc};

Point(8) = {xM,yM,0,lc};
Point(9) = {xM,h,0,lc};


pListBL[0] = 1;
For i In {1:N-1}
	g = Pi-Pi/2*(i/N);
	x = a*2*((Cos(g)^2 - 2*Cos(g) + 1)*l^3 + (2*Cos(g)^2 - 3*Cos(g) + 1)*l^2 + (Cos(g)^2 - 2*Cos(g) + 1)*l - Cos(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);
	y = a*2*(l^3*(Cos(g) - 1)*Sin(g) + 2*l^2*(Cos(g) - 1)*Sin(g) + l*(Cos(g) - 1)*Sin(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);
	
	pListBL[i] = newp;
	Point(pListBL[i]) = {x,y,0,lc};
EndFor
pListBL[N] = 8;


pListBR[0] = 8;
For i In {1:N-1}
	g = Pi/2-Pi/2*(i/N);
	x = a*2*((Cos(g)^2 - 2*Cos(g) + 1)*l^3 + (2*Cos(g)^2 - 3*Cos(g) + 1)*l^2 + (Cos(g)^2 - 2*Cos(g) + 1)*l - Cos(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);
	y = a*2*(l^3*(Cos(g) - 1)*Sin(g) + 2*l^2*(Cos(g) - 1)*Sin(g) + l*(Cos(g) - 1)*Sin(g))/(2*l^2*(Cos(g) - 1) + 2*l*(Cos(g) - 1) - 1);
	
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


Line Loop (4001) = {1001,1010,-1005,-1009};
Line Loop (4002) = {1002,1013,-1006,-1010};
Line Loop (4003) = {1003,1011,-1007,-1013};
Line Loop (4004) = {1004,1012,-1008,-1011};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};
Plane Surface(4004) = {4004};



If (MeshType == 0)
	Transfinite Line {1001:1008} = 1*2^(Refine)+1 Using Progression 1;
	Transfinite Line {1009:1013} = 1*2^(Refine)+1 Using Progression 1;

	Transfinite Surface{4001,4003} Right;
	Transfinite Surface{4002,4004};
ElseIf (MeshType == 1)
	Transfinite Line {1001:1002,1005,1006} = 1*2^(Refine)+1 Using Progression 1;
	Transfinite Line {1009:1010,1013} = 1*2^(Refine)+1 Using Progression 1;

	Transfinite Surface{4001} Right;
	Transfinite Surface{4002} Left;
EndIf


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
