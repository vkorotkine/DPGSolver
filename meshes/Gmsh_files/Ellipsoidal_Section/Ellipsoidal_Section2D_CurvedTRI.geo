// Modifiable Parameters

Refine = 0;

lc = 0.6/2.0^Refine;

EqIndex = 1; // Options: 0 (Poisson), 1 (Euler)

rIn  = 0.5;
rOut = 1.0;

aIn  = 0.5;
aOut = 1.0;

bIn  = 0.5;
bOut = 1.0;

t   = Pi/4.0;
r   = rIn;
r_e = Sqrt(1.0/((Cos(t)/aOut)^2.0+(Sin(t)/bOut)^2.0));


// Geometry Specification

Point(1) = {aIn,0,0,lc};
Point(2) = {aOut,0,0,lc};
Point(3) = {0,bIn,0,lc};
Point(4) = {0,bOut,0,lc};
Point(5) = {r*Cos(t),r*Sin(t),0,lc};
Point(6) = {r_e*Cos(t),r_e*Sin(t),0,lc};
Point(7) = {0,0,0,lc};

Line(1001)    = {1,2};
Line(1002)    = {3,4};
Line(1003)    = {5,6};
Ellipse(1004) = {5,7,3,1}; // start, centre, point on major-axis, end
Ellipse(1005) = {5,7,3,3};
Ellipse(1006) = {6,7,4,2};
Ellipse(1007) = {6,7,4,4};


Transfinite Line {1001:1007} = 1*2^(Refine)+1 Using Progression 1;

Line Loop (4001) = {1001,-1006,-1003,1004};
Line Loop (4002) = {-1005,1003,1007,-1002};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface{4001} Right;
Transfinite Surface{4002};

//Recombine Surface{4002};




// Physical Parameters for '.msh' file

If (EqIndex == 0) // Poisson
	Physical Line(10011) = {1002}; // Straight Dirichlet
	Physical Line(10012) = {1001}; // Straight Neumann
	Physical Line(20011) = {1004:1005}; // Curved Dirichlet
	Physical Line(20012) = {1006:1007}; // Curved Neumann
ElseIf (EqIndex == 1) // Euler
	Physical Line(10001) = {1001,1002}; // Straight Riemann Invariant
	Physical Line(20002) = {1004:1007}; // Curved SlipWall
EndIf

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
