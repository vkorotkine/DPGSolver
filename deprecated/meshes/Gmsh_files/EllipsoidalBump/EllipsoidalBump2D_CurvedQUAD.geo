// Modifiable Parameters

Refine = 0;
BoundaryOut = 1; // Options: 0 (Riemann), 1 (BackPressure)

lc = 0.5/2.0^Refine;


// Geometry Specification
a = 0.5;
b = a/3;

lL = 3*a;
lR = lL;
h  = 3*b;
lO = Tan(Pi/4)*h;

Point(0)  = {-lL-a,0,0,lc};
Point(1)  = {-a,0,0,lc};
Point(2) = {0,b,0,lc};
Point(3)  = {a,0,0,lc};
Point(4)  = {lR+a,0,0,lc};
Point(5)  = {-lL-a,h,0,lc};
Point(6)  = {-a-lO,h,0,lc};
Point(7)  = {0,h,0,lc};
Point(8)  = {a+lO,h,0,lc};
Point(9)  = {lR+a,h,0,lc};
Point(10) = {0,0,0,lc};

Line(1001)    = {0,1};
Ellipse(1002) = {1,10,1,2}; // start, centre, point on major-axis, end
Ellipse(1003) = {2,10,3,3}; // start, centre, point on major-axis, end
Line(1004)    = {3,4};
Line(1005)    = {5,6};
Line(1006)    = {6,7};
Line(1007)    = {7,8};
Line(1008)    = {8,9};

Line(1009)    = {0,5};
Line(1010)    = {1,6};
Line(1011)    = {2,7};
Line(1012)    = {3,8};
Line(1013)    = {4,9};

Transfinite Line {1001,1004:1005,1008} = 1*2^(Refine+1)+1 Using Progression 1;
Transfinite Line {1002:1003,1006:1007} = 1*2^(Refine)+1 Using Progression 1;
Transfinite Line {1009:1013} = 1*2^(Refine)+1 Using Progression 1;

Line Loop (4001) = {1001,1010,-1005,-1009};
Line Loop (4002) = {1002,1011,-1006,-1010};
Line Loop (4003) = {1003,1012,-1007,-1011};
Line Loop (4004) = {1004,1013,-1008,-1012};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};
Plane Surface(4004) = {4004};

Transfinite Surface{4001,4003} Right;
Transfinite Surface{4002,4004};
Recombine Surface{4001:4004};




// Physical Parameters for '.msh' file


If (BoundaryOut == 0)
	Physical Line(10001) = {1009,1013}; // Straight Riemann
ElseIf (BoundaryOut == 1)
	Physical Line(10004) = {1009};      // Total Pressure/Temperature
	Physical Line(10003) = {1013};      // Straight BackPressure
EndIf

Physical Line(10002) = {1001,1004,1005:1008}; // Straight SlipWall
Physical Line(20002) = {1002,1003};           // Curved   SlipWall

Physical Surface(9401) = {4001:4004};



// Visualization in gmsh

Color Black{ Surface{4001:4004}; }
Geometry.Color.Points = Black;
Coherence;