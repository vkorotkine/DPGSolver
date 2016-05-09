// Modifiable Parameters
lc = 1; // Not used.

rIn = 1;
rOut = 1.384;



// Geometry Specification

Point(1) = {rIn,0,0,lc};
Point(2) = {rOut,0,0,lc};
Point(3) = {0,rIn,0,lc};
Point(4) = {rIn,rIn,0,lc};
Point(5) = {0,rOut,0,lc};
Point(6) = {rOut,rOut,0,lc};

Line(1001) = {1,2};
Line(1002) = {3,5};
Line(1003) = {4,3};
Line(1004) = {4,1};
Line(1005) = {6,5};
Line(1006) = {6,2};
Line(1007) = {4,6};

Transfinite Line {1003:1006}      = 4 Using Progression 1;
Transfinite Line {1001,1002,1007} = 2 Using Progression 1;

//Line Loop (4001) = {1003,1002,-1005,-1007};
Line Loop (4001) = {1007,1005,-1002,-1003};
Line Loop (4002) = {-1007,1004,1001,-1006};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface{4001} Left;
Transfinite Surface{4002} Right;



// Physical Parameters for '.msh' file

Physical Line(10001) = {1001,1002}; // Riemann Invariant Inflow/Outflow
Physical Line(20002) = {1003:1006}; // Slip-wall (curved)

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;