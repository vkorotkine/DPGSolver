// Modifiable Parameters
lc = 1; // Not used.

l = 2; L = 4;
H = 4;

Point(1) = {-L,+0,-0,lc};
Point(2) = {-l,+0,-0,lc};
Point(3) = {+l,+0,-0,lc};
Point(4) = {+L,+0,-0,lc};
Point(5) = {-L,+H,-0,lc};
Point(6) = {-l,+H,-0,lc};
Point(7) = {+l,+H,-0,lc};
Point(8) = {+L,+H,-0,lc};

Line(1001) = {1,2};
Line(1002) = {2,3};
Line(1003) = {3,4};
Line(1004) = {5,6};
Line(1005) = {6,7};
Line(1006) = {7,8};
Line(1007) = {1,5};
Line(1008) = {2,6};
Line(1009) = {3,7};
Line(1010) = {4,8};

Transfinite Line{1001,1003,1004,1006} = 2 Using Progression 1;
Transfinite Line{1002,1005,1007:1010} = 3 Using Progression 1;

Line Loop (4001) = {1001,1008,-1004,-1007};
Line Loop (4002) = {1002,1009,-1005,-1008};
Line Loop (4003) = {1003,1010,-1006,-1009};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};

Transfinite Surface {4001:4003};
Recombine Surface{4001};
Recombine Surface{4002};
Recombine Surface{4003};



// Physical parameters for '.msh' file

Physical Line(10001) = {1007,1010};      // Riemann Invariant Inflow/Outflow
Physical Line(10002) = {1001,1003:1006}; // Slip-wall
Physical Line(20002) = {1002};           // Slip-wall (curved)

Physical Surface(9401) = {4001:4003};



// Visualization in gmsh

Color Black{ Surface{4001:4003}; }
Geometry.Color.Points = Black;


