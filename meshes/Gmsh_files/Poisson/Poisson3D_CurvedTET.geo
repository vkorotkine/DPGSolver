// Modifiable Parameters

Refine = 0;

lc = 0.6/2.0^Refine;

rIn = 0.5;
rOut = 1.0;



// Geometry Specification

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {rIn,0.0,0.0,lc};
Point(3) = {0.0,rIn,0.0,lc};
Point(4) = {0.0,0.0,rIn,lc};
Point(5) = {rOut,0.0,0.0,lc};
Point(6) = {0.0,rOut,0.0,lc};
Point(7) = {0.0,0.0,rOut,lc};

Circle(1001) = {2,1,3};
Circle(1002) = {4,1,3};
Circle(1003) = {2,1,4};
Circle(1004) = {5,1,6};
Circle(1005) = {7,1,6};
Circle(1006) = {5,1,7};
Line(1007) = {2,5};
Line(1008) = {3,6};
Line(1009) = {4,7};


Line Loop(4001) = {1003,1002,-1001};
Line Loop(4002) = {1006,1005,-1004};
Line Loop(4003) = {1003,1009,-1006,-1007};
Line Loop(4004) = {-1002,1009,1005,-1008};
Line Loop(4005) = {-1001,1007,1004,-1008};

Ruled Surface(4001) = {4001};
Ruled Surface(4002) = {4002};
Plane Surface(4003) = {4003};
Plane Surface(4004) = {4004};
Plane Surface(4005) = {4005};

Surface Loop (7001) = {4001:4005};
Volume(7001) = {7001};



// Physical parameters for '.msh' file

Physical Surface(10011) = {4003:4004}; // Straight Dirichlet
Physical Surface(10012) = {4005};      // Straight Neumann
Physical Surface(20011) = {4001};      // Curved Dirichlet
Physical Surface(20012) = {4002};      // Curved Neumann

Physical Volume(9701) = 7001;


/*
Line(1001) = {1,2};
Line(1002) = {3,5};
Circle(1003) = {4,7,3};
Circle(1004) = {4,7,1};
Circle(1005) = {6,7,5};
Circle(1006) = {6,7,2};
Line(1007) = {4,6};




Line Loop (4001) = {1007,1005,-1002,-1003};
Line Loop (4002) = {-1007,1004,1001,-1006};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

//Recombine Surface{4002};



// Physical Parameters for '.msh' file

Physical Line(10011) = {1002}; // Straight Dirichlet
Physical Line(10012) = {1001}; // Straight Neumann
Physical Line(20011) = {1003:1004}; // Curved Dirichlet
Physical Line(20012) = {1005:1006}; // Curved Neumann

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};
*/


// Visualization in gmsh

Color Black{ Surface{4001:4005}; }
Color Black{ Volume{7001}; }
Geometry.Color.Points = Black;





/*
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {1,0.0,0.0,lc};
Point(3) = {0,1,0.0,lc};
Point(4) = {0,0,1,lc};

Circle(1001) = {2,1,3};
Circle(1002) = {4,1,3};
Circle(1003) = {2,1,4};

Line(1004) = {1,2};
Line(1005) = {1,3};
Line(1006) = {1,4};

Line Loop(4008) = {1004,1003,-1006};
Line Loop(4010) = {1005,-1002,-1006};
Line Loop(4012) = {1004,1001,-1005};
Line Loop(4017) = {-1002,-1003,1001};

Plane Surface(4009) = {4008};
Plane Surface(4011) = {4010};
Plane Surface(4013) = {4012};
Ruled Surface(4018) = {4017};

Surface Loop(7029) = {4009,4011,4013,4018};
Volume(7030) = {7029};

Physical Surface(10011) = {4009,4011};
Physical Surface(10012) = {4013};
Physical Surface(20012) = {4018};

Physical Volume(9701) = {7030};
*/
