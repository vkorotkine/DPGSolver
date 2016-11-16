// Modifiable Parameters

Refine = 0;

lc = 0.6/2.0^Refine;

rIn = 0.5;
rOut = 1.0;



// Geometry Specification

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {rIn,0.0,0.0,lc/2};
Point(3) = {0.0,rIn,0.0,lc/2};
Point(4) = {0.0,0.0,rIn,lc/2};
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

Transfinite Surface{4001:4005};
Recombine Surface{4001:4005};

Surface Loop (7001) = {4001:4005};


Volume(7001) = {7001};
//Transfinite Volume{7001};
Recombine Volume{7001};



// Physical parameters for '.msh' file

Physical Surface(10011) = {4003:4004}; // Straight Dirichlet
Physical Surface(10012) = {4005};      // Straight Neumann
Physical Surface(20011) = {4001};      // Curved Dirichlet
Physical Surface(20012) = {4002};      // Curved Neumann

Physical Volume(9701) = 7001;



// Visualization in gmsh

Color Black{ Surface{4001:4005}; }
Color Black{ Volume{7001}; }
Geometry.Color.Points = Black;
