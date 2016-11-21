// Modifiable Parameters

Refine = 0;

lc = 0.6/2.0^Refine;

rIn = 0.5;
rOut = 1.0;



// Geometry Specification

Point(1) = {rIn,0,0,lc/2};
Point(2) = {rOut,0,0,lc};
Point(3) = {0,rIn,0,lc/2};
Point(4) = {Sqrt(0.5)*rIn,Sqrt(0.5)*rIn,0,lc};
Point(5) = {0,rOut,0,lc};
Point(6) = {Sqrt(0.5)*rOut,Sqrt(0.5)*rOut,0,lc};
Point(7) = {0,0,0,lc};

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



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
