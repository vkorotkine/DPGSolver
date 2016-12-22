// Modifiable Parameters

Refine = 0;

lc = 1.0/2.0^Refine;

h = 0.5;
l = 1.0;
h15 = 0.133974596215561353236276829247;

// Geometry Specification

Point(1) = {-l,0,0,lc};
Point(2) = {0,0,0,lc};
Point(3) = {0.5,h15,0,lc};
Point(4) = {-l,h,0,lc};
Point(5) = {0,h,0,lc};
Point(6) = {0.5,h,0,lc};

Line(1001) = {1,2};
Line(1002) = {2,3};
Line(1003) = {4,5};
Line(1004) = {5,6};
Line(1005) = {1,4};
Line(1006) = {2,5};
Line(1007) = {3,6};

Transfinite Line {1001,1003} = 2*2^(Refine)+1 Using Progression 1;
Transfinite Line {1002,1004} = 1*2^(Refine)+1 Using Progression 1;
Transfinite Line {1005:1007} = 1*2^(Refine)+1 Using Progression 1;

Line Loop (4001) = {1001,1006,-1003,-1005};
Line Loop (4002) = {1002,1007,-1004,-1006};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface{4001};
Transfinite Surface{4002};

//Recombine Surface{4001};
//Recombine Surface{4002};



// Physical Parameters for '.msh' file

Physical Line(10011) = {1003,1004,1001}; // Straight Dirichlet
Physical Line(10012) = {1005,1007,1002}; // Straight Neumann
//Physical Line(20011) = {1003:1004}; // Curved Dirichlet
//Physical Line(20012) = {1005:1006}; // Curved Neumann

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;