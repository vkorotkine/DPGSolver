// Modifiable Parameters

Refine = 0;

lc = 1.0/2.0^Refine;

L = 1.0;



// Geometry Specification

Point(1) = {0,0,0,lc/3};
Point(2) = {L,0,0,lc/2};
Point(3) = {0,L,0,lc};
Point(4) = {L,L,0,lc};

Line(1001) = {1,2};
Line(1002) = {3,4};
Line(1003) = {1,3};
Line(1004) = {2,4};


Line Loop (4001) = {1001,1004,-1002,-1003};

Plane Surface(4001) = {4001};

//Recombine Surface{4001};



// Physical Parameters for '.msh' file

Physical Line(10011) = {1001:1004}; // Straight Dirichlet
//Physical Line(10011) = {1002,1004}; // Straight Dirichlet
//Physical Line(10012) = {1001,1003}; // Straight Neumann

Physical Surface(9401) = {4001};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
