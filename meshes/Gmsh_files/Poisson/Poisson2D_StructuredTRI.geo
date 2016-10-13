// Modifiable Parameters
lc = 1; // Not used.

L = 1;
H = 1;



// Geometry Specification

Point(1) = {-L,-H,-0,lc};
Point(2) = {+L,-H,-0,lc};
Point(3) = {-L,+H,-0,lc};
Point(4) = {+L,+H,-0,lc};

Line(1001) = {1,2};
Line(1002) = {3,4};
Line(2001) = {1,3};
Line(2002) = {2,4};

Transfinite Line{1001:1002} = 2 Using Progression 1;
Transfinite Line{2001:2002} = 2 Using Progression 1;

Line Loop (4001) = {1001,2002,-1002,-2001};

Plane Surface(4001) = {4001};

Transfinite Surface {4001};
//Recombine Surface{4001};



// Physical parameters for '.msh' file

Physical Line(10011) = {2001,1001}; // Dirichlet
Physical Line(10012) = {2002,1002}; // Neumann

Physical Surface(9401) = 4001;




// Visualization in gmsh

Color Black{ Surface{4001}; }
Geometry.Color.Points = Black;
