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

Physical Point(20051) = {1,3}; // Periodic (-x)
Physical Point(20052) = {2,4}; // Periodic (+x)
Physical Point(20053) = {1,2}; // Periodic (-y)
Physical Point(20054) = {3,4}; // Periodic (+y)

Physical Line(20051) = {2001}; // Periodic (-x)
Physical Line(20052) = {2002}; // Periodic (+x)
Physical Line(20053) = {1001}; // Periodic (-y)
Physical Line(20054) = {1002}; // Periodic (+y)

Physical Surface(9401) = 4001;

// Periodic Indicator (Slave = Master)

Periodic Line {2002} = {2001}; // Periodic (x)
Periodic Line {1002} = {1001}; // Periodic (y)



// Visualization in gmsh

Color Black{ Surface{4001}; }
Geometry.Color.Points = Black;
