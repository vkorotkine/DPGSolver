/*********************************************************************
 *
 *  Gmsh
 *
 *  3D Periodic Vortex (ToBeCurved)
 *
 *********************************************************************/



// Modifiable Parameters
Refine = 1;

lc = 1; // Not used.

L = 1;
H = 1;
W = 1;



// Geometry Specification

Point(1) = {-L,-H,-W,lc};
Point(2) = {+L,-H,-W,lc};
Point(3) = {-L,+H,-W,lc};
Point(4) = {+L,+H,-W,lc};

Point(5) = {-L,-H,+W,lc};
Point(6) = {+L,-H,+W,lc};
Point(7) = {-L,+H,+W,lc};
Point(8) = {+L,+H,+W,lc};

Line(1001) = {1,2};
Line(1002) = {3,4};
Line(1003) = {1+4,2+4};
Line(1004) = {3+4,4+4};

Line(2001) = {1,3};
Line(2002) = {2,4};
Line(2003) = {1+4,3+4};
Line(2004) = {2+4,4+4};

Line(3001) = {1,1+4};
Line(3002) = {2,2+4};
Line(3003) = {3,3+4};
Line(3004) = {4,4+4};

Transfinite Line{1001:1004} = 2^(Refine)+1 Using Progression 1;
Transfinite Line{2001:2004} = 2^(Refine)+1 Using Progression 1;
Transfinite Line{3001:3004} = 2^(Refine)+1 Using Progression 1;

Line Loop (4001) = {1001,2002,-1002,-2001};
Line Loop (4002) = {1003,2004,-1004,-2003};
Line Loop (5001) = {1001,3002,-1003,-3001};
Line Loop (5002) = {1002,3004,-1004,-3003};
Line Loop (6001) = {2001,3003,-2003,-3001};
Line Loop (6002) = {2002,3004,-2004,-3002};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(5001) = {5001};
Plane Surface(5002) = {5002};
Plane Surface(6001) = {6001};
Plane Surface(6002) = {6002};

Transfinite Surface {4001:4002,5001:5002,6001:6002};
Recombine   Surface {4001:4002};
Recombine   Surface {5001:5002};
Recombine   Surface {6001:6002};

Surface Loop (7001) = {4001:4002,5001:5002,6001:6002};

Volume(7001) = {7001};

Transfinite Volume{7001};
Recombine Volume{7001};



// Physical parameters for '.msh' file

Physical Point(20051) = {1,3,5,7}; // Periodic (-x)
Physical Point(20052) = {2,4,6,8}; // Periodic (+x)
Physical Point(20053) = {1,2,5,6}; // Periodic (-y)
Physical Point(20054) = {3,4,7,8}; // Periodic (+y)
Physical Point(20055) = {1,2,3,4}; // Periodic (-z)
Physical Point(20056) = {5,6,7,8}; // Periodic (+z)

Physical Line(20051) = {2001,2003,3001,3003}; // Periodic (-x)
Physical Line(20052) = {2002,2004,3002,3004}; // Periodic (+x)
Physical Line(20053) = {1001,1003,3001,3002}; // Periodic (-y)
Physical Line(20054) = {1002,1004,3003,3004}; // Periodic (+y)
Physical Line(20055) = {1001,1002,2001,2002}; // Periodic (-z)
Physical Line(20056) = {1003,1004,2003,2004}; // Periodic (+z)

Physical Surface(20051) = {6001}; // Periodic (-x)
Physical Surface(20052) = {6002}; // Periodic (+x)
Physical Surface(20053) = {5001}; // Periodic (-y)
Physical Surface(20054) = {5002}; // Periodic (+y)
Physical Surface(20055) = {4001}; // Periodic (-z)
Physical Surface(20056) = {4002}; // Periodic (+z)

Physical Volume(9701) = 7001;


// Periodic Indicator (Slave = Master)

Periodic Line {2002} = {2001}; // Periodic (x)
Periodic Line {2004} = {2003};
Periodic Line {3002} = {3001};
Periodic Line {3004} = {3003};

Periodic Line {1002} = {1001}; // Periodic (y)
Periodic Line {1004} = {1003};
Periodic Line {3003} = {3001};
Periodic Line {3004} = {3002};

Periodic Line {1003} = {1001}; // Periodic (z)
Periodic Line {1004} = {1002};
Periodic Line {2003} = {2001};
Periodic Line {2004} = {2002};

Periodic Surface 6002 {2002,3004,-2004,-3002} = 6001 {2001,3003,-2003,-3001}; // Periodic (x)
Periodic Surface 5002 {1002,3004,-1004,-3003} = 5001 {1001,3002,-1003,-3001}; // Periodic (y)
Periodic Surface 4002 {1003,2004,-1004,-2003} = 4001 {1001,2002,-1002,-2001}; // Periodic (z)



// Visualization in gmsh

Color Black{ Surface{4001:4002,5001:5002,6001:6002}; }
Color Black{ Volume{7001}; }
Geometry.Color.Points = Black;