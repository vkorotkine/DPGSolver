// Modifiable Parameters
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

Transfinite Line{1001:1004} = 2 Using Progression 1;
Transfinite Line{2001:2004} = 2 Using Progression 1;
Transfinite Line{3001:3004} = 2 Using Progression 1;

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

Physical Surface(10001) = {6001:6002,5001:5002,4001:4002};

Physical Volume(9701) = 7001;



// Visualization in gmsh

Color Black{ Surface{4001:4002,5001:5002,6001:6002}; }
Color Black{ Volume{7001}; }
Geometry.Color.Points = Black;