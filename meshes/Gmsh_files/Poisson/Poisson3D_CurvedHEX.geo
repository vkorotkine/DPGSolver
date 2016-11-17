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
Point(5) = {1/Sqrt(2.0)*rIn,1/Sqrt(2.0)*rIn,0.0,lc};
Point(6) = {1/Sqrt(2.0)*rIn,0.0,1/Sqrt(2.0)*rIn,lc};
Point(7) = {0.0,1/Sqrt(2.0)*rIn,1/Sqrt(2.0)*rIn,lc};
Point(8) = {1/Sqrt(3.0)*rIn,1/Sqrt(3.0)*rIn,1/Sqrt(3.0)*rIn,lc};

Point(12) = {rOut,0.0,0.0,lc/2};
Point(13) = {0.0,rOut,0.0,lc/2};
Point(14) = {0.0,0.0,rOut,lc/2};
Point(15) = {1/Sqrt(2.0)*rOut,1/Sqrt(2.0)*rOut,0.0,lc/2};
Point(16) = {1/Sqrt(2.0)*rOut,0.0,1/Sqrt(2.0)*rOut,lc/2};
Point(17) = {0.0,1/Sqrt(2.0)*rOut,1/Sqrt(2.0)*rOut,lc/2};
Point(18) = {1/Sqrt(3.0)*rOut,1/Sqrt(3.0)*rOut,1/Sqrt(3.0)*rOut,lc/2};


Circle(1001) = {2,1,5};
Circle(1002) = {2,1,6};
Circle(1003) = {3,1,5};
Circle(1004) = {3,1,7};
Circle(1005) = {4,1,6};
Circle(1006) = {4,1,7};
Circle(1007) = {5,1,8};
Circle(1008) = {6,1,8};
Circle(1009) = {7,1,8};

Circle(1011) = {12,1,15};
Circle(1012) = {12,1,16};
Circle(1013) = {13,1,15};
Circle(1014) = {13,1,17};
Circle(1015) = {14,1,16};
Circle(1016) = {14,1,17};
Circle(1017) = {15,1,18};
Circle(1018) = {16,1,18};
Circle(1019) = {17,1,18};

Line(1022) = {2,12};
Line(1023) = {3,13};
Line(1024) = {4,14};
Line(1025) = {5,15};
Line(1026) = {6,16};
Line(1027) = {7,17};
Line(1028) = {8,18};


Line Loop(4001) = {1001,1007,-1008,-1002}; Ruled Surface(4001) = {4001};
Line Loop(4002) = {1003,1007,-1009,-1004}; Ruled Surface(4002) = {4002};
Line Loop(4003) = {1005,1008,-1009,-1006}; Ruled Surface(4003) = {4003};
Line Loop(4004) = {1011,1017,-1018,-1012}; Ruled Surface(4004) = {4004};
Line Loop(4005) = {1013,1017,-1019,-1014}; Ruled Surface(4005) = {4005};
Line Loop(4006) = {1015,1018,-1019,-1016}; Ruled Surface(4006) = {4006};

Line Loop(4007) = {1001,1025,-1011,-1022}; Plane Surface(4007) = {4007};
Line Loop(4008) = {1002,1026,-1012,-1022}; Plane Surface(4008) = {4008};
Line Loop(4009) = {1003,1025,-1013,-1023}; Plane Surface(4009) = {4009};
Line Loop(4010) = {1004,1027,-1014,-1023}; Plane Surface(4010) = {4010};
Line Loop(4011) = {1005,1026,-1015,-1024}; Plane Surface(4011) = {4011};
Line Loop(4012) = {1006,1027,-1016,-1024}; Plane Surface(4012) = {4012};

Line Loop(4013) = {1007,1028,-1017,-1025}; Plane Surface(4013) = {4013};
Line Loop(4014) = {1008,1028,-1018,-1026}; Plane Surface(4014) = {4014};
Line Loop(4015) = {1009,1028,-1019,-1027}; Plane Surface(4015) = {4015};


Surface Loop (7001) = {4001,4004,4007,4008,4013,4014}; Volume(7001) = {7001};
Surface Loop (7002) = {4002,4005,4009,4010,4013,4015}; Volume(7002) = {7002};
Surface Loop (7003) = {4003,4006,4011,4012,4014,4015}; Volume(7003) = {7003};


Transfinite Line {1001:1009,1011:1019} = 2*2^Refine+1 Using Progression 1;
Transfinite Line {1022:1028}           = 3*2^Refine+1 Using Progression 1;

Transfinite Surface "*";
Transfinite Volume "*";

// HEX
Recombine Surface "*";

// WEDGE
//Recombine Surface {4007:4015};                                         // TRI
//Recombine Surface {4001:4006,4008,4009,4010,4011,4013,4015};           // QUAD
//Recombine Surface {4001,4003,4004,4006,4008,4009,4010,4011,4013,4015}; // Mixed

// PYR
// Comment all “Transfinite” options and uncomment the recombine option below.
//Recombine Surface {4004:4012};

// TET
// No Recombine



// Physical parameters for '.msh' file

Physical Surface(10011) = {4007:4010}; // Straight Dirichlet
Physical Surface(10012) = {4011:4012}; // Straight Neumann
Physical Surface(20011) = {4001:4003}; // Curved Dirichlet
Physical Surface(20012) = {4004:4006}; // Curved Neumann

Physical Volume(9701) = {7001:7003};



// Visualization in gmsh

Color Black{ Surface{4001:4014}; }
Color Black{ Volume{7001:7003}; }
Geometry.Color.Points = Black;