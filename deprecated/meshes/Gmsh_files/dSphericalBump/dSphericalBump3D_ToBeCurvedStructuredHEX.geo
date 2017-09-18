// Modifiable Parameters
lc = 1; // Not used.

rFactor = 4;
r = 1/rFactor;
H = rFactor*r;
L = rFactor*r;
W = rFactor*r;



// Geometry Specification

Point(1)  = {-L,0,-W,lc};
Point(2)  = { 0,0,-W,lc};
Point(3)  = { L,0,-W,lc};
Point(4)  = {-L,H,-W,lc};
Point(5)  = { 0,H,-W,lc};
Point(6)  = { L,H,-W,lc};
Point(7)  = {-r,0,-r,lc};
Point(8)  = { 0,0,-r,lc};
Point(9)  = { r,0,-r,lc};
Point(10) = {-r,r,-r,lc};
Point(11) = { 0,r,-r,lc};
Point(12) = { r,r,-r,lc};
Point(13) = {-L,0, 0,lc};
Point(14) = {-r,0, 0,lc};
Point(15) = { r,0, 0,lc};
Point(16) = { L,0, 0,lc};
Point(17) = {-r,r, 0,lc};
Point(18) = { 0,r, 0,lc};
Point(19) = { r,r, 0,lc};
Point(20) = {-L,H, 0,lc};
Point(21) = { 0,H, 0,lc};
Point(22) = { L,H, 0,lc};

Line(1001) = {1,2};
Line(1002) = {2,3};
Line(1003) = {4,5};
Line(1004) = {5,6};
Line(1005) = {1,4};
Line(1006) = {2,5};
Line(1007) = {3,6};

Line(1008) = {7,8};
Line(1009) = {8,9};
Line(1010) = {10,11};
Line(1011) = {11,12};
Line(1012) = {7,10};
Line(1013) = {8,11};
Line(1014) = {9,12};

Line(1015) = {13,14};
Line(1016) = {15,16};
Line(1017) = {17,18};
Line(1018) = {18,19};
Line(1019) = {20,21};
Line(1020) = {21,22};
Line(1021) = {13,20};
Line(1022) = {14,17};
Line(1023) = {15,19};
Line(1024) = {16,22};
Line(1025) = {18,21};

Line(1026) = {1,13};
Line(1027) = {1,7};
Line(1028) = {2,8};
Line(1029) = {3,9};
Line(1030) = {3,16};
Line(1031) = {7,14};
Line(1032) = {9,15};

Line(1033) = {10,17};
Line(1034) = {11,18};
Line(1035) = {12,19};

Line(1036) = {4,20};
Line(1037) = {4,10};
Line(1038) = {5,21};
Line(1039) = {5,11};
Line(1040) = {6,22};
Line(1041) = {6,12};

Line(1042) = {20,17};
Line(1043) = {22,19};

Transfinite Line {1001:1014,1017:1024,1026,1030:1036,1038,1040} = 2 Using Progression 1;
Transfinite Line {-1015,1016,1025,-1027,-1028,-1029,
                  -1037,-1039,-1041,-1042,-1043} = rFactor Using Progression 1;

Line Loop (4001) = {1001, 1006,-1003,-1005};
Line Loop (4002) = {1002, 1007,-1004,-1006};
Line Loop (4003) = {1008, 1013,-1010,-1012};
Line Loop (4004) = {1009, 1014,-1011,-1013};
Line Loop (4005) = {1015, 1022,-1042,-1021};
Line Loop (4006) = {1017, 1025,-1019, 1042};
Line Loop (4007) = {1018,-1043,-1020,-1025};
Line Loop (4008) = {1016, 1024, 1043,-1023};
Line Loop (4009) = {1026, 1021,-1036,-1005};
Line Loop (4010) = {1031, 1022,-1033,-1012};
Line Loop (4011) = {1028, 1013,-1039,-1006};
Line Loop (4012) = {1039, 1034, 1025,-1038};
Line Loop (4013) = {1032, 1023,-1035,-1014};
Line Loop (4014) = {1030, 1024,-1040,-1007};
Line Loop (4015) = {1027, 1031,-1015,-1026};
Line Loop (4016) = {1027, 1008,-1028,-1001};
Line Loop (4017) = {1028, 1009,-1029,-1002};
Line Loop (4018) = {1029, 1032, 1016,-1030};
Line Loop (4019) = {1003, 1038,-1019,-1036};
Line Loop (4020) = {1004, 1040,-1020,-1038};
Line Loop (4021) = {1005, 1037,-1012,-1027};
Line Loop (4022) = {1007, 1041,-1014,-1029};
Line Loop (4023) = {1037, 1033,-1042,-1036};
Line Loop (4024) = {1039,-1010,-1037, 1003};
Line Loop (4025) = {1039, 1011,-1041,-1004};
Line Loop (4026) = {1041, 1035,-1043,-1040};
Line Loop (4027) = {1010, 1034,-1017,-1033};
Line Loop (4028) = {1011, 1035,-1018,-1034};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};
Plane Surface(4004) = {4004};
Plane Surface(4005) = {4005};
Plane Surface(4006) = {4006};
Plane Surface(4007) = {4007};
Plane Surface(4008) = {4008};
Plane Surface(4009) = {4009};
Plane Surface(4010) = {4010};
Plane Surface(4011) = {4011};
Plane Surface(4012) = {4012};
Plane Surface(4013) = {4013};
Plane Surface(4014) = {4014};
Plane Surface(4015) = {4015};
Plane Surface(4016) = {4016};
Plane Surface(4017) = {4017};
Plane Surface(4018) = {4018};
Plane Surface(4019) = {4019};
Plane Surface(4020) = {4020};
Plane Surface(4021) = {4021};
Plane Surface(4022) = {4022};
Plane Surface(4023) = {4023};
Plane Surface(4024) = {4024};
Plane Surface(4025) = {4025};
Plane Surface(4026) = {4026};
Plane Surface(4027) = {4027};
Plane Surface(4028) = {4028};

Transfinite Surface {4001:4028};
Recombine Surface{4001:4028};

Surface Loop (7001) = {4005,4009,4010,4015,4021,4023};
Surface Loop (7002) = {4001,4003,4011,4016,4021,4024};
Surface Loop (7003) = {4006,4012,4019,4023,4024,4027};
Surface Loop (7004) = {4002,4004,4011,4017,4022,4025};
Surface Loop (7005) = {4008,4013,4014,4018,4022,4026};
Surface Loop (7006) = {4007,4012,4020,4025,4026,4028};

Volume(7001) = {7001};
Volume(7002) = {7002};
Volume(7003) = {7003};
Volume(7004) = {7004};
Volume(7005) = {7005};
Volume(7006) = {7006};

Transfinite Volume{7001:7006};
Recombine Volume{7001:7006};



// Physical parameters for '.msh' file

Physical Surface(10001) = {4001,4002,4009,4014,4019,4020}; // Riemann Invariant Inflow/Outflow
//Physical Surface(10003) = {4014};                        // Outflow Mach
Physical Surface(10002) = {4005:4008,4015:4018};           // Slip-wall
Physical Surface(20002) = {4003,4004,4010,4013,4027,4028}; // Slip-wall (curved)

Physical Volume(9701) = {7001:7006};



// Visualization in gmsh

Color Black{ Volume{7001:7006}; }
Color Black{ Surface{4001:4028}; }
Geometry.Color.Points = Black;




