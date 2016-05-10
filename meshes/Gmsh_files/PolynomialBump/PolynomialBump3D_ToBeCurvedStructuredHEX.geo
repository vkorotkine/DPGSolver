// Modifiable Parameters
lc = 1; // Not used.

l = 2; L = 4;
H = 4;
w = 2; W = 4;

Point(1)  = {-L,+0,-0,lc};
Point(2)  = {-l,+0,-0,lc};
Point(3)  = {+l,+0,-0,lc};
Point(4)  = {+L,+0,-0,lc};
Point(5)  = {-L,+H,-0,lc};
Point(6)  = {-l,+H,-0,lc};
Point(7)  = {+l,+H,-0,lc};
Point(8)  = {+L,+H,-0,lc};

Point(9)  = {-L,+0,-w,lc};
Point(10) = {-l,+0,-w,lc};
Point(11) = {+l,+0,-w,lc};
Point(12) = {+L,+0,-w,lc};
Point(13) = {-L,+H,-w,lc};
Point(14) = {-l,+H,-w,lc};
Point(15) = {+l,+H,-w,lc};
Point(16) = {+L,+H,-w,lc};

Point(17) = {-L,+0,-W,lc};
Point(18) = {-l,+0,-W,lc};
Point(19) = {+l,+0,-W,lc};
Point(20) = {+L,+0,-W,lc};
Point(21) = {-L,+H,-W,lc};
Point(22) = {-l,+H,-W,lc};
Point(23) = {+l,+H,-W,lc};
Point(24) = {+L,+H,-W,lc};

Line(1001) = {1,2};
Line(1002) = {2,3};
Line(1003) = {3,4};
Line(1004) = {5,6};
Line(1005) = {6,7};
Line(1006) = {7,8};
Line(1007) = {1,5};
Line(1008) = {2,6};
Line(1009) = {3,7};
Line(1010) = {4,8};

Line(1011) = {1+8,2+8};
Line(1012) = {2+8,3+8};
Line(1013) = {3+8,4+8};
Line(1014) = {5+8,6+8};
Line(1015) = {6+8,7+8};
Line(1016) = {7+8,8+8};
Line(1017) = {1+8,5+8};
Line(1018) = {2+8,6+8};
Line(1019) = {3+8,7+8};
Line(1020) = {4+8,8+8};

Line(1021) = {1+16,2+16};
Line(1022) = {2+16,3+16};
Line(1023) = {3+16,4+16};
Line(1024) = {5+16,6+16};
Line(1025) = {6+16,7+16};
Line(1026) = {7+16,8+16};
Line(1027) = {1+16,5+16};
Line(1028) = {2+16,6+16};
Line(1029) = {3+16,7+16};
Line(1030) = {4+16,8+16};

Line(1031) = {1,9};
Line(1032) = {2,10};
Line(1033) = {3,11};
Line(1034) = {4,12};
Line(1035) = {5,13};
Line(1036) = {6,14};
Line(1037) = {7,15};
Line(1038) = {8,16};

Line(1039) = {1+8,9+8};
Line(1040) = {2+8,10+8};
Line(1041) = {3+8,11+8};
Line(1042) = {4+8,12+8};
Line(1043) = {5+8,13+8};
Line(1044) = {6+8,14+8};
Line(1045) = {7+8,15+8};
Line(1046) = {8+8,16+8};

Transfinite Line{1001,1003,1004,1006,1011,1013,1014,1016,1021,1023,1024,1026,1031:1046} = 2 Using Progression 1;
Transfinite Line{1002,1005,1012,1015,1022,1025,1007:1010,1017:1020,1027:1030} = 3 Using Progression 1;

Line Loop (4001) = {1001,1008,-1004,-1007};
Line Loop (4002) = {1002,1009,-1005,-1008};
Line Loop (4003) = {1003,1010,-1006,-1009};
Line Loop (4004) = {1011,1018,-1014,-1017};
Line Loop (4005) = {1012,1019,-1015,-1018};
Line Loop (4006) = {1013,1020,-1016,-1019};
Line Loop (4007) = {1021,1028,-1024,-1027};
Line Loop (4008) = {1022,1029,-1025,-1028};
Line Loop (4009) = {1023,1030,-1026,-1029};

Line Loop (4010) = {1001,1032,-1011,-1031};
Line Loop (4011) = {1002,1033,-1012,-1032};
Line Loop (4012) = {1003,1034,-1013,-1033};
Line Loop (4013) = {1011,1040,-1021,-1039};
Line Loop (4014) = {1012,1041,-1022,-1040};
Line Loop (4015) = {1013,1042,-1023,-1041};

Line Loop (4016) = {1004,1036,-1014,-1035};
Line Loop (4017) = {1005,1037,-1015,-1036};
Line Loop (4018) = {1006,1038,-1016,-1037};
Line Loop (4019) = {1014,1044,-1024,-1043};
Line Loop (4020) = {1015,1045,-1025,-1044};
Line Loop (4021) = {1016,1046,-1026,-1045};

Line Loop (4022) = {1031,1017,-1035,-1007};
Line Loop (4023) = {1039,1027,-1043,-1017};
Line Loop (4024) = {1032,1018,-1036,-1008};
Line Loop (4025) = {1040,1028,-1044,-1018};
Line Loop (4026) = {1033,1019,-1037,-1009};
Line Loop (4027) = {1041,1029,-1045,-1019};
Line Loop (4028) = {1034,1020,-1038,-1010};
Line Loop (4029) = {1042,1030,-1046,-1020};


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
Plane Surface(4029) = {4029};


Transfinite Surface {4001:4029};
Recombine Surface{4001:4029};


Surface Loop (7001) = {4001,4004,4010,4016,4022,4024};
Surface Loop (7002) = {4002,4005,4011,4017,4024,4026};
Surface Loop (7003) = {4003,4006,4012,4018,4026,4028};
Surface Loop (7004) = {4004,4007,4013,4019,4023,4025};
Surface Loop (7005) = {4005,4008,4014,4020,4025,4027};
Surface Loop (7006) = {4006,4009,4015,4021,4027,4029};

Volume(7001) = {7001};
Volume(7002) = {7002};
Volume(7003) = {7003};
Volume(7004) = {7004};
Volume(7005) = {7005};
Volume(7006) = {7006};

Transfinite Volume{7001:7006};
Recombine Volume{7001:7006};



// Physical parameters for '.msh' file

Physical Surface(10001) = {4022,4023,4028,4029}; // Riemann Invariant Inflow/Outflow
Physical Surface(10002) = {4001,4003,4007:4021}; // Slip-wall
Physical Surface(20002) = {4002,4011};           // Slip-wall (curved)

Physical Volume(9701) = {7001:7006};



// Visualization in gmsh

Color Black{ Volume{7001:7006}; }
Color Black{ Surface{4001:4029}; }
Geometry.Color.Points = Black;




