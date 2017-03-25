Include "../Parameters.geo";
MeshLevel = 0; MeshType = WEDGE; PDEName = EULER;

// Geometry Specification
rIn  = 1.000;
rOut = 1.384;
dz   = 0.384;

Point(1) = {rIn,0,0,lc};
Point(2) = {rOut,0,0,lc};
Point(3) = {0,rIn,0,lc};
Point(4) = {rIn,rIn,0,lc};
Point(5) = {0,rOut,0,lc};
Point(6) = {rOut,rOut,0,lc};

Point(7) = {rIn,0,dz,lc};
Point(8) = {rOut,0,dz,lc};
Point(9) = {0,rIn,dz,lc};
Point(10) = {rIn,rIn,dz,lc};
Point(11) = {0,rOut,dz,lc};
Point(12) = {rOut,rOut,dz,lc};

Line(1001) = {1,2};
Line(1002) = {3,5};
Line(1003) = {4,3};
Line(1004) = {4,1};
Line(1005) = {6,5};
Line(1006) = {6,2};
Line(1007) = {4,6};

Line(1008) = {7,8};
Line(1009) = {9,11};
Line(1010) = {10,9};
Line(1011) = {10,7};
Line(1012) = {12,11};
Line(1013) = {12,8};
Line(1014) = {10,12};

Line(1015) = {1,7};
Line(1016) = {2,8};
Line(1017) = {3,9};
Line(1018) = {4,10};
Line(1019) = {5,11};
Line(1020) = {6,12};

Transfinite Line {1003:1006,1010:1013}           = 2^(MeshLevel)+1 Using Progression 1;
Transfinite Line {1001,1002,1007,1008,1009,1014} = 2^(MeshLevel)+1 Using Progression 1;
Transfinite Line {1015:1020}                     = 2^(MeshLevel)+1 Using Progression 1;

Line Loop (4001) = {1002,1019,-1009,-1017};
Line Loop (4002) = {1007,1020,-1014,-1018};
Line Loop (4003) = {1001,1016,-1008,-1015};

Line Loop (4004) = {1003,1017,-1010,-1018};
Line Loop (4005) = {1004,1015,-1011,-1018};
Line Loop (4006) = {1005,1019,-1012,-1020};
Line Loop (4007) = {1006,1016,-1013,-1020};

Line Loop (4008) = {1007,1005,-1002,-1003};
//Line Loop (4008) = {1003,1002,-1005,-1007};
Line Loop (4009) = {-1007,1004,1001,-1006};
Line Loop (4010) = {1014,1012,-1009,-1010};
//Line Loop (4010) = {1010,1009,-1012,-1014};
Line Loop (4011) = {-1014,1011,1008,-1013};

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

If (MeshType == TET)
	// Do nothing
ElseIf (MeshType == HEX)
	Transfinite Surface{4001:4011};
	Recombine Surface{4001:4011};
ElseIf (MeshType == WEDGE)
	Transfinite Surface{4001:4007,4009,4011};
	Transfinite Surface{4008,4010} Right;
	Recombine Surface{4001:4007};
ElseIf (MeshType == MIXED3D_TP)
	Recombine Surface{4001:4011};
ElseIf (MeshType == MIXED3D_HW)
	Transfinite Surface{4001:4011};
//	Recombine Surface{4001:4007,4009,4011};
//	Recombine Surface{4001:4003,4005,4007:4011};
	Recombine Surface{4001:4008,4010};
EndIf

Surface Loop (7001) = {4001,4002,4004,4006,4008,4010};
Surface Loop (7002) = {4002,4003,4005,4007,4009,4011};

Volume(7001) = {7001};
Volume(7002) = {7002};

If (MeshType == TET)
    // Do nothing
ElseIf (MeshType == HEX)
	Transfinite Volume{7001:7002};
//	Recombine Volume{7001,7002};
ElseIf (MeshType == WEDGE)
	Transfinite Volume{7001:7002};
	Recombine Volume{7001,7002};
ElseIf (MeshType == MIXED3D_TP)
	Recombine Volume{7001};
ElseIf (MeshType == MIXED3D_HW)
	Transfinite Volume{7001:7002};
	Recombine Volume{7001,7002};
EndIf



// Physical parameters for '.msh' file
BC_Straight =   BC_STEP_SC;
BC_Curved   = 2*BC_STEP_SC;

If (PDEName == EULER)
	Physical Line    (BC_Curved  +BC_SLIPWALL) = {1003:1006,1010:1013};
	Physical Surface (BC_Straight+BC_RIEMANN)  = {4001,4003};
	Physical Surface (BC_Curved  +BC_SLIPWALL) = {4004:4007};
	Physical Surface (BC_Straight+BC_SLIPWALL) = {4008:4011};
EndIf

Physical Volume(9701) = {7001:7002};



// Visualization in gmsh

Color Black{ Volume{7001:7002}; }
Color Black{ Surface{4001:4011}; }
Geometry.Color.Points = Black;
