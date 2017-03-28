Include "../Parameters.geo";
//MeshCurving = CURVED; MeshType = TRI; PDEName = NAVIERSTOKES; MeshLevel = 0;


// Geometry Specification
rIn  = 0.5;
rOut = 1.0;

If (MeshCurving == TOBECURVED)
	IndP = 1; r = rIn;
	Point(IndP++) = {-r,-r,0,lc};
	Point(IndP++) = { r,-r,0,lc};
	Point(IndP++) = {-r, r,0,lc};
	Point(IndP++) = { r, r,0,lc}; r = rOut;
	Point(IndP++) = {-r,-r,0,lc};
	Point(IndP++) = { r,-r,0,lc};
	Point(IndP++) = {-r, r,0,lc};
	Point(IndP++) = { r, r,0,lc};


	IndL = 1001; IndP = 0;
	Line(IndL++) = {1+IndP,2+IndP};
	Line(IndL++) = {3+IndP,4+IndP};
	Line(IndL++) = {1+IndP,3+IndP};
	Line(IndL++) = {2+IndP,4+IndP}; IndP = 4;
	Line(IndL++) = {1+IndP,2+IndP};
	Line(IndL++) = {3+IndP,4+IndP};
	Line(IndL++) = {1+IndP,3+IndP};
	Line(IndL++) = {2+IndP,4+IndP}; IndP = 1;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
ElseIf (MeshCurving == CURVED)
	IndP = 0;
	Point(IndP++) = {0,0,0,lc};
	r = rIn;
	t = 5.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};
	t = 7.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};
	t = 3.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};
	t = 1.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};
	r = rOut;
	t = 5.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};
	t = 7.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};
	t = 3.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};
	t = 1.0/4.0*Pi; Point(IndP++) = {r*Cos(t),r*Sin(t),0,lc};


	IndL = 1001; IndP = 0;
	Circle(IndL++) = {1+IndP,0,2+IndP};
	Circle(IndL++) = {3+IndP,0,4+IndP};
	Circle(IndL++) = {1+IndP,0,3+IndP};
	Circle(IndL++) = {2+IndP,0,4+IndP}; IndP = 4;
	Circle(IndL++) = {1+IndP,0,2+IndP};
	Circle(IndL++) = {3+IndP,0,4+IndP};
	Circle(IndL++) = {1+IndP,0,3+IndP};
	Circle(IndL++) = {2+IndP,0,4+IndP}; IndP = 1;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
	Line(IndL++) = {IndP,IndP+4};   IndP++;
EndIf


// Include something for aspect ratio: 1.0, 2.5, 5.0, 20.0
Transfinite Line {1001:1008} = 2^(MeshLevel+1)+1 Using Progression 1;
Transfinite Line {1009:1012} = 2^(MeshLevel)+1 Using Progression 1;

IndL = 4001;
//Line Loop (IndL++) = {1001,1010,-1005,-1009};
Line Loop (IndL++) = {1005,-1010,-1001,1009};
Line Loop (IndL++) = {1002,1012,-1006,-1011};
Line Loop (IndL++) = {1003,1011,-1007,-1009};
//Line Loop (IndL++) = {1004,1012,-1008,-1010};
Line Loop (IndL++) = {1008,-1012,-1004,1010};

IndL = 4001;
For (1:4)
//	Plane Surface(IndL) = {IndL}; IndL++; // Makes a more equilateral mesh
	Plane Surface(IndL) = {IndL}; Transfinite Surface{IndL}; IndL++;
EndFor


If (MeshType == MIXED2D)
	Recombine Surface{4001,4003};
ElseIf (MeshType == QUAD)
	Recombine Surface{4001:4004};
EndIf



// Physical parameters for '.msh' file
BC_Curved   = 2*BC_STEP_SC;

If (PDEName == NAVIERSTOKES)
	// Use: BC_EXACT and BC_NOSLIP_ADIABATIC
	Physical Line (BC_Curved+BC_DIRICHLET)        = {1001:1004};
	Physical Line (BC_Curved+BC_NOSLIP_ADIABATIC) = {1005:1008};
EndIf

IndV = 9401; IndS = 4001;
For (1:4)
	Physical Surface(IndV++) = {IndS++};
EndFor



// Visualization in gmsh

Color Black{ Surface{4001:4004}; }
Geometry.Color.Points = Black;
