Include "../parameters.geo";
//MESH_LEVEL = 2; MESH_TYPE = TRI; MESH_DOMAIN = STRAIGHT; PDE_NAME = ADVECTION;

// Geometry Specification
L = 1;
H = 1;

Point(1) = {-L,-H,-0,lc};
Point(2) = {+L,-H,-0,lc};
Point(3) = {-L,+H,-0,lc};
Point(4) = {+L,+H,-0,lc};

Line(1001) = {1,2};
Line(1002) = {3,4};
Line(2001) = {1,3};
Line(2002) = {2,4};

Transfinite Line{1001:1002} = 2^(MESH_LEVEL)+1 Using Progression 1;
Transfinite Line{2001:2002} = 2^(MESH_LEVEL)+1 Using Progression 1;

Line Loop (4001) = {1001,2002,-1002,-2001};

Plane Surface(4001) = {4001};

Transfinite Surface {4001};
If (MESH_TYPE == QUAD)
	Recombine Surface{4001};
EndIf



// Physical parameters for '.msh' file
BC_Straight =   BC_STEP_SC;
BC_Curved   = 2*BC_STEP_SC;
If (MESH_DOMAIN == STRAIGHT)
	BC_Base = BC_Straight;
Else
	BC_Base = BC_Curved;
EndIf

If (PDE_NAME == ADVECTION)
	If (Geom_Adv == GEOM_ADV_YL)
		Physical Line(BC_Base+BC_INFLOW)  = {1001};
		Physical Line(BC_Base+BC_OUTFLOW) = {2001,1002,2002};
	EndIf
ElseIf (PDE_NAME == POISSON)
	Physical Line(BC_Base+BC_DIRICHLET) = {1001:1002,2001:2002};
ElseIf (PDE_NAME == EULER)
	Physical Line(BC_Base+PERIODIC_XL) = {2001};
	Physical Line(BC_Base+PERIODIC_XR) = {2002};
	Physical Line(BC_Base+PERIODIC_YL) = {1001};
	Physical Line(BC_Base+PERIODIC_YR) = {1002};

	// Periodic Indicator (Slave = Master)
	Periodic Line {2002} = {2001}; // Periodic (x)
	Periodic Line {1002} = {1001}; // Periodic (y)
ElseIf (PDE_NAME == NAVIERSTOKES)
	Physical Line(BC_Base+PERIODIC_XL) = {2001};
	Physical Line(BC_Base+PERIODIC_XR) = {2002};

	Periodic Line {2002} = {2001}; // Periodic (x)

	// No slip for remaining boundaries
	Physical Line (BC_Base+BC_NOSLIP_T)         = {1001};
	Physical Line (BC_Base+BC_NOSLIP_ADIABATIC) = {1002};
EndIf

Physical Surface(9401) = 4001;



// Visualization in gmsh

Color Black{ Surface{4001}; }
Geometry.Color.Points = Black;
