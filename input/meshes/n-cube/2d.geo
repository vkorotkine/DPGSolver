Include "../parameters.geo";
//mesh_level = 2; mesh_type = TRI; mesh_domain = STRAIGHT; pde_name = ADVECTION;

// Geometry Specification
l = 1;
h = 1;

Point(1) = {-l,-h,-0,lc};
Point(2) = {+l,-h,-0,lc};
Point(3) = {-l,+h,-0,lc};
Point(4) = {+l,+h,-0,lc};

Line(1001) = {1,2};
Line(1002) = {3,4};
Line(2001) = {1,3};
Line(2002) = {2,4};

Transfinite Line{1001:1002} = 2^(mesh_level)+1 Using Progression 1;
Transfinite Line{2001:2002} = 2^(mesh_level)+1 Using Progression 1;

Line Loop (4001) = {1001,2002,-1002,-2001};

Plane Surface(4001) = {4001};

Transfinite Surface {4001};
If (mesh_type == TRI)
	// Do nothing.
ElseIf (mesh_type == QUAD)
	Recombine Surface{4001};
Else
	Error("Unsupported mesh_type: %d",mesh_type); Exit;
EndIf



// Physical parameters for '.msh' file
bc_straight =   BC_STEP_SC;
bc_curved   = 2*BC_STEP_SC;
If (mesh_domain == STRAIGHT)
	bc_base = bc_straight;
Else
	bc_base = bc_curved;
EndIf

If (pde_name == ADVECTION)
	If (geom_adv == GEOM_ADV_YL)
		Physical Line(bc_base+BC_INFLOW)  = {1001};
		Physical Line(bc_base+BC_OUTFLOW) = {2001,1002,2002};
	ElseIf (geom_adv == GEOM_ADV_XYL)
		Physical Line(bc_base+BC_INFLOW)  = {1001,2001};
		Physical Line(bc_base+BC_OUTFLOW) = {1002,2002};
	Else
		Error("Unsupported geom_adv: %d",geom_adv); Exit;
	EndIf
ElseIf (pde_name == POISSON)
	Physical Line(bc_base+BC_DIRICHLET) = {1001:1002,2001:2002};
ElseIf (pde_name == EULER)
	Physical Line(bc_base+PERIODIC_XL) = {2001};
	Physical Line(bc_base+PERIODIC_XR) = {2002};
	Physical Line(bc_base+PERIODIC_YL) = {1001};
	Physical Line(bc_base+PERIODIC_YR) = {1002};

	// Periodic Indicator (Slave = Master)
	Periodic Line {2002} = {2001}; // Periodic (x)
	Periodic Line {1002} = {1001}; // Periodic (y)
ElseIf (pde_name == NAVIERSTOKES)
	Physical Line(bc_base+PERIODIC_XL) = {2001};
	Physical Line(bc_base+PERIODIC_XR) = {2002};

	Periodic Line {2002} = {2001}; // Periodic (x)

	// No slip for remaining boundaries
	Physical Line (bc_base+BC_NOSLIP_T)         = {1001};
	Physical Line (bc_base+BC_NOSLIP_ADIABATIC) = {1002};
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

Physical Surface(9401) = 4001;



// Visualization in gmsh

Color Black{ Surface{4001}; }
Geometry.Color.Points = Black;
