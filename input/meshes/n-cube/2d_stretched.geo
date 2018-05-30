Include "../parameters.geo";
//mesh_level = 2; mesh_type = MIXED; mesh_domain = STRAIGHT; pde_name = ADVECTION; geom_adv = GEOM_ADV_XL; geom_unaligned = 1;
//mesh_level = 2; mesh_type = QUAD; mesh_domain = PARAMETRIC; pde_name = EULER; geom_adv = GEOM_ADV_INTERNAL; geom_unaligned = 0; geom_ar = 0.5;

// Geometry Specification
l = 1;
h = 1;

Point(1) = {-l,-h,-0,lc};
Point(2) = {+l,-h,-0,lc};
Point(3) = {-l,+h,-0,lc};
Point(4) = {+l,+h,-0,lc};
Point(5) = {-0,-h,-0,lc};
Point(6) = {+0,+h,-0,lc};

Line(1001) = {1,5};
Line(1002) = {5,2};
Line(1003) = {3,6};
Line(1004) = {6,4};
Line(2001) = {1,3};
Line(2002) = {2,4};
Line(2003) = {5,6};

If (!geom_unaligned)
	prog_spec = 1;
Else
	prog_spec = 1.1;
EndIf

aspect_ratio = geom_ar;

Printf("aspect_ratio ~= %g.",aspect_ratio);
If (aspect_ratio == GMSH_DUMMY)
	Transfinite Line{1001:1002} = 2^(mesh_level)*2+1   Using Progression 1;
	Transfinite Line{2001:2002} = 2^(mesh_level+1)+1 Using Progression 1.3;
	Transfinite Line{2003}      = 2^(mesh_level+1)+1 Using Progression 1.3;
	Transfinite Line{1003:1004} = 2^(mesh_level)*2+1   Using Progression prog_spec;
Else
	If (aspect_ratio >= 1.0)
		Transfinite Line{1001:1002} = aspect_ratio*2^(mesh_level)*2+1 Using Progression 1;
		Transfinite Line{2001:2002} = 1*2^(mesh_level+1)+1          Using Progression 1.3;
		Transfinite Line{2003}      = 1*2^(mesh_level+1)+1          Using Progression 1.3;
		Transfinite Line{1003:1004} = aspect_ratio*2^(mesh_level)*2+1 Using Progression prog_spec;
	Else
		Transfinite Line{1001:1002} = 2^(mesh_level)*2+1 Using Progression 1;
		Transfinite Line{2001:2002} = 1.0/aspect_ratio*2^(mesh_level+1)+1 Using Progression 1.3;
		Transfinite Line{2003}      = 1.0/aspect_ratio*2^(mesh_level+1)+1 Using Progression 1.3;
		Transfinite Line{1003:1004} = 2^(mesh_level)*2+1 Using Progression prog_spec;
	EndIf
EndIf

Line Loop (4001) = {1001,2003,-1003,-2001};
Line Loop (4002) = {1002,2002,-1004,-2003};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface {4001};
Transfinite Surface {4002};
If (mesh_type == TRI)
	// Do nothing.
ElseIf (mesh_type == QUAD)
	Recombine Surface{4001,4002};
ElseIf (mesh_type == MIXED)
	Recombine Surface{4002};
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
	If (geom_adv == GEOM_ADV_XL)
		Physical Line(bc_base+BC_INFLOW)       = {2001};
		Physical Line(bc_base+BC_OUTFLOW)      = {2002};
		Physical Line(bc_base+BC_OUTFLOW_ALT1) = {1001,1002};
		Physical Line(bc_base+BC_OUTFLOW_ALT2) = {1003,1004};
	ElseIf (geom_adv == GEOM_ADV_YL)
		Physical Line(bc_base+BC_INFLOW)       = {1001,1002};
		Physical Line(bc_base+BC_OUTFLOW)      = {2001};
		Physical Line(bc_base+BC_OUTFLOW_ALT1) = {2002};
		Physical Line(bc_base+BC_OUTFLOW_ALT2) = {1003,1004};
	ElseIf (geom_adv == GEOM_ADV_XYL)
		Physical Line(bc_base+BC_INFLOW)       = {1001:1002};
		Physical Line(bc_base+BC_INFLOW_ALT1)  = {2001};
		Physical Line(bc_base+BC_OUTFLOW)      = {1003:1004};
		Physical Line(bc_base+BC_OUTFLOW_ALT1) = {2002};
	ElseIf (geom_adv == GEOM_ADV_XL_YR)
		Physical Line(bc_base+BC_INFLOW)       = {2001};
		Physical Line(bc_base+BC_INFLOW_ALT1)  = {1003:1004};
		Physical Line(bc_base+BC_OUTFLOW)      = {2002};
		Physical Line(bc_base+BC_OUTFLOW_ALT1) = {1001:1002};
	ElseIf (geom_adv == GEOM_ADV_UPWIND)
		Physical Line(bc_base+BC_UPWIND)      = {2001};
		Physical Line(bc_base+BC_UPWIND_ALT1) = {1003:1004};
		Physical Line(bc_base+BC_UPWIND_ALT2) = {2002};
		Physical Line(bc_base+BC_UPWIND_ALT3) = {1001:1002};
	Else
		Error("Unsupported geom_adv: %d",geom_adv); Exit;
	EndIf
ElseIf (pde_name == DIFFUSION)
	Physical Line(bc_base+BC_DIRICHLET)      = {1001:1002};
	Physical Line(bc_base+BC_DIRICHLET_ALT1) = {2001};
	Physical Line(bc_base+BC_NEUMANN)        = {1003:1004};
	Physical Line(bc_base+BC_NEUMANN_ALT1)   = {2002};
ElseIf (pde_name == EULER)
	If (geom_adv == GEOM_ADV_PERIODIC)
		Physical Line(bc_base+PERIODIC_XL) = {2001};
		Physical Line(bc_base+PERIODIC_XR) = {2002};
		Physical Line(bc_base+PERIODIC_YL) = {1001:1002};
		Physical Line(bc_base+PERIODIC_YR) = {1003:1004};

		// Periodic Indicator (Slave = Master)
		Periodic Line {2002} = {2001}; // Periodic (x)
		Periodic Line {1003} = {1001}; // Periodic (y)
		Periodic Line {1004} = {1002}; // Periodic (y)
	ElseIf (geom_adv == GEOM_ADV_XL)
		Physical Line(bc_base+BC_SUPERSONIC_IN)  = {2001};
		Physical Line(bc_base+BC_SUPERSONIC_OUT) = {2002};
		Physical Line(bc_base+BC_SLIPWALL)       = {1001:1004}; // Can also try with Riemann
	ElseIf (geom_adv == GEOM_ADV_XLR)
		Physical Line(bc_base+BC_RIEMANN)  = {2001,2002};
		Physical Line(bc_base+BC_SLIPWALL) = {1001:1004};
//		Physical Line(bc_straight+BC_RIEMANN) = {1001:1004};
	ElseIf (geom_adv == GEOM_ADV_INTERNAL)
		Physical Line(bc_straight+BC_BACKPRESSURE) = {2002};
		Physical Line(bc_straight+BC_TOTAL_TP)     = {2001};
		Physical Line(bc_base+BC_SLIPWALL)         = {1001:1002};
		Physical Line(bc_straight+BC_SLIPWALL)     = {1003:1004};
	ElseIf (geom_adv == GEOM_NURBS_AIRFOIL_O_GRID)
		// Periodic left and right x face
		//Physical Line(bc_base+BC_BACKPRESSURE) = {2001};
		//Physical Line(bc_base+BC_TOTAL_TP) = {2002};
		
		Physical Line(bc_base+PERIODIC_XL) = {2001};
	 	Physical Line(bc_base+PERIODIC_XR) = {2002};
		Periodic Line {2002} = {2001}; // Periodic (x)

		// Slipwall for bottom face and reimann for top
		Physical Line (bc_base+BC_SLIPWALL)         = {1001:1002};
		Physical Line (bc_base+BC_RIEMANN) = {1003:1004};
	Else
		Error("Unsupported geom_adv: %d",geom_adv); Exit;
	EndIf
ElseIf (pde_name == NAVIER_STOKES)
	Physical Line(bc_base+PERIODIC_XL) = {2001};
	Physical Line(bc_base+PERIODIC_XR) = {2002};

	Periodic Line {2002} = {2001}; // Periodic (x)

	// No slip for remaining boundaries
	Physical Line (bc_base+BC_NOSLIP_T)         = {1001:1002};
	Physical Line (bc_base+BC_NOSLIP_ADIABATIC) = {1003:1004};
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

Physical Surface(9401) = 4001;
Physical Surface(9402) = 4002;



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
