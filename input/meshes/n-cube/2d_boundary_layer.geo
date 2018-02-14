Include "../parameters.geo";
//mesh_level = 0; mesh_type = MIXED; mesh_domain = PARAMETRIC; pde_name = EULER; geom_conformal = GEOM_CONFORMAL_HALF;

// Geometry Specification
l = 1;
h = 1;
bl = 0.1*h;

Point(0) = {-l,-h,-0,lc};
Point(1) = {+l,-h,-0,lc};
Point(2) = {-l,+h,-0,lc};
Point(3) = {+l,+h,-0,lc};
Point(4) = {-l,-h+bl,-0,lc};
Point(5) = {+l,-h+bl,-0,lc};

Line(1000) = {0,1};
Line(1001) = {2,3};
Line(1002) = {4,5};
Line(2000) = {0,4};
Line(2001) = {4,2};
Line(2002) = {1,5};
Line(2003) = {5,3};

Transfinite Line{2000,2002} = 2^(mesh_level+2)+1 Using Progression 1.2;
Transfinite Line{1000,1002} = 2^(mesh_level+2)+1 Using Bump 0.2;
Transfinite Line{2001,2003} = 2^(mesh_level+2)+0 Using Progression 1.4;

Line Loop (4000) = {1000,2002,-1002,-2000};
Line Loop (4001) = {1002,2003,-1001,-2001};

Plane Surface(4000) = {4000};
Plane Surface(4001) = {4001};

Transfinite Surface {4000};
If (mesh_type == MIXED)
	Recombine Surface{4000};
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

Physical Line(bc_base+BC_RIEMANN) = {1001};
If (pde_name == EULER)
	Physical Line(bc_base+BC_SLIPWALL)     = {1000};
	If (geom_conformal == GEOM_CONFORMAL_HALF)
		Physical Line(bc_straight+BC_SLIPWALL) = {2000:2003};
	ElseIf (geom_conformal == GEOM_CONFORMAL_FULL)
		Physical Line(bc_straight+PERIODIC_XL) = {2000:2001};
		Physical Line(bc_straight+PERIODIC_XR) = {2002:2003};

		// Periodic Indicator (Slave = Master)
		Periodic Line {2000} = {2002}; // Periodic (x)
		Periodic Line {2001} = {2003};
	Else
		Error("Unsupported geom_conformal: %d",geom_conformal); Exit;
	EndIf
ElseIf (pde_name == NAVIER_STOKES)
	Physical Line(bc_base+BC_NOSLIP_ADIABATIC) = {1000};
	If (geom_conformal == GEOM_CONFORMAL_HALF)
		Physical Line(bc_straight+BC_SLIPWALL) = {2000:2003};
	ElseIf (geom_conformal == GEOM_CONFORMAL_FULL)
		Physical Line(bc_straight+PERIODIC_XL) = {2000:2001};
		Physical Line(bc_straight+PERIODIC_XR) = {2002:2003};

		// Periodic Indicator (Slave = Master)
		Periodic Line {2000} = {2002}; // Periodic (x)
		Periodic Line {2001} = {2003};
	Else
		Error("Unsupported geom_conformal: %d",geom_conformal); Exit;
	EndIf
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

Physical Surface(9400) = 4000;
Physical Surface(9401) = 4001;



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
