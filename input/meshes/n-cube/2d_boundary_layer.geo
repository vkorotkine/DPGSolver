Include "../parameters.geo";
//mesh_level = 0; mesh_type = MIXED; mesh_domain = PARAMETRIC; pde_name = EULER; geom_conformal = GEOM_CONFORMAL_HALF; pde_spec = STEADY_JOUKOWSKI;

// Geometry Specification
If (pde_name == EULER && pde_spec == STEADY_JOUKOWSKI)
	Include "../../input_files/euler/steady/joukowski/geometry_parameters.geo";
Else
	Error("Unsupported pde_name, pde_spec: %d, %d",pde_name,pde_spec); Exit;
EndIf
Printf("s_offset = %g.",s_offset);

// Geometry Specification
b = s_offset;
l = 1;
h = 1;
bl = 0.1*h;

Point(0) = {-l,b,-0,lc};
Point(1) = {0, b,-0,lc};
Point(2) = {+l,b,-0,lc};
Point(3) = {-l,b+bl,-0,lc};
Point(4) = {0, b+bl,-0,lc};
Point(5) = {+l,b+bl,-0,lc};
Point(6) = {-l,b+h,-0,lc};
Point(7) = {0, b+h,-0,lc};
Point(8) = {+l,b+h,-0,lc};

Line(1000) = {0,1};
Line(1001) = {1,2};
Line(1002) = {3,4};
Line(1003) = {4,5};
Line(1004) = {6,7};
Line(1005) = {7,8};

Line(2000) = {0,3};
Line(2001) = {3,6};
Line(2002) = {1,4};
Line(2003) = {4,7};
Line(2004) = {2,5};
Line(2005) = {5,8};


Line Loop(4000) = {1000,2002,-1002,-2000};
Line Loop(4001) = {1002,2003,-1004,-2001};
Line Loop(4002) = {1001,2004,-1003,-2002};
Line Loop(4003) = {1003,2005,-1005,-2003};

Plane Surface(4000) = {4000};
Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};

all_transfinite = 1;
If (all_transfinite)
	Transfinite Line{2000,2002,2004} = 2^(mesh_level+1)+1 Using Progression 1.2;
	Transfinite Line{2001,2003,2005} = 2^(mesh_level+1)+0 Using Progression 1.4;
	Transfinite Line{1000,1002,1004} = 2^(mesh_level+1)+1 Using Bump 0.2;
	Transfinite Line{1001,1003,1005} = 2^(mesh_level+1)+1 Using Progression 1.4;
	Transfinite Surface {4000:4003};
Else
	Transfinite Line{2000,2002} = 2^(mesh_level+1)+1 Using Progression 1.2;
	Transfinite Line{2001,2003} = 2^(mesh_level+1)+0 Using Progression 1.4;
	Transfinite Line{1000,1002,1004} = 2^(mesh_level+1)+1 Using Bump 0.2;
	Transfinite Line{1001,1003} = 2^(mesh_level+1)+0 Using Progression 1.4;
	Transfinite Surface {4000:4001};
EndIf


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

Physical Line(bc_base+BC_RIEMANN)     = {1004};
Physical Line(bc_straight+BC_RIEMANN) = {1005,2004,2005};
If (pde_name == EULER)
	Physical Line(bc_base+BC_SLIPWALL) = {1000};
	If (geom_conformal == GEOM_CONFORMAL_HALF)
		Physical Line(bc_straight+BC_SLIPWALL) = {2000,2001,1001};
	ElseIf (geom_conformal == GEOM_CONFORMAL_FULL)
        Error("Put this option in a different geo file"); Exit;
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
		Physical Line(bc_straight+BC_SLIPWALL) = {2000,2001,1001};
	ElseIf (geom_conformal == GEOM_CONFORMAL_FULL)
		Error("Put this option in a different geo file"); Exit;
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
Physical Surface(9402) = 4002;
Physical Surface(9403) = 4003;



// Visualization in gmsh

Color Black{ Surface{4000:4003}; }
Geometry.Color.Points = Black;
